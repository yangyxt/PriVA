#!/usr/bin/env python3
"""
gpt_mavedb_assay_mapper.py

Iterate over MaveDB experiment metadata JSON files, prompt GPT to map each
experiment to a standardized assay spec (family/subtype/directionality/scale/units)
and select the corresponding interpretation function. Writes a TSV summary.

Usage:
  python gpt_mavedb_assay_mapper.py \
      --assay-tsv assay_catalog.tsv \
      --metadata-json meta_formatted.json \
      --prompt-template assay_mapper_prompt.txt \
      --out mapped_specs.tsv \
      --model gpt-5  # override as needed
"""

import os
import sys
import json
import time
import argparse
import logging
import pandas as pd
from typing import List, Dict, Any, Optional


KEEP_KEYS = {
    "urn",
    "title",
    "shortDescription",
    "methodText",
    "abstractText",
    "numVariants",
    "datasetColumns",
    "scoreRanges",
    "targetGenes",
    "keywords",
    "processingState",
    "mappingErrors",
}


logger = logging.getLogger()
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(processName)s:%(funcName)s:%(lineno)s:%(name)s:%(exc_info)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)
logger.setLevel(os.environ.get("LOGLEVEL", "INFO"))

try:
    from openai import OpenAI
    _HAS_OPENAI = True
except Exception:
    _HAS_OPENAI = False
    exit("Please install the OpenAI Python client: pip install openai")

def read_text(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        return f.read()

def read_json(path: str) -> Any:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def build_prompt(template: str, meta_obj: Dict[str, Any]) -> str:
    return (template
            .replace("{{METADATA_JSON}}", json.dumps(meta_obj, ensure_ascii=False, indent=2)))


def chat_with_openai(messages):
    client = OpenAI(
        api_key=os.environ.get("OPENAI_API_KEY", "sk-8WrRbv0xCYM7hPs_N4vmHjDLgYoB-TkuD6wZpLSciEoi90ekiXIzUx-eSsQ"),
        base_url=os.environ.get("OPENAI_BASE_URL", "https://api.zmon.me/v1")
    )

    try:
        response = client.chat.completions.create(
            model="gpt-5",
            messages=[{"role": "user", "content": messages}],
            temperature=0.0,
            max_tokens=5000
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Error: {e}"
    

def call_gpt(prompt: str, model: str) -> str:
    """
    Calls the Responses API via the OpenAI Python client.
    """
    client = OpenAI(
        api_key=os.environ.get("OPENAI_API_KEY", "sk-8WrRbv0xCYM7hPs_N4vmHjDLgYoB-TkuD6wZpLSciEoi90ekiXIzUx-eSsQ"),
        base_url=os.environ.get("OPENAI_BASE_URL", "https://api.zmon.me/v1")
    )
    
    resp = client.responses.create(
        model=model,
        input=prompt,
        temperature=0.0,
        max_output_tokens=5000,
        reasoning={"effort": "high"}
    )

    # Prefer the SDK’s convenience property if available
    if hasattr(resp, "output_text") and resp.output_text:
        logger.info(f"Received output text from Responses API, which looks like:\n{resp.output_text}\n")
        return resp.output_text
    
    logger.warning(f"Falling back to manual extraction from response object, what deos the response look like? \n{resp}\n")

    # Robust extraction via dict dump (avoids typing/union issues)
    try:
        data = resp.model_dump()
    except Exception:
        data = getattr(resp, "to_dict", lambda: None)() or getattr(resp, "__dict__", {}) or {}

    # 1) Responses API shape
    parts = []
    for item in data.get("output", []) or []:
        for c in item.get("content", []) or []:
            if c.get("type") in ("output_text", "text"):
                parts.append(c.get("text", ""))
    if parts:
        return "".join(parts)

    # 2) Chat Completions-like fallback (if proxy returns choices)
    if "choices" in data and data["choices"]:
        msg = data["choices"][0].get("message", {})
        if "content" in msg and msg["content"]:
            return msg["content"]

    # Last resort
    return str(data)


def validate_and_row(urn: str, obj: Dict[str, Any]) -> Dict[str, Any]:
    def get(o, k, default="NA"):
        v = o.get(k, default)
        if v in (None, "", []):
            return default
        return v
    
    flags = obj.get("flags", {}) if isinstance(obj.get("flags"), dict) else {}
    row = {
        "urn": urn,
        "assay_family": get(obj, "assay_family"),
        "subtype_enum": get(obj, "subtype_enum"),
        "subtype_value": get(obj, "subtype_value"),
        "directionality": get(obj, "directionality"),
        "scale": get(obj, "scale"),
        "baseline": get(obj, "baseline"),
        "units": get(obj, "units"),
        "interpreter": get(obj, "interpreter"),
        "neutral_bands": get(obj, "neutral_bands", "[NA, NA]"),
        "flexibility_requirement": bool(get(obj, "flexibility_requirement", False)),
        "needs_gene_specific_rules": bool(flags.get("needs_gene_specific_rules", False)),
        "nonmonotonic": bool(flags.get("nonmonotonic", False)),
        "ambiguous": bool(flags.get("ambiguous", False)),
        "rationale": get(obj, "rationale", ""),
    }
    
    if row["assay_family"] == "Other":
        row["interpreter"] = "NA"
    return row

def append_tsv_row(row: Dict[str, Any], out_path: str, write_header: bool = True) -> None:
    """Append a single row to TSV file in real-time"""
    cols = ["urn","assay_family","subtype_enum","subtype_value","directionality",
            "scale","baseline","units","interpreter","neutral_bands",
            "flexibility_requirement",
            "needs_gene_specific_rules","nonmonotonic","ambiguous","rationale"]
    
    mode = 'w' if write_header else 'a'
    with open(out_path, mode, encoding="utf-8") as f:
        if write_header:
            f.write("\t".join(cols) + "\n")
        
        vals = []
        for c in cols:
            v = row.get(c, "")
            if isinstance(v, bool):
                v = "true" if v else "false"
            elif isinstance(v, float):
                v = f"{v:.3f}"
            else:
                v = str(v)
            v = v.replace("\t", " ").replace("\n", " ").strip()
            vals.append(v)
        f.write("\t".join(vals) + "\n")


def load_processed_urns(tsv_path: str) -> set:
    """
    Return a set of URNs from the first column of an existing TSV file.
    Handles both headered (with 'urn' column) and headerless files.
    """
    processed: set = set()
    if not tsv_path or not os.path.exists(tsv_path) or os.path.getsize(tsv_path) == 0:
        return processed

    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)
        if "urn" in df.columns:
            col = df["urn"].dropna().astype(str)
        else:
            # No header; read first column explicitly
            df = pd.read_csv(tsv_path, sep="\t", header=None, dtype=str, usecols=[0])
            col = df.iloc[:, 0].dropna().astype(str)
        vals = [v.strip() for v in col.tolist() if v and v.strip().lower() != "urn"]
        return set(vals)
    except Exception as e:
        logger.warning(f"Pandas failed to read existing TSV {tsv_path}: {e}. Falling back to plain parsing.")

    return processed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata-json", required=True, help="JSON file with experimentSets")
    ap.add_argument("--prompt-template", required=True, help="Prompt template file")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--model", default=os.environ.get("GPT_MODEL", "gpt-5"))
    ap.add_argument("--sleep", type=float, default=0.0, help="Sleep seconds between requests")
    args = ap.parse_args()

    template = read_text(args.prompt_template)
    
    # Load the MaveDB JSON file
    mavedb_data = read_json(args.metadata_json)
    
    # Extract all experiments from experimentSets
    all_experiment_scoresets = []
    for exp_set in mavedb_data.get("experimentSets", []):
        for experiment in exp_set.get("experiments", []):
            for scoresets in experiment.get("scoreSets", []):
                all_experiment_scoresets.append({k: scoresets.get(k) for k in KEEP_KEYS})
    
    logger.info(f"Found {len(all_experiment_scoresets)} experiment scoresets to process")
    processed_urns = load_processed_urns(args.out)
    if processed_urns:
        logger.info(f"Found {len(processed_urns)} already processed URNs in {args.out}, skipping them")
        all_experiment_scoresets = [s for s in all_experiment_scoresets if s.get("urn") not in processed_urns]
        logger.info(f"{len(all_experiment_scoresets)} experiment scoresets remain to process")
    
    out_rows: List[Dict[str, Any]] = []
    for i, scoreset in enumerate(all_experiment_scoresets, len(processed_urns) + 1):
        urn = scoreset.get("urn", f"unknown_{i}")
        logger.info(f"Processing {i}/{len(all_experiment_scoresets)}: {urn}, the metadata json block looks like:\n{scoreset}\n")
        
        prompt = build_prompt(template, scoreset)

        try:
            text = chat_with_openai(messages=prompt)
            # Text is already json format, convert to dict object
            logger.info(f"The returning text is: \n{text}\n")
            obj = json.loads(text)
            assert isinstance(obj, dict), f"Output is not a JSON object: \n{text}\n"
        except Exception as e:
            logger.warning(f"Error processing {urn}: {e}, the input json block is \n{scoreset}\n")
            raise e
            

        row = validate_and_row(urn, obj if isinstance(obj, dict) else {})
        # Write row immediately after processing 
        if i == 1: 
            append_tsv_row(row, args.out, write_header=True)
        else:
            append_tsv_row(row, args.out, write_header=False)
        logger.info(f"✓ Wrote {urn} to {args.out}\n")

        # if i > 1: break # DEBUG limit to first 3 experiments
        
        if args.sleep > 0:
            time.sleep(args.sleep)

    logger.info(f"Wrote {i} rows to {args.out}")

if __name__ == "__main__":
    main()