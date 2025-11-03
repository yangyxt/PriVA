# -*- coding: utf-8 -*-
"""
mavedb_interpretation.py

A lightweight, hierarchical pool of interpretation functions for MaveDB-style
assay outputs. This module avoids per-dataset calibration and provides
conservative, widely-used heuristics for *standardized* score spaces.
If an experiments semantics are unclear, functions return an NA call.

Core idea:
- LLM API assign each experiment a categorical (assay_family, subtype, directionality, scale, baseline, units,
  neutral_band) record and save it in a metadata table keyed by mavedb_id.
- Pass the per-experiment metadata + per-variant score into this module.
- The module standardizes the score to a centered "effect space" and emits:
    call ∈ {"LoF","GoF","Neutral","Unknown","NA"}
    standardized_score (float | None)
    reason (str), notes (str)

Design constraints:
- Heuristics only for broadly standardized assays; otherwise NA.
- Optional gene-level overrides (e.g., function_prefers_lower_stability).

This file intentionally DOES NOT convert calls into PS3/BS3 strengths.
Use these calls as one ingredient in downstream logic; if a dataset needs
tailored thresholds from a paper, skip it (NA).

Author: generated for PriVA
"""

from dataclasses import dataclass
from typing import Optional, Tuple, Dict, Any, Type
import math
import numpy as np
from enum import Enum


# ----------------------------- Families (string constants) -----------------------------

class AssayFamily:
    STABILITY_K50 = "Stability_K50"
    STABILITY_DG  = "Stability_DG"
    ABUNDANCE     = "Abundance"
    ACTIVITY      = "Activity"
    BINDING       = "Binding"
    FITNESS       = "Fitness"
    SPLICING      = "Splicing"
    EPHYS         = "Electrophysiology"
    COMPLEMENT    = "Complementation"
    OTHER         = "Other"


# ----------------------------- Categorical Subtypes (Enums) -----------------------------

class StabilityK50Subtype(Enum):
    LOG10_K50_STANDARDIZED = "LOG10_K50_STANDARDIZED"   # centered/log transformed, WT~=0
    LOG10_K50_RELATIVE     = "LOG10_K50_RELATIVE"       # WT~=1:log2; centered to WT~=0
    K50_RAW                = "K50_RAW"                  # raw K50:not interpreted here


class StabilityDGSubtype(Enum):
    DDG_KCAL_MOL = "DDG_KCAL_MOL"    # ΔΔG (kcal/mol), positive = destabilizing
    DG_KCAL_MOL  = "DG_KCAL_MOL"     # ΔG (kcal/mol), more negative = more stable


class AbundanceSubtype(Enum):
    VAMP_SEQ             = "VAMP_SEQ"
    CELL_SURFACE_TRAFF   = "CELL_SURFACE_TRAFF"
    TOTAL_EXPRESSION     = "TOTAL_EXPRESSION"


class ActivitySubtype(Enum):
    ENZYME_ACTIVITY_REL   = "ENZYME_ACTIVITY_REL"       # relative to WT=1
    REPORTER_ACTIVITY_REL = "REPORTER_ACTIVITY_REL"     # relative to WT=1


class BindingSubtype(Enum):
    KD          = "KD"             # lower = stronger (converted to ~pKd)
    PKD         = "PKD"            # higher = stronger
    KA          = "KA"             # higher = stronger (treated like pKd)
    KD_APPARENT = "KD_APPARENT"    # treat like KD if units valid


class FitnessSubtype(Enum):
    LOG2_ENRICHMENT = "LOG2_ENRICHMENT"


class SplicingSubtype(Enum):
    DELTA_PSI_PERCENT  = "DELTA_PSI_PERCENT"    # -100..100 or -50..50
    DELTA_PSI_FRACTION = "DELTA_PSI_FRACTION"   # -1..1


class EphysSubtype(Enum):
    V1_2            = "V1_2"
    CURRENT_DENSITY = "CURRENT_DENSITY"


class ComplementationSubtype(Enum):
    YEAST          = "YEAST"
    MAMMALIAN_CELL = "MAMMALIAN_CELL"


# ----------------------------- Other enums -----------------------------

class Scale(Enum):
    CENTERED_0 = "CENTERED_0"   # WT~=0 (e.g., log2 or z-score)
    WT_EQ_1    = "WT_EQ_1"      # WT~=1:log2 transform to center
    RAW        = "RAW"          # raw units (kcal/mol, Kd)
    PERCENT    = "PERCENT"      # e.g., ΔPSI in percentage points
    FRACTION   = "FRACTION"     # e.g., ΔPSI as fraction


class Directionality(Enum):
    HIGHER_IS_FUNCTION = "HIGHER_IS_FUNCTION"
    LOWER_IS_FUNCTION  = "LOWER_IS_FUNCTION"
    NON_MONOTONIC      = "NON_MONOTONIC"


# family:required subtype enum (for validation)
FAMILY_TO_SUBTYPE_ENUM: Dict[str, Type[Enum]] = {
    AssayFamily.STABILITY_K50: StabilityK50Subtype,
    AssayFamily.STABILITY_DG:  StabilityDGSubtype,
    AssayFamily.ABUNDANCE:     AbundanceSubtype,
    AssayFamily.ACTIVITY:      ActivitySubtype,
    AssayFamily.BINDING:       BindingSubtype,
    AssayFamily.FITNESS:       FitnessSubtype,
    AssayFamily.SPLICING:      SplicingSubtype,
    AssayFamily.EPHYS:         EphysSubtype,
    AssayFamily.COMPLEMENT:    ComplementationSubtype,
}


@dataclass
class AssaySpec:
    assay_family: str
    subtype: Enum                      # categorical subtype
    directionality: Directionality
    scale: Scale
    baseline: str                      # short human description
    units: str
    neutral_band: Tuple[float, float]  # inclusive (lo, hi)
    flexibility_requirement: bool = False  # If True, function prefers lower stability


# ----------------------------- Utilities -----------------------------

def _to_log2(x: float) -> Optional[float]:
    try:
        if x is None or x <= 0:
            return None
        return math.log(x, 2)
    except Exception:
        return None


def validate_spec(spec: AssaySpec) -> bool:
    enum_cls = FAMILY_TO_SUBTYPE_ENUM.get(spec.assay_family)
    if not enum_cls or not isinstance(spec.subtype, enum_cls):
        return False
    lo, hi = spec.neutral_band
    return lo < hi


# ----------------------------- Base Interpreter -----------------------------

class BaseInterpreter:
    def standardize(self, score: float, spec: AssaySpec) -> Optional[float]:
        if score is None:
            return None
        if spec.scale == Scale.CENTERED_0:
            return float(score)
        if spec.scale == Scale.WT_EQ_1:
            return _to_log2(float(score))
        if spec.scale in (Scale.PERCENT, Scale.FRACTION):
            return float(score)
        # RAW handled by specific interpreters where applicable
        return None

    def call_from_band(self, s: float, spec: AssaySpec) -> str:
        lo, hi = spec.neutral_band
        if s is None:
            return "Unknown"
        if lo <= s <= hi:
            return "Neutral"

        if spec.directionality == Directionality.HIGHER_IS_FUNCTION:
            return "GoF" if s > hi else "LoF"
        elif spec.directionality == Directionality.LOWER_IS_FUNCTION:
            # Use flexibility_requirement from spec
            if spec.flexibility_requirement:
                # Function prefers lower stability: lower values are GoF
                return "GoF" if s < lo else "LoF"
            else:
                # Standard interpretation: lower values are LoF
                return "LoF" if s < lo else "GoF"
        else:
            return "Unknown"

    def interpret(self, score: float, spec: AssaySpec, **kwargs) -> Dict[str, Any]:
        raise NotImplementedError


# ----------------------------- Interpreters -----------------------------

class StabilityK50Interpreter(BaseInterpreter):
    def interpret(self, score: float, spec: AssaySpec, **kwargs):
        if not validate_spec(spec):
            return {"call":"NA","standardized_score":None,"reason":"Invalid spec","notes":""}

        if spec.subtype in (StabilityK50Subtype.LOG10_K50_STANDARDIZED, StabilityK50Subtype.LOG10_K50_RELATIVE):
            s = self.standardize(score, spec)
            if s is None:
                return {"call":"NA","standardized_score":None,"reason":"K50 not standardized","notes":spec.units}
            call = self.call_from_band(s, spec)
            flex_note = " | flexibility_req=True" if spec.flexibility_requirement else ""
            return {"call":call,"standardized_score":s,"reason":"K50 heuristic band","notes":f"{spec.neutral_band} | {spec.units}{flex_note}"}

        return {"call":"NA","standardized_score":None,"reason":"K50 RAW requires transform","notes":spec.subtype.value}


class StabilityDGInterpreter(BaseInterpreter):
    def interpret(self, score: float, spec: AssaySpec, **kwargs):
        if not validate_spec(spec):
            return {"call":"NA","standardized_score":None,"reason":"Invalid spec","notes":""}

        s = self.standardize(score, spec)
        raw_used = False
        if s is None:
            raw_used = True
            if spec.subtype == StabilityDGSubtype.DDG_KCAL_MOL:
                s = -float(score)           # ddG: higher (+) = less stable:invert so higher = more function
            elif spec.subtype == StabilityDGSubtype.DG_KCAL_MOL:
                s = float(score)            # dG: more negative = more stable (already lower→more function)
            else:
                return {"call":"NA","standardized_score":None,"reason":"Unknown ΔG subtype for RAW","notes":spec.units}

        call = self.call_from_band(s, spec)
        flex_note = " | flexibility_req=True" if spec.flexibility_requirement else ""
        return {"call":call,"standardized_score":s,"reason":"ΔG/ΔΔG transform" if raw_used else "ΔG/ΔΔG standardized","notes":f"{spec.neutral_band} | {spec.units}{flex_note}"}


class AbundanceInterpreter(BaseInterpreter):
    def interpret(self, score: float, spec: AssaySpec, **kwargs):
        if not validate_spec(spec):
            return {"call":"NA","standardized_score":None,"reason":"Invalid spec","notes":""}
        s = self.standardize(score, spec)
        if s is None:
            return {"call":"NA","standardized_score":None,"reason":"Abundance not standardized","notes":spec.units}
        return {"call":self.call_from_band(s, spec),"standardized_score":s,"reason":"Abundance band","notes":f"{spec.neutral_band} | {spec.subtype.value}"}


class ActivityInterpreter(BaseInterpreter):
    def interpret(self, score: float, spec: AssaySpec, **kwargs):
        if not validate_spec(spec):
            return {"call":"NA","standardized_score":None,"reason":"Invalid spec","notes":""}
        s = self.standardize(score, spec)
        if s is None:
            return {"call":"NA","standardized_score":None,"reason":"Activity not standardized","notes":spec.units}
        return {"call":self.call_from_band(s, spec),"standardized_score":s,"reason":"Activity band","notes":f"{spec.neutral_band} | {spec.subtype.value}"}


class BindingInterpreter(BaseInterpreter):
    def interpret(self, score: float, spec: AssaySpec, **kwargs):
        if not validate_spec(spec):
            return {"call":"NA","standardized_score":None,"reason":"Invalid spec","notes":""}

        s = self.standardize(score, spec)
        if s is None:
            # RAW: convert KD/KD_APPARENT:~pKd; PKD/KA already "higher=stronger"
            if spec.subtype in (BindingSubtype.KD, BindingSubtype.KD_APPARENT):
                kd = float(score)
                if kd <= 0:
                    return {"call":"NA","standardized_score":None,"reason":"Kd <= 0 invalid","notes":""}
                s = math.log(1.0/kd, 10)  # ~ -log10(Kd)
                if spec.directionality == Directionality.NON_MONOTONIC:
                    spec.directionality = Directionality.HIGHER_IS_FUNCTION
            elif spec.subtype in (BindingSubtype.PKD, BindingSubtype.KA):
                s = float(score)
                if spec.directionality == Directionality.NON_MONOTONIC:
                    spec.directionality = Directionality.HIGHER_IS_FUNCTION
            else:
                return {"call":"NA","standardized_score":None,"reason":"Unsupported Binding subtype for RAW","notes":spec.subtype.value}

        if spec.directionality == Directionality.NON_MONOTONIC:
            return {"call":"NA","standardized_score":s,"reason":"Binding directionality unknown","notes":spec.units}

        return {"call":self.call_from_band(s, spec),"standardized_score":s,"reason":"Binding band","notes":f"{spec.neutral_band} | {spec.subtype.value}"}


class FitnessInterpreter(BaseInterpreter):
    def interpret(self, score: float, spec: AssaySpec, **kwargs):
        if not validate_spec(spec):
            return {"call":"NA","standardized_score":None,"reason":"Invalid spec","notes":""}
        s = self.standardize(score, spec)
        if s is None:
            return {"call":"NA","standardized_score":None,"reason":"Fitness not standardized","notes":spec.units}
        return {"call":self.call_from_band(s, spec),"standardized_score":s,"reason":"Fitness band","notes":f"{spec.neutral_band}"}


class SplicingInterpreter(BaseInterpreter):
    def interpret(self, score: float, spec: AssaySpec, **kwargs):
        if not validate_spec(spec):
            return {"call":"NA","standardized_score":None,"reason":"Invalid spec","notes":""}
        s = self.standardize(score, spec)
        if s is None:
            return {"call":"NA","standardized_score":None,"reason":"Splicing not standardized","notes":spec.units}
        return {"call":self.call_from_band(s, spec),"standardized_score":s,"reason":"Splicing band","notes":f"{spec.neutral_band} | ΔPSI mapping is gene/context-specific"}


class EphysInterpreter(BaseInterpreter):
    def interpret(self, score: float, spec: AssaySpec, gene_specific_rules: Optional[Dict[str, Any]]=None, **kwargs):
        if not validate_spec(spec):
            return {"call":"NA","standardized_score":None,"reason":"Invalid spec","notes":""}
        if not gene_specific_rules:
            return {"call":"NA","standardized_score":None,"reason":"Gene-specific rule required","notes":""}
        s = self.standardize(score, spec)
        if s is None:
            return {"call":"NA","standardized_score":None,"reason":"Ephys not standardized","notes":spec.units}

        # Map V1/2 direction to GoF/LoF using gene rule
        gof_dir = gene_specific_rules.get("gof_dir", "")
        lo, hi = spec.neutral_band
        if lo <= s <= hi:
            return {"call":"Neutral","standardized_score":s,"reason":"Ephys neutral band","notes":spec.subtype.value}

        if gof_dir == "more_negative":
            call = "GoF" if s < lo else "LoF"
        elif gof_dir == "more_positive":
            call = "GoF" if s > hi else "LoF"
        else:
            call = "Unknown"

        return {"call":call,"standardized_score":s,"reason":"Ephys rule","notes":f"{spec.neutral_band} | {spec.subtype.value} | {gof_dir}"}


class ComplementationInterpreter(BaseInterpreter):
    def interpret(self, score: float, spec: AssaySpec, **kwargs):
        if not validate_spec(spec):
            return {"call":"NA","standardized_score":None,"reason":"Invalid spec","notes":""}
        s = self.standardize(score, spec)
        if s is None:
            return {"call":"NA","standardized_score":None,"reason":"Complementation not standardized","notes":spec.units}
        return {"call":self.call_from_band(s, spec),"standardized_score":s,"reason":"Complementation band","notes":f"{spec.neutral_band} | {spec.subtype.value}"}


# ----------------------------- Registry & API -----------------------------

REGISTRY = {
    AssayFamily.STABILITY_K50: StabilityK50Interpreter(),
    AssayFamily.STABILITY_DG:  StabilityDGInterpreter(),
    AssayFamily.ABUNDANCE:     AbundanceInterpreter(),
    AssayFamily.ACTIVITY:      ActivityInterpreter(),
    AssayFamily.BINDING:       BindingInterpreter(),
    AssayFamily.FITNESS:       FitnessInterpreter(),
    AssayFamily.SPLICING:      SplicingInterpreter(),
    AssayFamily.EPHYS:         EphysInterpreter(),
    AssayFamily.COMPLEMENT:    ComplementationInterpreter(),
}


def interpret_score(score: float,
                    spec: AssaySpec,
                    gene_specific_rules: Optional[Dict[str, Any]]=None) -> Dict[str, Any]:
    """
    Interpret a functional score using the appropriate assay family interpreter.
    
    Args:
        score: The variant functional score
        spec: AssaySpec containing all interpretation parameters (including flexibility_requirement)
        gene_specific_rules: Optional gene-specific rules (only for Ephys)
    
    Returns:
        Dict with keys: call, standardized_score, reason, notes
    """
    interpreter = REGISTRY.get(spec.assay_family)
    if not interpreter:
        return {"call":"NA","standardized_score":None,"reason":"No interpreter","notes":spec.assay_family}

    kwargs = {}
    if spec.assay_family == AssayFamily.EPHYS:
        kwargs["gene_specific_rules"] = gene_specific_rules or {}

    return interpreter.interpret(score, spec, **kwargs)


# ----------------------------- Defaults (using categorical subtypes) -----------------------------

DEFAULT_SPECS = {
    "Stability_K50__default": AssaySpec(
        assay_family=AssayFamily.STABILITY_K50,
        subtype=StabilityK50Subtype.LOG10_K50_STANDARDIZED,
        directionality=Directionality.HIGHER_IS_FUNCTION,
        scale=Scale.CENTERED_0,
        baseline="WT~=0 (after log2 from relative)",
        units="log2_from_relative",
        neutral_band=(-0.5, 0.5),
    ),
    "Stability_DG__ddg_default": AssaySpec(
        assay_family=AssayFamily.STABILITY_DG,
        subtype=StabilityDGSubtype.DDG_KCAL_MOL,
        directionality=Directionality.LOWER_IS_FUNCTION,
        scale=Scale.CENTERED_0,  # RAW allowed with minimal transform
        baseline="WT delta~=0",
        units="kcal/mol",
        neutral_band=(-0.5, 0.5),
    ),
    "Stability_DG__dg_default": AssaySpec(
        assay_family=AssayFamily.STABILITY_DG,
        subtype=StabilityDGSubtype.DG_KCAL_MOL,
        directionality=Directionality.LOWER_IS_FUNCTION,
        scale=Scale.CENTERED_0,
        baseline="WT delta~=0",
        units="kcal/mol",
        neutral_band=(-0.5, 0.5),
    ),
    "Abundance__VAMP_seq": AssaySpec(
        assay_family=AssayFamily.ABUNDANCE,
        subtype=AbundanceSubtype.VAMP_SEQ,
        directionality=Directionality.HIGHER_IS_FUNCTION,
        scale=Scale.WT_EQ_1,
        baseline="WT=1:log2",
        units="relative",
        neutral_band=(-0.32, 0.32),
    ),
    "Abundance__CellSurface": AssaySpec(
        assay_family=AssayFamily.ABUNDANCE,
        subtype=AbundanceSubtype.CELL_SURFACE_TRAFF,
        directionality=Directionality.HIGHER_IS_FUNCTION,
        scale=Scale.WT_EQ_1,
        baseline="WT=1:log2",
        units="relative",
        neutral_band=(-0.32, 0.32),
    ),
    "Activity__default": AssaySpec(
        assay_family=AssayFamily.ACTIVITY,
        subtype=ActivitySubtype.ENZYME_ACTIVITY_REL,
        directionality=Directionality.HIGHER_IS_FUNCTION,
        scale=Scale.WT_EQ_1,
        baseline="WT=1:log2",
        units="relative",
        neutral_band=(-0.32, 0.32),
    ),
    "Binding__pKd_default": AssaySpec(
        assay_family=AssayFamily.BINDING,
        subtype=BindingSubtype.PKD,
        directionality=Directionality.HIGHER_IS_FUNCTION,
        scale=Scale.CENTERED_0,
        baseline="WT delta~=0",
        units="pKd",
        neutral_band=(-0.3, 0.3),
    ),
    "Fitness__default": AssaySpec(
        assay_family=AssayFamily.FITNESS,
        subtype=FitnessSubtype.LOG2_ENRICHMENT,
        directionality=Directionality.HIGHER_IS_FUNCTION,
        scale=Scale.CENTERED_0,
        baseline="WT~=0",
        units="log2",
        neutral_band=(-0.5, 0.5),
    ),
    "Splicing__deltaPSI_percent": AssaySpec(
        assay_family=AssayFamily.SPLICING,
        subtype=SplicingSubtype.DELTA_PSI_PERCENT,
        directionality=Directionality.HIGHER_IS_FUNCTION,
        scale=Scale.CENTERED_0,
        baseline="0~=WT inclusion",
        units="ΔPSI%",
        neutral_band=(-10.0, 10.0),
    ),
    "Ephys__default": AssaySpec(
        assay_family=AssayFamily.EPHYS,
        subtype=EphysSubtype.V1_2,
        directionality=Directionality.NON_MONOTONIC,
        scale=Scale.CENTERED_0,
        baseline="WT~=0",
        units="mV",
        neutral_band=(-2.0, 2.0),
    ),
    "Complementation__default": AssaySpec(
        assay_family=AssayFamily.COMPLEMENT,
        subtype=ComplementationSubtype.YEAST,
        directionality=Directionality.HIGHER_IS_FUNCTION,
        scale=Scale.WT_EQ_1,
        baseline="WT=1:log2",
        units="relative",
        neutral_band=(-0.32, 0.32),
    ),
}


def parse_neutral_band_from_metadata(meta: Dict[str, Any]) -> Tuple[Optional[float], Optional[float]]:
    """
    Parse neutral band from metadata with the following priority:
    1. neutral_bands column (may have one-sided intervals like "[0, NA]" or "[NA, 0.5]")
    2. q05/q95 columns (percentile-based bounds)
    3. Return (None, None) if unparseable
    
    Args:
        meta: Dictionary with keys: neutral_bands, q05, q95, flexibility_requirement
    
    Returns:
        Tuple of (lo, hi) where either or both can be None for one-sided intervals
    """    
    # Try parsing neutral_bands column first
    neutral_band_str = str(meta.get('neutral_bands', '[NA, NA]')).strip(" []") 
    if neutral_band_str not in ('[NA, NA]', 'NA', '', 'None', 'nan'):
        try:
            # Parse "[lo, hi]" format, handling NA values
            parsed = [ x.strip() for x in neutral_band_str.split(',') ]
            parsed = [np.nan if x in ('NA', 'nan', 'None', '') else float(x) for x in parsed]
            if isinstance(parsed, (list, tuple)) and len(parsed) == 2:
                lo, hi = parsed
                
                # Convert "NA" strings or None to None
                if lo in ('NA', 'nan', None) or (isinstance(lo, float) and math.isnan(lo)):
                    lo = None
                else:
                    lo = float(lo)
                
                if hi in ('NA', 'nan', None) or (isinstance(hi, float) and math.isnan(hi)):
                    hi = None
                else:
                    hi = float(hi)
                
                # One-sided interval interpretation with flexibility_requirement
                flexibility_req = meta.get('flexibility_requirement', False)
                
                if lo is None and hi is not None:
                    # Pattern: [NA, 0.5] or [NA, 0]
                    if flexibility_req:
                        # All values <= hi are neutral/flexible; above hi is LoF
                        return (-float('inf'), hi)
                    else:
                        # All values <= hi are LoF; above hi is neutral
                        return (hi, float('inf'))
                
                elif hi is None and lo is not None:
                    # Pattern: [0, NA] or [0.5, NA]
                    if flexibility_req:
                        # All values >= lo are neutral/flexible; below lo is LoF
                        return (lo, float('inf'))
                    else:
                        # All values >= lo are LoF; below lo is neutral
                        return (-float('inf'), lo)
                
                elif lo is not None and hi is not None:
                    # Two-sided: [lo, hi]
                    if lo < hi:
                        return (lo, hi)
                    else:
                        # Invalid interval, fall through to q05/q95
                        pass
                
        except (ValueError, SyntaxError, TypeError):
            pass  # Fall through to q05/q95
    
    # Fallback to q05/q95 percentiles
    try:
        q05 = meta.get('q05')
        q95 = meta.get('q95')
        
        if q05 is not None and q95 is not None:
            q05 = float(q05)
            q95 = float(q95)
            if not (math.isnan(q05) or math.isnan(q95)) and q05 < q95:
                return (q05, q95)
    except (ValueError, TypeError):
        pass
    
    # Unable to determine neutral band
    return (None, None)


def pick_interpreter_from_metadata(
    urn: str,
    score: float,
    metadata_df,  # pandas DataFrame or dict lookup
) -> Dict[str, Any]:
    """
    Pick the correct interpretation function and arguments based on metadata table.
    
    Args:
        urn: MaveDB URN (e.g., "urn:mavedb:00000001-a-1")
        score: The variant score from VEP annotation
        metadata_df: DataFrame or dict with columns matching your metadata table
        function_prefers_lower_stability: For stability assays where lower may be GoF
    
    Returns:
        Dict with keys: call, standardized_score, reason, notes, metadata_flags, urn
    """
    # Lookup metadata row
    if hasattr(metadata_df, 'loc'):  # pandas DataFrame
        try:
            meta = metadata_df[metadata_df['urn'] == urn].iloc[0].to_dict()
        except (IndexError, KeyError):
            return {
                "call": "NA",
                "standardized_score": None,
                "reason": "URN not found in metadata",
                "notes": urn,
                "metadata_flags": {},
                "urn": urn
            }
    else:  # dict lookup
        meta = metadata_df.get(urn)
        if not meta:
            return {
                "call": "NA",
                "standardized_score": None,
                "reason": "URN not found in metadata",
                "notes": urn,
                "metadata_flags": {},
                "urn": urn
            }
    
    # Check quality flags - reject problematic datasets immediately
    flags = {
        "ambiguous": meta.get('ambiguous', False),
        "nonmonotonic": meta.get('nonmonotonic', False),
        "needs_gene_rules": meta.get('needs_gene_specific_rules', False),
        "flexibility_required": meta.get('flexibility_requirement', False)
    }
    
    if flags["ambiguous"]:
        return {
            "call": "NA",
            "standardized_score": None,
            "reason": "Dataset marked ambiguous",
            "notes": meta.get('rationale', ''),
            "metadata_flags": flags,
            "urn": urn
        }
    
    if flags["nonmonotonic"] or flags["needs_gene_rules"]:
        return {
            "call": "NA",
            "standardized_score": None,
            "reason": "Gene-specific rules required (not supported)",
            "notes": meta.get('rationale', ''),
            "metadata_flags": flags,
            "urn": urn
        }
    
    # Parse neutral band from metadata (handles one-sided intervals)
    lo, hi = parse_neutral_band_from_metadata(meta)
    
    if lo is None or hi is None:
        return {
            "call": "NA",
            "standardized_score": None,
            "reason": "Unable to determine neutral band",
            "notes": f"neutral_bands={meta.get('neutral_bands')}, q05={meta.get('q05')}, q95={meta.get('q95')}",
            "metadata_flags": flags,
            "urn": urn
        }
    
    neutral_band = (lo, hi)
    
    # Build AssaySpec from metadata
    try:
        family = meta.get('assay_family')
        subtype_value = meta.get('subtype_value')
        
        # Get subtype enum class and member
        subtype_enum_cls = FAMILY_TO_SUBTYPE_ENUM.get(family)
        if not subtype_enum_cls:
            return {
                "call": "NA",
                "standardized_score": None,
                "reason": f"Unknown assay family: {family}",
                "notes": "",
                "metadata_flags": flags,
                "urn": urn
            }
        
        try:
            subtype = subtype_enum_cls[subtype_value]
        except KeyError:
            return {
                "call": "NA",
                "standardized_score": None,
                "reason": f"Unknown subtype {subtype_value} for {family}",
                "notes": "",
                "metadata_flags": flags,
                "urn": urn
            }
        
        # Parse directionality and scale
        directionality_str = meta.get('directionality', 'HIGHER_IS_FUNCTION')
        scale_str = meta.get('scale', 'CENTERED_0')
        
        try:
            directionality = Directionality[directionality_str]
            scale = Scale[scale_str]
        except KeyError as e:
            return {
                "call": "NA",
                "standardized_score": None,
                "reason": f"Invalid enum value: {e}",
                "notes": f"directionality={directionality_str}, scale={scale_str}",
                "metadata_flags": flags,
                "urn": urn
            }
        
        # Extract flexibility_requirement from metadata
        flexibility_requirement = meta.get('flexibility_requirement', False)
        
        spec = AssaySpec(
            assay_family=family,
            subtype=subtype,
            directionality=directionality,
            scale=scale,
            baseline=meta.get('baseline', ''),
            units=meta.get('units', ''),
            neutral_band=neutral_band,
            flexibility_requirement=flexibility_requirement  # Added from metadata
        )
    except (KeyError, AttributeError, TypeError) as e:
        return {
            "call": "NA",
            "standardized_score": None,
            "reason": f"Failed to construct AssaySpec: {e}",
            "notes": str(meta),
            "metadata_flags": flags,
            "urn": urn
        }
    
    # Call interpret_score (no longer needs function_prefers_lower_stability param)
    result = interpret_score(
        score=score,
        spec=spec,
    )
    
    # Add metadata to result
    result["metadata_flags"] = flags
    result["urn"] = urn
    result["n_syn_used"] = meta.get('n_syn_used', 0)
    result["band_source"] = meta.get('band_source', '')
    result["assay_family"] = family
    result["subtype"] = subtype_value
    result["neutral_band_parsed"] = neutral_band
    result["MaveDB_BS3"] = False
    result["MaveDB_PS3"] = False

        # === PS3 Logic: Any functional impact (LoF OR GoF) ===
    if result["call"] in ("LoF", "GoF"):
        result["MaveDB_PS3"] = True
    # === BS3 Logic: No functional impact (Neutral) ===
    elif result["call"] == "Neutral":
        result["MaveDB_BS3"] = True
    
    return result

