#!/usr/bin/env python3
"""Compare `priva gofcards-clinvar-injection` and `priva gofcards-cache-validate`
with inject_clinvar_into_gofcards.py and validate_gofcards_cache.py.

  python gofcards_clinvar_injection_compare.py check --binary PRIVA --repo CHECKOUT --work-root NEW_DIR [--jobs N]
  python gofcards_clinvar_injection_compare.py real --binary PRIVA --repo CHECKOUT --work-root NEW_DIR \
      --gofcards-cache JSON_GZ --clinvar-vcv XML_GZ --min-review-stars N \
      [--deployed JSON_GZ --audit TSV --stats JSON | --python-output DIR]

Each `check` case writes the same small hand-made VCV XML and GoFCards cache
into two folders and runs one program in each, with relative paths. Exit status,
failure category and message, log lines, stdout and every file left in the
folder are compared; gzip files by their bytes. The order of the confirmed_by_*
counts follows Python's string hash, so the statistics are compared as a
mapping and those keys are left out of the order check.

`real` runs the native command on the full XML and compares its cache, audit
and statistics with the output of a Python run: the deployed files, or the
folder written by a separate Python job.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import copy
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

SHA = "0" * 64


def variant(chrom, pos, ref, alt, hg38=None, status="eligible", variation_id=None, **extra):
    assemblies = {"hg19": {"genomic": {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt},
                           "transcripts": {"NM_1.1": {"by_hgvsc": {"c.1A>G": {}}, "by_hgvsp": {"p.K1E": ["c.1A>G"]}}}}}
    if hg38:
        assemblies["hg38"] = {"genomic": {"chrom": chrom, "pos": hg38, "ref": ref, "alt": alt}}
    eligibility = {"status": status}
    if status != "eligible":
        eligibility["reason"] = "reviewed"
    record = {"eligibility": eligibility, "evidence": [{"pmid": "1"}]}
    if variation_id is not None:
        record["annotations"] = {"clinvar_variation_id": variation_id}
    return {"record": record, "assemblies": assemblies, **extra}


def vid(chrom, pos, ref, alt):
    return f"loc_{chrom}:{pos}:{ref}->{alt}_grch37"


def base_cache():
    return {
        "metadata": {"source": "GoFCards", "mechanism": "GOF", "version": "é 1"},
        "genes": {
            "GENE1": {"hgnc_id": "HGNC:1", "variants": {
                vid(1, 100, "A", "G"): variant("1", 100, "A", "G", hg38=200),
                vid(1, 300, "C", "T"): variant("1", 300, "C", "T")}},
            "GENE2": {"variants": {vid(2, 500, "G", "A"): variant("2", 500, "G", "A", hg38=600, variation_id="3003")}},
            "GENE3": {"variants": {vid(3, 700, "T", "C"): variant("3", 700, "T", "C", status="quarantined_reviewed_lof")}},
            "GENE4": {"variants": {vid(1, 100, "A", "G"): variant("1", 100, "A", "G", hg38=200)}},
        },
    }


def location(assembly="GRCh37", chrom="1", pos=100, ref="A", alt="G"):
    return (f'<SequenceLocation Assembly="{assembly}" Chr="{chrom}" positionVCF="{pos}" '
            f'referenceAlleleVCF="{ref}" alternateAlleleVCF="{alt}"/>')


def allele(locations=None, genes=("GENE1",), hgvs=True, attributes='AlleleID="11" VariationID="1001"'):
    locations = [location()] if locations is None else locations
    gene_list = "".join(f'<Gene Symbol="{gene}" GeneID="5" HGNC_ID="HGNC:1" RelationshipType="within"/>' for gene in genes)
    expressions = ('<HGVSlist><HGVS Type="coding"><NucleotideExpression sequenceAccessionVersion="NM_1.1">'
                   '<Expression> NM_1.1:c.1A&gt;G </Expression></NucleotideExpression></HGVS>'
                   '<HGVS Type="protein"><ProteinExpression sequenceAccessionVersion="NP_1.1">'
                   '<Expression>NP_1.1:p.Lys1Glu</Expression></ProteinExpression></HGVS>'
                   '<HGVS Type="coding"><NucleotideExpression><Expression>NM_1.1:c.1A&gt;G</Expression>'
                   '</NucleotideExpression></HGVS><HGVS Type="empty"/></HGVSlist>') if hgvs else ""
    return (f'<SimpleAllele {attributes}><Name>a name</Name><VariantType>single nucleotide variant</VariantType>'
            f'<CanonicalSPDI>NC:1:A:G</CanonicalSPDI><GeneList>{gene_list}</GeneList>'
            f'<Location>{"".join(locations)}</Location>{expressions}</SimpleAllele>')


def rcv(accession="RCV1", status="criteria provided, multiple submitters, no conflicts", name="Disease One",
        condition_id="C0001", db="MedGen", version='Version="2"', submission='SubmissionCount="3"'):
    return (f'<RCVAccession Accession="{accession}" {version} Title="title &amp; more">'
            f'<ClassifiedConditionList TraitSetID="77"><ClassifiedCondition DB="{db}" ID="{condition_id}"> {name} '
            f'</ClassifiedCondition></ClassifiedConditionList><RCVClassifications><GermlineClassification>'
            f'<ReviewStatus>{status}</ReviewStatus><Description DateLastEvaluated="2024-01-01" {submission}>'
            f'Pathogenic</Description></GermlineClassification></RCVClassifications></RCVAccession>')


OBSERVED = ('<ObservedInList><ObservedIn><Sample><Origin>germline</Origin><AffectedStatus>yes</AffectedStatus>'
            '<Sex>female</Sex><NumberTested>12</NumberTested><NumberMales>x1</NumberMales><NumberFemales>0</NumberFemales>'
            '<NumberChrTested> 7 </NumberChrTested></Sample><ObservedData><Attribute Type="Description">note</Attribute>'
            '</ObservedData><ObservedData/><Co-occurrenceSet><Zygosity>Heterozygote</Zygosity><Count>3</Count>'
            '<AlleleDescSet><Name>b</Name><RelativeOrientation>trans</RelativeOrientation><Zygosity>Het</Zygosity>'
            '<ClinicalSignificance><Description>Benign</Description></ClinicalSignificance></AlleleDescSet>'
            '<AlleleDescSet><ClinicalSignificance>plain text</ClinicalSignificance></AlleleDescSet></Co-occurrenceSet>'
            '<Co-occurrenceSet><Zygosity>Heterozygote</Zygosity><Count>many</Count></Co-occurrenceSet>'
            '<Co-occurrenceSet><Zygosity>Homozygote</Zygosity></Co-occurrenceSet>'
            '<Co-occurrenceSet><Count>9</Count></Co-occurrenceSet></ObservedIn><ObservedIn/></ObservedInList>')


def assertion(identifier="9", accession="SCV1", status="current", contributes="true", germline=True,
              traits='<XRef DB="MedGen" ID="C0001"/>', name="Other", observed=OBSERVED, version='Version="4"'):
    classification = (f'<Classification DateLastEvaluated="2023-02-02"><ReviewStatus>criteria provided</ReviewStatus>'
                      f'{"<GermlineClassification>Pathogenic</GermlineClassification>" if germline else ""}</Classification>')
    return (f'<ClinicalAssertion ID="{identifier}" ContributesToAggregateClassification="{contributes}">'
            f'<ClinVarAccession Accession="{accession}" {version} SubmitterName="Lab &#233;" OrgID="55" '
            f'OrganizationCategory="laboratory"/><RecordStatus>{status}</RecordStatus>{classification}'
            f'<AttributeSet><Attribute Type="ModeOfInheritance">Dominant</Attribute></AttributeSet>'
            f'<AttributeSet><Attribute Type="ModeOfInheritance">Dominant</Attribute></AttributeSet>'
            f'<AttributeSet><Attribute Type="ModeOfInheritance">Autosomal</Attribute></AttributeSet>'
            f'<AttributeSet><Attribute Type="Penetrance">low</Attribute></AttributeSet><AttributeSet/>'
            f'{observed}<TraitSet><Trait>{traits}<Name><ElementValue>  {name}  </ElementValue></Name></Trait></TraitSet>'
            f'</ClinicalAssertion>')


MAPPING = ('<TraitMappingList><TraitMapping ClinicalAssertionID="9" TraitType="Disease" MappingType="Name" '
           'MappingValue="v" MappingRef="r"><MedGen Name="Disease   one" CUI="C9999"/></TraitMapping>'
           '<TraitMapping TraitType="Disease"/><TraitMapping ClinicalAssertionID="8"/></TraitMappingList>')


def archive(identifier="1001", accession="VCV1", status="current", body=None, alleles=None, rcvs=None,
            assertions=None, mapping=MAPPING, tag="ClassifiedRecord"):
    alleles = allele() if alleles is None else alleles
    rcvs = rcv() if rcvs is None else rcvs
    assertions = assertion() if assertions is None else assertions
    if body is None:
        body = (f'<{tag}>{alleles}<RCVList>{rcvs}</RCVList><ClinicalAssertionList>{assertions}'
                f'</ClinicalAssertionList>{mapping}</{tag}>')
    return (f'<VariationArchive VariationID="{identifier}" Accession="{accession}" Version="3" '
            f'VariationName=" name {identifier} " VariationType="SNV" RecordType="classified">'
            f'<RecordStatus>{status}</RecordStatus>{body}</VariationArchive>\n')


def document(*archives, head='<?xml version="1.0"?>\n<ClinVarVariationRelease ReleaseDate="2026-01-01">',
             tail="</ClinVarVariationRelease>\n"):
    return head + "\n" + "".join(archives) + tail


def base_spec():
    return {"cache": base_cache(), "cache_name": "cache.json.gz", "xml": document(archive()), "xml_name": "vcv.xml",
            "out": "out.json.gz", "args": ["--audit-tsv", "audit.tsv", "--stats-json", "stats.json"]}


CASES = []


def case(function):
    CASES.append(function)
    return function


def set_xml(spec, *archives):
    spec["xml"] = document(*archives)


@case
def base(spec): pass
@case
def gz_xml(spec): spec["xml_name"] = "vcv.xml.gz"
@case
def plain_cache(spec): spec["cache_name"] = "cache.json"
@case
def plain_output(spec): spec["out"] = "out.json"
@case
def output_without_suffix(spec): spec["out"] = "deployed"
@case
def output_in_new_folder(spec): spec["out"] = "new/deep/out.json.gz"
@case
def no_audit_or_stats(spec): spec["args"] = []
@case
def audit_in_new_folder(spec): spec["args"] = ["--audit-tsv", "reports/deep/audit.tsv"]
@case
def stats_parent_missing(spec): spec["args"] = ["--stats-json", "reports/stats.json"]
@case
def unnormalized_paths(spec):
    spec["xml_path"] = "./sub//vcv.xml"
    spec["xml_name"] = "sub/vcv.xml"
@case
def stars_3_rejects(spec): spec["args"] += ["--min-review-stars", "3"]
@case
def stars_3_accepts(spec):
    set_xml(spec, archive(rcvs=rcv(status="reviewed by expert panel")))
    spec["args"] += ["--min-review-stars", "3"]
@case
def stars_4_rejects_panel(spec):
    set_xml(spec, archive(rcvs=rcv(status="reviewed by expert panel")))
    spec["args"] += ["--min-review-stars", "4"]
@case
def stars_4_accepts(spec):
    set_xml(spec, archive(rcvs=rcv(status="practice guideline")))
    spec["args"] += ["--min-review-stars=4"]
@case
def stars_abbreviated_option(spec): spec["args"] += ["--min-review", "2"]
@case
def stars_low_status(spec): set_xml(spec, archive(rcvs=rcv(status="criteria provided, single submitter")))
@case
def stars_unknown_status(spec): set_xml(spec, archive(rcvs=rcv(status="something else")))
@case
def stars_status_spacing_and_case(spec):
    set_xml(spec, archive(rcvs=rcv(status="  Criteria Provided,\n  MULTIPLE submitters,\tno conflicts ")))
@case
def several_rcvs_mixed_stars(spec):
    set_xml(spec, archive(rcvs=rcv("RCV1") + rcv("RCV2", status="no assertion criteria provided", name="Two")
                          + rcv("RCV3", name="Disease Three", condition_id="C0003")))
@case
def rcv_without_germline(spec):
    set_xml(spec, archive(rcvs='<RCVAccession Accession="RCV0"/>' + rcv()))
@case
def rcv_without_conditions(spec):
    set_xml(spec, archive(rcvs='<RCVAccession Accession="RCV0" Title="t"><RCVClassifications><GermlineClassification>'
                               '<ReviewStatus>practice guideline</ReviewStatus></GermlineClassification>'
                               '</RCVClassifications></RCVAccession>'))
@case
def rcv_version_text_and_absent(spec):
    set_xml(spec, archive(rcvs=rcv("RCV1", version='Version="v2"', submission='SubmissionCount=" 8 "')
                          + rcv("RCV2", name="Two", version="", submission="")))
@case
def hg38_match_only(spec): set_xml(spec, archive(alleles=allele([location("GRCh38", pos=200)])))
@case
def hg19_match_when_liftover_failed(spec):
    set_xml(spec, archive(alleles=allele([location(pos=300, ref="C", alt="T"), location("GRCh38", pos=999)])))
@case
def both_builds_match(spec): set_xml(spec, archive(alleles=allele([location(), location("GRCh38", pos=200)])))
@case
def unknown_assembly_and_missing_fields(spec):
    set_xml(spec, archive(alleles=allele(['<SequenceLocation Assembly="NCBI36" Chr="1" positionVCF="100" referenceAlleleVCF="A" alternateAlleleVCF="G"/>',
                                          '<SequenceLocation Assembly="GRCh37" Chr="1" positionVCF="100" referenceAlleleVCF="A"/>',
                                          location()])))
@case
def later_location_of_an_assembly_wins(spec):
    set_xml(spec, archive(alleles=allele([location(), location(pos=999)])))
@case
def nested_location_is_found(spec):
    set_xml(spec, archive(alleles=allele([f"<Wrapped>{location()}</Wrapped>"])))
@case
def chr_prefix_mt_and_lowercase(spec):
    spec["cache"]["genes"]["GENE1"]["variants"][vid("M", 50, "A", "G")] = variant("chrMT", 50, "a", "g")
    set_xml(spec, archive(alleles=allele([location(chrom=" MT ", pos=" 50 ", ref="a", alt=" g ")])))
@case
def variation_id_route_only(spec):
    set_xml(spec, archive("3003", alleles=allele([location(pos=999)], genes=("GENE2",))))
@case
def variation_id_with_spaces(spec):
    set_xml(spec, archive(" 3003 ", alleles=allele([location(pos=999)], genes=("GENE2",))))
@case
def both_routes_confirm(spec):
    set_xml(spec, archive("3003", alleles=allele([location(chrom="2", pos=500, ref="G", alt="A")], genes=("GENE2",))))
@case
def variation_id_number_in_cache(spec):
    spec["cache"]["genes"]["GENE2"]["variants"][vid(2, 500, "G", "A")]["record"]["annotations"]["clinvar_variation_id"] = 3003
    set_xml(spec, archive("3003", alleles=allele([location(pos=999)], genes=("GENE2",))))
@case
def variation_id_gene_discordant(spec):
    set_xml(spec, archive("3003", alleles=allele([location(pos=999)], genes=("OTHER",))))
@case
def quarantined_variant_never_matches(spec):
    set_xml(spec, archive(alleles=allele([location(chrom="3", pos=700, ref="T", alt="C")], genes=("GENE3",))))
@case
def quarantined_variation_id_never_matches(spec):
    spec["cache"]["genes"]["GENE3"]["variants"][vid(3, 700, "T", "C")]["record"]["annotations"] = {"clinvar_variation_id": "7"}
    set_xml(spec, archive("7", alleles=allele([location(pos=999)], genes=("GENE3",))))
@case
def two_genes_at_one_allele(spec): set_xml(spec, archive(alleles=allele(genes=("GENE1", "GENE4"))))
@case
def second_gene_only(spec): set_xml(spec, archive(alleles=allele(genes=("GENE4",))))
@case
def no_gene_list_keeps_every_match(spec): set_xml(spec, archive(alleles=allele(genes=())))
@case
def gene_discordant(spec): set_xml(spec, archive(alleles=allele(genes=("OTHER",))))
@case
def no_match(spec): set_xml(spec, archive(alleles=allele([location(pos=999)])))
@case
def noncurrent_archive(spec): set_xml(spec, archive(status="previous"))
@case
def included_record(spec): set_xml(spec, archive(tag="IncludedRecord"))
@case
def status_spacing(spec): set_xml(spec, archive(status="  current\n"))
@case
def no_alleles(spec): set_xml(spec, archive(body="<ClassifiedRecord/>"))
@case
def no_eligible_rcv(spec): set_xml(spec, archive(rcvs=""))
@case
def haplotype(spec):
    set_xml(spec, archive(body=f'<ClassifiedRecord><Haplotype>{allele()}{allele([location(pos=300, ref="C", alt="T")])}'
                               f'{allele([location(pos=999)])}</Haplotype><RCVList>{rcv()}</RCVList></ClassifiedRecord>'))
@case
def genotype(spec):
    set_xml(spec, archive(body=f'<ClassifiedRecord><Genotype>{allele()}<Haplotype>{allele([location(pos=300, ref="C", alt="T")])}'
                               f'</Haplotype><Other/></Genotype><RCVList>{rcv()}</RCVList></ClassifiedRecord>'))
@case
def simple_allele_beats_haplotype(spec):
    set_xml(spec, archive(body=f'<ClassifiedRecord><Haplotype>{allele([location(pos=300, ref="C", alt="T")])}</Haplotype>'
                               f'{allele()}<RCVList>{rcv()}</RCVList></ClassifiedRecord>'))
@case
def two_records_one_variant(spec):
    set_xml(spec, archive(), archive("1002", "VCV2", rcvs=rcv("RCV9")), archive("1003", "VCV3"), archive("1002", "VCV2"))
@case
def same_record_twice(spec): set_xml(spec, archive(), archive())
@case
def existing_clinvar_block(spec):
    spec["cache"]["genes"]["GENE4"]["variants"][vid(1, 100, "A", "G")]["clinvar"] = {"vcv_accession": "VCV0"}
@case
def existing_additional(spec):
    block = spec["cache"]["genes"]["GENE4"]["variants"][vid(1, 100, "A", "G")]
    block["clinvar"] = {"vcv_accession": "VCV0"}
    block["clinvar_additional"] = [{"vcv_accession": "VCV1"}]
@case
def existing_null_clinvar(spec):
    spec["cache"]["genes"]["GENE4"]["variants"][vid(1, 100, "A", "G")]["clinvar"] = None
@case
def existing_same_accession(spec):
    spec["cache"]["genes"]["GENE4"]["variants"][vid(1, 100, "A", "G")]["clinvar"] = {"vcv_accession": "VCV1"}
@case
def existing_clinvar_not_a_dict(spec):
    spec["cache"]["genes"]["GENE4"]["variants"][vid(1, 100, "A", "G")]["clinvar"] = "text"
@case
def variant_key_with_spaces(spec):
    genes = spec["cache"]["genes"]
    genes["GENE1"]["variants"] = {" spaced ": variant("1", 100, "A", "G")}
    del genes["GENE4"]
@case
def trait_mapping_links(spec):
    set_xml(spec, archive(rcvs=rcv(name="Disease One", condition_id=""),
                          assertions=assertion(traits="", name="Nothing")))
@case
def trait_mapping_by_cui(spec):
    set_xml(spec, archive(rcvs=rcv(condition_id="c9999", db=""), assertions=assertion(traits="", name="Nothing")))
@case
def trait_name_links(spec):
    set_xml(spec, archive(rcvs=rcv(name="Trait   Name", condition_id=""), mapping="",
                          assertions=assertion(traits="", name="trait\n name")))
@case
def trait_database_identifier_links(spec):
    set_xml(spec, archive(rcvs=rcv(condition_id="OMIM:1", db="OMIM"), mapping="",
                          assertions=assertion(traits='<XRef DB="OMIM" ID="omim:1"/>')))
@case
def scv_unlinked(spec):
    set_xml(spec, archive(mapping="", assertions=assertion(traits='<XRef DB="X" ID="Z"/>')))
@case
def scv_filters(spec):
    set_xml(spec, archive(assertions=assertion(status="replaced", identifier="1") + assertion(contributes="false", identifier="2")
                          + assertion(germline=False, identifier="3") + assertion(contributes="TRUE", identifier="4")
                          + '<ClinicalAssertion ID="5" ContributesToAggregateClassification="true"><RecordStatus>current</RecordStatus></ClinicalAssertion>'))
@case
def scv_ambiguous_across_rcvs(spec):
    set_xml(spec, archive(rcvs=rcv("RCV2") + rcv("RCV1"), assertions=assertion() + assertion("10", "SCV2")))
@case
def scv_without_identifier(spec):
    set_xml(spec, archive(rcvs=rcv("RCV2") + rcv("RCV1"), assertions=assertion("", "") + assertion("", "SCV3")))
@case
def scv_version_text(spec): set_xml(spec, archive(assertions=assertion(version='Version="x"')))
@case
def scv_minimal(spec):
    set_xml(spec, archive(assertions='<ClinicalAssertion ID="9" ContributesToAggregateClassification="true">'
                                     '<RecordStatus>current</RecordStatus><Classification><GermlineClassification/>'
                                     '</Classification><TraitSet><Trait><XRef ID="C0001"/></Trait></TraitSet></ClinicalAssertion>'))
@case
def unicode_and_entities(spec):
    set_xml(spec, archive(rcvs=rcv(name="Maladie &#233;&#x1F600; &lt;x&gt; <![CDATA[a & b]]>", condition_id="Ⅷ")))
@case
def text_interrupted_by_child_and_comment(spec):
    set_xml(spec, archive(rcvs=rcv(name='before <!-- note --> middle <Sub>inner</Sub> tail')))
@case
def nonascii_space_in_text(spec): set_xml(spec, archive(rcvs=rcv(name="\u00a0Disease\u2003One\u00a0")))
@case
def default_namespace(spec):
    spec["xml"] = document(archive(), head='<?xml version="1.0"?>\n<ClinVarVariationRelease xmlns="urn:x">')
@case
def prefixed_tags(spec):
    body = archive().replace("<", "<p:").replace("<p:/", "</p:")
    spec["xml"] = '<p:Release xmlns:p="urn:x">' + body + "</p:Release>"
@case
def windows_1252_declaration(spec):
    spec["xml_bytes"] = document(archive(rcvs=rcv(name="Caf\u00e9")),
                                 head='<?xml version="1.0" encoding="windows-1252"?>\n<R>', tail="</R>").encode("cp1252")
@case
def utf16_document(spec):
    spec["xml_bytes"] = document(archive(), head='<?xml version="1.0" encoding="UTF-16"?>\n<R>', tail="</R>").encode("utf-16")
@case
def internal_entity(spec):
    spec["xml"] = document(archive(rcvs=rcv(name="&word;")),
                           head='<?xml version="1.0"?>\n<!DOCTYPE R [<!ENTITY word "Disease One">]>\n<R>', tail="</R>")
@case
def root_is_the_archive(spec): spec["xml"] = archive()
@case
def nested_archive(spec):
    inner = archive("1002", "VCV2")
    spec["xml"] = document(archive(body=inner + f'<ClassifiedRecord>{allele()}<RCVList>{rcv()}</RCVList></ClassifiedRecord>'))
@case
def no_archives(spec): spec["xml"] = document()
@case
def archive_spans_feeds(spec):
    filler = "<!--" + "x" * 40000 + "-->"
    set_xml(spec, archive(), filler, archive("1002", "VCV2", rcvs=rcv(name=" " * 20000 + "Disease One")),
            '<Note>' + "text " * 9000 + '</Note>', archive("1003", "VCV3", alleles=allele([location(pos=999)])))
@case
def progress_line(spec):
    miss = archive("5", alleles=allele([location(pos=999)], hgvs=False), rcvs="", assertions="", mapping="")
    set_xml(spec, *([miss] * 100_001 + [archive()]))
@case
def malformed_after_a_match(spec): spec["xml"] = document(archive())[:-40]
@case
def mismatched_tag(spec): spec["xml"] = "<R>\n  <A></B>\n</R>"
@case
def undefined_entity(spec): set_xml(spec, archive(rcvs=rcv(name="&nope;")))
@case
def trailing_content(spec): spec["xml"] = document(archive()) + "<Second/>"
@case
def not_xml(spec): spec["xml"] = "just some text\n"
@case
def empty_document_in_gzip(spec):
    spec["xml"], spec["xml_name"] = "", "vcv.xml.gz"
@case
def gz_name_not_gzip(spec):
    spec["raw_xml"], spec["xml_name"] = True, "vcv.xml.gz"
@case
def truncated_gzip(spec):
    spec["xml_name"], spec["truncate"] = "vcv.xml.gz", 20
@case
def missing_xml(spec): spec["no_xml"] = True
@case
def empty_xml(spec): spec["xml_bytes"], spec["raw_xml"] = b"", True
@case
def xml_is_a_directory(spec): spec["xml_dir"] = True
@case
def digest_accepted(spec): spec["digest"] = "right"
@case
def digest_recorded_for_gz(spec):
    spec["digest"], spec["xml_name"] = "right", "vcv.xml.gz"
@case
def digest_mismatch(spec): spec["args"] += ["--clinvar-vcv-sha256", SHA]
@case
def digest_uppercase(spec): spec["args"] += ["--clinvar-vcv-sha256", "A" * 64]
@case
def digest_short(spec): spec["args"] += ["--clinvar-vcv-sha256", "abc"]
@case
def digest_bad_and_xml_missing(spec):
    spec["no_xml"] = True
    spec["args"] += ["--clinvar-vcv-sha256", "abc"]
@case
def missing_cache(spec): spec["no_cache"] = True
@case
def cache_not_json(spec): spec["cache_text"] = "not json"
@case
def cache_no_genes(spec): del spec["cache"]["genes"]
@case
def cache_list(spec): spec["cache"] = ["genes"]
@case
def cache_number(spec): spec["cache"] = 3
@case
def cache_wrong_source(spec): spec["cache"]["metadata"]["source"] = "Other"
@case
def cache_wrong_mechanism(spec): spec["cache"]["metadata"]["mechanism"] = "LOF"
@case
def cache_no_metadata(spec): del spec["cache"]["metadata"]
@case
def cache_genes_list(spec): spec["cache"]["genes"] = ["x"]
@case
def cache_genes_empty(spec): spec["cache"]["genes"] = {}
@case
def cache_gene_not_dict(spec): spec["cache"]["genes"]["GENE1"] = 3
@case
def cache_variants_list(spec): spec["cache"]["genes"]["GENE1"]["variants"] = [1]
@case
def cache_variant_not_dict(spec): spec["cache"]["genes"]["GENE1"]["variants"]["x"] = 3
@case
def cache_assemblies_list(spec):
    spec["cache"]["genes"]["GENE1"]["variants"][vid(1, 100, "A", "G")]["assemblies"] = [1]
@case
def cache_block_not_dict(spec):
    spec["cache"]["genes"]["GENE1"]["variants"][vid(1, 100, "A", "G")]["assemblies"]["hg19"] = 3
@case
def cache_pipe_in_allele(spec):
    spec["cache"]["genes"]["GENE1"]["variants"][vid(1, 100, "A", "G")]["assemblies"]["hg19"]["genomic"]["alt"] = "G|T"
@case
def cache_numeric_and_null_alleles(spec):
    genomic = spec["cache"]["genes"]["GENE1"]["variants"][vid(1, 100, "A", "G")]["assemblies"]["hg19"]["genomic"]
    genomic.update(chrom=1, ref=None)
    del genomic["alt"]
@case
def cache_no_position(spec):
    del spec["cache"]["genes"]["GENE1"]["variants"][vid(1, 100, "A", "G")]["assemblies"]["hg19"]["genomic"]["pos"]
@case
def cache_float_position(spec):
    spec["cache"]["genes"]["GENE1"]["variants"][vid(1, 100, "A", "G")]["assemblies"]["hg19"]["genomic"]["pos"] = 100.0
@case
def cache_record_not_dict(spec):
    spec["cache"]["genes"]["GENE1"]["variants"][vid(1, 100, "A", "G")]["record"] = 3
@case
def cache_nan_kept(spec):
    spec["cache_text"] = json.dumps(base_cache()).replace('"é 1"', "NaN").replace('"hgnc_id": "HGNC:1"', '"hgnc_id": 1.5e300')
@case
def argument_stars_text(spec): spec["args"] += ["--min-review-stars", "two"]
@case
def argument_stars_choice(spec): spec["args"] += ["--min-review-stars", "5"]
@case
def argument_unknown(spec): spec["args"] += ["--unknown", "x"]
@case
def argument_positional(spec): spec["args"] += ["extra"]
@case
def argument_no_value(spec): spec["args"] += ["--audit-tsv"]
@case
def argument_required(spec): spec["drop"] = ["--out-json"]
@case
def argument_required_two(spec): spec["drop"] = ["--clinvar-vcv", "--gofcards-cache"]


def write_inputs(folder, spec):
    folder.mkdir(parents=True)
    cache = folder / spec["cache_name"]
    if not spec.get("no_cache"):
        text = spec.get("cache_text") or json.dumps(spec["cache"], ensure_ascii=False)
        if spec["cache_name"].endswith(".gz"):
            cache.write_bytes(gzip.compress(text.encode(), mtime=0))
        else:
            cache.write_text(text, encoding="utf-8")
    xml = folder / spec["xml_name"]
    xml.parent.mkdir(parents=True, exist_ok=True)
    if spec.get("xml_dir"):
        xml.mkdir()
    elif not spec.get("no_xml"):
        data = spec.get("xml_bytes", spec["xml"].encode() if "xml_bytes" not in spec else b"")
        if spec["xml_name"].endswith(".gz") and not spec.get("raw_xml"):
            data = gzip.compress(data, mtime=0)
        if spec.get("truncate"):
            data = data[:spec["truncate"]]
        xml.write_bytes(data)


def arguments(spec, folder):
    values = {"--gofcards-cache": spec["cache_name"], "--clinvar-vcv": spec.get("xml_path", spec["xml_name"]),
              "--out-json": spec["out"]}
    for name in spec.get("drop", []):
        del values[name]
    command = [part for pair in values.items() for part in pair] + list(spec["args"])
    if spec.get("digest") == "right":
        command += ["--clinvar-vcv-sha256", hashlib.sha256((folder / spec["xml_name"]).read_bytes()).hexdigest()]
    return command


def outcome(command, folder, env):
    result = subprocess.run([str(part) for part in command], cwd=folder, capture_output=True, text=True, env=env,
                            errors="replace")
    lines = result.stderr.splitlines()
    value = {"exit": result.returncode, "stdout": result.stdout}
    log = [line for line in lines if line.startswith(("[clinvar-inject]", "clinvar_vcv_progress"))]
    value["log"] = [line if not line.startswith("[clinvar-inject] injected: ")
                    else "injected: " + "; ".join(sorted(line.split("injected: ", 1)[1].split("; "))) for line in log]
    if result.returncode == 2 and lines:
        value["message"] = lines[-1].split("error: ", 1)[-1].removeprefix("priva: ")
    if result.returncode == 1 and lines:
        last = lines[-1]
        native = re.fullmatch(r"priva: ([A-Za-z]+): (.*)", last)
        if native:
            value["category"], value["message"] = native[1], native[2]
        elif "Traceback (most recent call last)" in result.stderr:
            name, _, message = last.partition(": ")
            value["category"], value["message"] = name.rsplit(".", 1)[-1], message
        else:
            value["category"], value["message"] = "SystemExit", last
    return value


def snapshot(folder):
    files = {}
    for path in sorted(folder.rglob("*")):
        if path.is_file():
            relative = str(path.relative_to(folder))
            if path.name == "stats.json":
                pairs = json.loads(path.read_text(), object_pairs_hook=list)
                files[relative] = {"mapping": sorted(pairs),
                                   "order": [key for key, _ in pairs if not key.startswith("confirmed_by_")],
                                   "indented": "\n  " in path.read_text()}
            else:
                files[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    return files


def environment():
    env = dict(os.environ)
    env["PYTHONHASHSEED"] = "0"
    return env


def compare(name, python_value, native_value, folders):
    same = python_value == native_value
    if same:
        return None
    detail = {}
    for key in sorted(set(python_value) | set(native_value)):
        if python_value.get(key) != native_value.get(key):
            detail[key] = [python_value.get(key), native_value.get(key)]
    return {"case": name, "differences": detail, "folders": [str(folder) for folder in folders]}


def run_injection(function, binary, repo, root):
    spec = base_spec()
    function(spec)
    left, right = root / function.__name__ / "python", root / function.__name__ / "native"
    write_inputs(left, spec)
    write_inputs(right, spec)
    env = environment()
    script = repo / "scripts" / "inject_clinvar_into_gofcards.py"
    first = outcome([sys.executable, script, *arguments(spec, left)], left, env)
    second = outcome([binary, "gofcards-clinvar-injection", "--source-root", repo, *arguments(spec, right)], right, env)
    first["files"], second["files"] = snapshot(left), snapshot(right)
    failure = compare(function.__name__, first, second, (left, right))
    if failure is None:
        shutil.rmtree(left.parent)
    return failure, first["exit"]


def good_cache():
    return base_cache()


def validator_cases():
    variants = lambda cache: cache["genes"]["GENE1"]["variants"]
    first = vid(1, 100, "A", "G")
    def edit(function):
        cache = good_cache()
        function(cache)
        return cache
    cases = {
        "good": (good_cache(), []),
        "provenance_skip_empty": (good_cache(), ["version=", "nothing"]),
        "provenance_match": (good_cache(), ["version=é 1", "source=GoFCards"]),
        "provenance_mismatch": (good_cache(), ["version=other"]),
        "provenance_missing": (good_cache(), ["clinvar.injector_sha256=abc"]),
        "provenance_nested": (edit(lambda c: c["metadata"].update(clinvar={"min_review_stars": 2})), ["clinvar.min_review_stars=2"]),
        "provenance_nested_mismatch": (edit(lambda c: c["metadata"].update(clinvar={"a": [1, "x", None]})), ["clinvar.a=2"]),
        "provenance_through_text": (good_cache(), ["version.deep=2"]),
        "provenance_equals_in_value": (edit(lambda c: c["metadata"].update(k="a=b")), ["k=a=b"]),
        "no_metadata": (edit(lambda c: c.pop("metadata")), []),
        "no_genes": (edit(lambda c: c.pop("genes")), []),
        "list_cache": (["metadata", "genes"], []),
        "list_without_keys": ([1], []),
        "text_cache": ("metadata genes", []),
        "number_cache": (3, []),
        "wrong_source": (edit(lambda c: c["metadata"].update(source="X")), []),
        "wrong_mechanism": (edit(lambda c: c["metadata"].update(mechanism="X")), []),
        "metadata_list": (edit(lambda c: c.update(metadata=[])), []),
        "genes_list": (edit(lambda c: c.update(genes=[])), []),
        "blank_symbol": (edit(lambda c: c["genes"].update({"  ": {}})), []),
        "gene_not_dict": (edit(lambda c: c["genes"].update(GENE1=3)), []),
        "no_variants": (edit(lambda c: c.update(genes={"G": {}})), []),
        "no_eligible": (edit(lambda c: c.update(genes={"G": {"variants": {vid(3, 700, "T", "C"): variant("3", 700, "T", "C", status="quarantined_reviewed_lof")}}})), []),
        "no_protein": (edit(lambda c: [block.update(transcripts={}) for gene in c["genes"].values() for v in gene["variants"].values() for block in v["assemblies"].values()]), []),
        "bad_variant_id": (edit(lambda c: variants(c).update(x=variant("1", 1, "A", "G"))), []),
        "variant_id_newline": (edit(lambda c: variants(c).update({vid(1, 5, "A", "G") + "\n": variant("1", 5, "A", "G")})), []),
        "variant_id_arabic_digit": (edit(lambda c: variants(c).update({"loc_1:\u0661:A->G_grch37": variant("1", 5, "A", "G")})), []),
        "variant_id_superscript": (edit(lambda c: variants(c).update({"loc_1:\u00b2:A->G_grch37": variant("1", 5, "A", "G")})), []),
        "variant_id_hyphen_alleles": (edit(lambda c: variants(c).update({"loc_X:9:--->TA_grch37": variant("1", 5, "A", "G")})), []),
        "variant_id_no_reference": (edit(lambda c: variants(c).update({"loc_X:9:->TA_grch37": variant("1", 5, "A", "G")})), []),
        "variant_id_two_arrows": (edit(lambda c: variants(c).update({"loc_X:9:A->T->G_grch37": variant("1", 5, "A", "G")})), []),
        "unknown_state": (edit(lambda c: variants(c)[first]["record"]["eligibility"].update(status="odd")), []),
        "state_number": (edit(lambda c: variants(c)[first]["record"]["eligibility"].update(status=3)), []),
        "state_absent": (edit(lambda c: variants(c)[first]["record"].pop("eligibility")), []),
        "no_assembly": (edit(lambda c: variants(c)[first].update(assemblies={})), []),
        "unknown_assembly": (edit(lambda c: variants(c)[first]["assemblies"].update(hg17={})), []),
        "no_coordinates": (edit(lambda c: variants(c)[first]["assemblies"]["hg19"].pop("genomic")), []),
        "position_text": (edit(lambda c: variants(c)[first]["assemblies"]["hg19"]["genomic"].update(pos="1")), []),
        "position_bool": (edit(lambda c: variants(c)[first]["assemblies"]["hg19"]["genomic"].update(pos=True)), []),
        "no_alternate": (edit(lambda c: variants(c)[first]["assemblies"]["hg19"]["genomic"].update(alt="")), []),
        "protein_key_empty": (edit(lambda c: variants(c)[first]["assemblies"]["hg19"]["transcripts"]["NM_1.1"]["by_hgvsp"].update({"p.X": []})), []),
        "protein_key_text": (edit(lambda c: variants(c)[first]["assemblies"]["hg19"]["transcripts"]["NM_1.1"]["by_hgvsp"].update({"p.X": "c.1"})), []),
        "protein_key_unknown": (edit(lambda c: variants(c)[first]["assemblies"]["hg19"]["transcripts"]["NM_1.1"]["by_hgvsp"].update({"p.X": ["c.9"]})), []),
        "quarantined_without_reason": (edit(lambda c: variants(c)[first]["record"].update(eligibility={"status": "quarantined_reviewed_lof"})), []),
        "evidence_number": (edit(lambda c: variants(c)[first]["record"].update(evidence=3)), []),
        "evidence_text": (edit(lambda c: variants(c)[first]["record"].update(evidence="abc")), []),
        "failed_liftover_counts": (good_cache(), []),
    }
    return cases


def run_validator(name, content, provenance, binary, repo, root, suffix):
    outcomes, folders = [], []
    for side, command in (("python", [sys.executable, repo / "scripts" / "validate_gofcards_cache.py"]),
                          ("native", [binary, "gofcards-cache-validate"])):
        folder = root / f"validate_{name}{suffix.replace('.', '_')}" / side
        folder.mkdir(parents=True)
        path = folder / f"cache{suffix}"
        text = json.dumps(content, ensure_ascii=False)
        if suffix.endswith(".gz"):
            path.write_bytes(gzip.compress(text.encode()))
        else:
            path.write_text(text, encoding="utf-8")
        outcomes.append(outcome([*command, path.name, *provenance], folder, environment()))
        folders.append(folder)
    failure = compare(f"validate_{name}{suffix}", outcomes[0], outcomes[1], folders)
    if failure is None:
        shutil.rmtree(folders[0].parent)
    return failure, outcomes[0]["exit"]


def check(args):
    binary, repo = Path(args.binary).resolve(), Path(args.repo).resolve()
    root = Path(args.work_root)
    root.mkdir(parents=True)
    root = root.resolve()
    wanted = set(filter(None, args.cases.split(",")))
    failures, statuses, count = [], {}, 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        jobs = [pool.submit(run_injection, function, binary, repo, root)
                for function in CASES if not wanted or function.__name__ in wanted]
        for name, (content, provenance) in validator_cases().items():
            if not wanted or f"validate_{name}" in wanted:
                for suffix in (".json.gz", ".json"):
                    jobs.append(pool.submit(run_validator, name, content, provenance, binary, repo, root, suffix))
        for job in jobs:
            failure, status = job.result()
            count += 1
            statuses[str(status)] = statuses.get(str(status), 0) + 1
            if failure:
                failures.append(failure)
    for failure in failures:
        print(json.dumps(failure, ensure_ascii=False)[:1500])
    print(json.dumps({"cases": count, "exit_statuses": statuses, "failed": [item["case"] for item in failures]}))
    return 1 if failures else 0


def load(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def real(args):
    binary, repo = Path(args.binary).resolve(), Path(args.repo).resolve()
    root = Path(args.work_root)
    root.mkdir(parents=True)
    out, audit, stats = root / "out.json.gz", root / "audit.tsv", root / "stats.json"
    command = [binary, "gofcards-clinvar-injection", "--source-root", repo, "--gofcards-cache", args.gofcards_cache,
               "--clinvar-vcv", args.clinvar_vcv, "--out-json", out, "--min-review-stars", args.min_review_stars,
               "--audit-tsv", audit, "--stats-json", stats]
    if args.digest:
        command += ["--clinvar-vcv-sha256", args.digest]
    started = time.time()
    subprocess.run(["/usr/bin/time", "-v", *map(str, command)], check=True, stderr=(root / "native.log").open("w"))
    report = {"native_seconds": round(time.time() - started, 1)}
    if args.python_output:
        base = Path(args.python_output)
        reference = {"cache": base / "out.json.gz", "audit": base / "audit.tsv", "stats": base / "stats.json"}
    else:
        reference = {"cache": Path(args.deployed), "audit": Path(args.audit), "stats": Path(args.stats)}
    left, right = load(reference["cache"]), load(out)
    report["cache_same"] = left == right
    if not report["cache_same"]:
        lc, rc = left["metadata"].get("clinvar", {}), right["metadata"].get("clinvar", {})
        report["metadata_differences"] = {key: [lc.get(key), rc.get(key)] for key in sorted(set(lc) | set(rc)) if lc.get(key) != rc.get(key)}
        left["metadata"].pop("clinvar", None)
        right["metadata"].pop("clinvar", None)
        report["cache_same_apart_from_clinvar_metadata"] = left == right
    report["gzip_bytes_same"] = reference["cache"].read_bytes() == out.read_bytes()
    report["audit_same"] = reference["audit"].read_bytes() == audit.read_bytes()
    report["stats_same"] = json.loads(reference["stats"].read_text()) == json.loads(stats.read_text())
    report["stats"] = json.loads(stats.read_text())
    report["validate_exit"] = subprocess.run([binary, "gofcards-cache-validate", out], capture_output=True).returncode
    report["same"] = bool(report["cache_same"] and report["audit_same"] and report["stats_same"])
    print(json.dumps(report))
    return 0 if report["same"] else 1


def main():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("check", "real"):
        part = commands.add_parser(name)
        part.add_argument("--binary", required=True)
        part.add_argument("--repo", required=True)
        part.add_argument("--work-root", required=True)
        if name == "check":
            part.add_argument("--jobs", type=int, default=4)
            part.add_argument("--cases", default="", help="comma-separated case names (default: all)")
        else:
            for option in ("gofcards-cache", "clinvar-vcv", "min-review-stars"):
                part.add_argument("--" + option, required=True)
            for option in ("digest", "deployed", "audit", "stats", "python-output"):
                part.add_argument("--" + option)
    args = parser.parse_args()
    return check(args) if args.command == "check" else real(args)


if __name__ == "__main__":
    sys.exit(main())
