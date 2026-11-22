#!/usr/bin/env python3
"""Writes the small made-up inputs of the Rust-versus-Python tests.

Run from this folder: python3 make_fixtures.py. Every value is invented; the
files only copy the layout of the public files they stand for.
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
RAW.mkdir(exist_ok=True)


def write(path: Path, lines, *, encoding="utf-8", newline="\n"):
    path.write_bytes(newline.join(lines).encode(encoding) + newline.encode(encoding))


def tsv(path: Path, rows, **kwargs):
    write(path, ["\t".join(row) for row in rows], **kwargs)


# --- HPO phenotype to genes (with ancestors, as the HPO file is) ---------------
genes = {
    "BTK": "695", "CD19": "930", "TNFRSF13B": "23495", "ADA": "100", "IL7R": "3575",
    "STAT3": "6774", "DOCK8": "81704", "FLG": "2312", "CFTR": "1080", "IL5": "3567",
    "JAK2": "3717", "WAS": "7454", "KRT14": "3861", "SCN1A": "6323", "KCNQ2": "3785",
}
terms = {
    "HP:0004313": ("Decreased antibody level in blood", ["BTK", "CD19", "TNFRSF13B"]),
    "HP:0002721": ("Immunodeficiency", ["BTK", "CD19", "TNFRSF13B", "ADA", "IL7R", "STAT3"]),
    "HP:0002715": ("Abnormality of the immune system",
                   ["BTK", "CD19", "TNFRSF13B", "ADA", "IL7R", "STAT3", "DOCK8", "FLG", "WAS", "IL5", "JAK2", "CFTR"]),
    "HP:0002090": ("Pneumonia", ["BTK", "STAT3", "DOCK8", "CFTR"]),
    "HP:0002719": ("Recurrent infections", ["BTK", "STAT3", "DOCK8", "CFTR", "ADA"]),
    "HP:0001880": ("Eosinophilia", ["DOCK8", "STAT3", "IL5"]),
    "HP:0001881": ("Abnormal leukocyte morphology", ["DOCK8", "STAT3", "IL5", "JAK2"]),
    "HP:0001871": ("Abnormality of blood and blood-forming tissues", ["DOCK8", "STAT3", "IL5", "JAK2", "ADA"]),
    "HP:0000964": ("Eczema", ["DOCK8", "STAT3", "FLG", "WAS"]),
    "HP:0011123": ('Inflammatory "abnormality" of the skin', ["FLG", "WAS", "DOCK8"]),
    "HP:0001574": ("Abnormality of the integument", ["FLG", "WAS", "DOCK8", "KRT14"]),
    "HP:0001250": ("Seizure", ["SCN1A", "KCNQ2"]),
    "HP:0002069": ("Bilateral tonic-clonic seizure", ["SCN1A"]),
    "HP:0000707": ("Abnormality of the nervous system", ["SCN1A", "KCNQ2", "DOCK8"]),
}
rows = [["hpo_id", "hpo_name", "ncbi_gene_id", "gene_symbol", "disease_id"]]
p_root_genes = []
for hp, (name, symbols) in terms.items():
    for index, symbol in enumerate(symbols):
        rows.append([hp, name, genes[symbol], symbol, f"OMIM:{600000 + index}"])
        if symbol not in p_root_genes:
            p_root_genes.append(symbol)
rows.append(["HP:0004313", "Decreased antibody level in blood", "930", "CD19", "ORPHA:1"])
rows.append(["HP:0001250", "Seizure", "", "SCN1A", "OMIM:607208"])
for index, symbol in enumerate(p_root_genes):
    rows.append(["HP:0000118", "Phenotypic abnormality", genes[symbol], symbol, f"OMIM:{610000 + index}"])
tsv(RAW / "phenotype_to_genes.txt", rows)

# --- GO annotations (GAF 2.2) --------------------------------------------------
def gaf(db_id, symbol, qualifier, go, evidence, name="protein", extra=""):
    return ["UniProtKB", db_id, symbol, qualifier, go, "PMID:1", evidence, extra, "P", name,
            symbol, "protein", "taxon:9606", "20240101", "UniProt", "", ""]

gaf_rows = [
    ["!gaf-version: 2.2"], ["!generated-by: made up for tests"],
    gaf("Q06187", "BTK", "involved_in", "GO:0002376", "IDA", "Tyrosine-protein kinase BTK"),
    gaf("Q06187", "BTK", "involved_in", "GO:0002376", "IMP"),
    gaf("P15391", "CD19", "NOT|involved_in", "GO:0002376", "IDA"),
    gaf("P40763", "STAT3", "involved_in", "GO:0002376", "IEA", 'Signal "transducer" 3'),
    gaf("Q9NZQ7", "PDL1ALIAS", "involved_in", "GO:0002376", "IMP"),
    gaf("P05231", "IL6", "acts_upstream_of", "GO:0002376", "TAS"),
    gaf("Q99999", "UNMAPPED", "involved_in", "GO:0002376", "EXP"),
    gaf("P04637", "TP53", "involved_in", "GO:0000001", "IDA"),
    ["UniProtKB", "P12345", "SHORT", "involved_in", "GO:0002376", "PMID:2", "IDA"],
]
tsv(RAW / "goa_human.gaf", gaf_rows)

# --- NCBI gene info ------------------------------------------------------------
info = [["#tax_id", "GeneID", "Symbol", "LocusTag", "Synonyms", "dbXrefs", "chromosome", "map_location",
         "description", "type_of_gene", "Symbol_from_nomenclature_authority",
         "Full_name_from_nomenclature_authority", "Nomenclature_status", "Other_designations",
         "Modification_date", "Feature_type"]]
for symbol, gid in list(genes.items()) + [("IL6", "3569"), ("DUPSYM", "900001"), ("DUPSYM", "900002"),
                                          ("OTHER", "900001")]:
    description = 'gene "with" quotes' if symbol == "BTK" else f"{symbol} protein"
    info.append(["9606", gid, symbol, "-", "-", f"HGNC:{gid}", "1", "1p1", description, "protein-coding",
                 symbol, description, "O", "-", "20240101", "-"])
info.append(["9606", "123", "SHORTROW"])
tsv(RAW / "Homo_sapiens.gene_info", info)

# --- Reactome -----------------------------------------------------------------
reactome = [
    ["695", "R-HSA-168256", "https://reactome.org/PathwayBrowser/#/R-HSA-168256", "Immune System", "TAS", "Homo sapiens"],
    ["3575", "R-HSA-1280215", "https://reactome.org/PathwayBrowser/#/R-HSA-1280215", "Cytokine Signaling in Immune system", "TAS", "Homo sapiens"],
    ["7454", "R-HSA-9999999", "https://reactome.org/PathwayBrowser/#/R-HSA-9999999", "Made-up child pathway", "IEA", "Homo sapiens"],
    ["2312", "R-HSA-1474244", "https://reactome.org/PathwayBrowser/#/R-HSA-1474244", "Extracellular matrix organization", "TAS", "Homo sapiens"],
    ["12229", "R-MMU-168256", "https://reactome.org/PathwayBrowser/#/R-MMU-168256", "Immune System", "IEA", "Mus musculus"],
    ["1080", "R-HSA-5678895"],
]
tsv(RAW / "NCBI2Reactome_All_Levels.txt", reactome)
tsv(RAW / "ReactomePathwaysRelation.txt", [
    ["R-HSA-168256", "R-HSA-1280215"], ["R-HSA-1280215", "R-HSA-9999999"], ["R-MMU-168256", "R-MMU-1280215"],
])

# --- GO DAG and gene2go ----------------------------------------------------------
write(RAW / "go-basic.obo", [
    "format-version: 1.2", "", "[Term]", "id: GO:0002376", "name: immune system process", "",
    "[Term]", "id: GO:0006952", "name: defense response", "",
    "[Term]", "id: GO:0045087", "name: innate immune response", "is_a: GO:0006952 ! defense response", "",
    "[Term]", "id: GO:0050776", "name: regulation of immune response",
    "relationship: part_of GO:0002376 ! immune system process", "",
    "[Term]", "id: GO:0000001", "name: unrelated process", "",
    "[Typedef]", "id: part_of", "name: part of",
])
tsv(RAW / "gene2go", [
    ["#tax_id", "GeneID", "GO_ID", "Evidence", "Qualifier", "GO_term", "PubMed", "Category"],
    ["9606", "6774", "GO:0002376", "IDA", "involved_in", "immune system process", "1", "Process"],
    ["9606", "3861", "GO:0045087", "IEA", "involved_in", "innate immune response", "-", "Process"],
    ["9606", "6323", "GO:0050776", "IEA", "involved_in", "regulation of immune response", "-", "Process"],
    ["10090", "3785", "GO:0002376", "IEA", "involved_in", "immune system process", "-", "Process"],
    ["9606", "3717", "GO:0000001", "IEA", "involved_in", "unrelated process", "-", "Process"],
])

# --- Protein Ontology ---------------------------------------------------------
write(RAW / "pr.obo", [
    "format-version: 1.2", "data-version: made-up", "",
    "[Term]", "id: PR:000001393", "name: interleukin-6", "xref: UniProtKB:P05231 ! IL6_HUMAN",
    "xref: NCBIGene:3569", "xref: Reactome:R-HSA-1", "",
    "[Term]", "id: PR:000000001", "name: protein", "xref: UniProtKB/Swiss-Prot:Q06187-2", "",
    "[Term]", "id: PR:000002222", "xref: UniProtKB:p99999", "xref: NCBIGene:900001", "",
    "[Term]", "id: CHEBI:1", "name: not a PR term", "xref: UniProtKB:P05231", "",
    "[Term]", "", "[Typedef]", "id: has_part", "name: has part",
])

# --- MGI ------------------------------------------------------------------------
tsv(RAW / "MGI_GenePheno.rpt", [
    ["Btk<tm1>/Btk<tm1>", "Btk<tm1>", "MGI:1", "involves: 129", "MP:0002491", "12345", "MGI:88216|MGI:99999", "MGI:2"],
    ["Cd19<tm1>/Cd19<tm1>", "Cd19<tm1>", "MGI:3", "involves: C57BL/6", "MP:0002491", "222", "MGI:88217", "MGI:4"],
    ["Cd19<tm1>/Cd19<tm1>", "Cd19<tm1>", "MGI:3", "involves: C57BL/6", "MP:0002491", "111", " MGI:88217 ", "MGI:5"],
    ["x", "y", "z", "w", "", "1", "MGI:88216", "MGI:6"],
    ["short", "row"],
])
tsv(RAW / "HOM_MouseHumanSequence.rpt", [
    ["DB Class Key", "Common Organism Name", "NCBI Taxon ID", "Symbol", "EntrezGene ID", "Mouse MGI ID",
     "HGNC ID", "OMIM Gene ID", "Genetic Location", "Genome Coordinates (mouse: GRCm39 human: GRCh38)",
     "Nucleotide RefSeq IDs", "Protein RefSeq IDs", "SWISS_PROT IDs"],
    ["1001", "mouse, laboratory", "10090", "Btk", "12229", "MGI:88216", "", "", "X", "", "", "", ""],
    ["1001", "human", "9606", "BTK", "695", "", "HGNC:1133", "300300", "X", "", "", "", ""],
    ["1002", "mouse, laboratory", "10090", "Cd19", "12478", "MGI:88217", "", "", "7", "", "", "", ""],
    ["1002", "human", "9606", "CD19", "930", "", "HGNC:1633", "107265", "16", "", "", "", ""],
    ["1002", "human", "9606", "", "", "", "", "", "", "", "", "", ""],
])

# --- UniProt --------------------------------------------------------------------
tsv(RAW / "uniprot_human_reviewed.tsv", [
    ["Entry", "Gene Names (primary)", "Gene Names", "Protein names", "Function [CC]",
     "Subcellular location [CC]", "Tissue specificity", "Involvement in disease", "GeneID"],
    ["P05231", "IL6", "IL6 IFNB2", "Interleukin-6", "FUNCTION: cytokine.", "Secreted", "", "", "3569;"],
    ["Q06187", "BTK", "BTK AGMX1", '"Weird" protein', "FUNCTION: kinase.", "", "", "DISEASE: XLA", "695;"],
    ["Q9NZQ7", "CD274", "CD274 PDL1", "PD-L1", "", "", "", "", "29126;"],
    ["P15391", "CD19", "CD19", "B-lymphocyte antigen CD19", "", "", "", "", "930;"],
    ["P40763", "STAT3", "STAT3", "Signal transducer", "", "", "", "", ";6774;6775;"],
    [],
    ["P00001", "SHORT"],
])

# --- Human Protein Atlas (one of the three files left out on purpose) -------------
tsv(RAW / "rna_tissue_consensus.tsv", [
    ["Gene", "Gene name", "Tissue", "nTPM"],
    ["ENSG00000010671", "BTK", "spleen", "120.5"],
    ["ENSG00000177455", "CD19", "lymph node", "88"],
])
tsv(RAW / "rna_single_cell_type.tsv", [
    ["Gene", "Gene name", "Cell type", "nTPM"],
    ["ENSG00000010671", "BTK", "B-cells", "300.1"],
    ["ENSG00000177455", "CD19"],
    [],
    ["ENSG00000111537", "IFNG", "T-cells", '"9.5"'],
])

# --- Orphadata product 6 (Latin-1, like the real file) ---------------------------
orpha = """<?xml version="1.0" encoding="ISO-8859-1"?>
<JDBOR date="2024-12-01 00:00:00" version="1.3.31">
  <Availability><Licence><FullName lang="en">Creative Commons Attribution 4.0</FullName></Licence></Availability>
  <DisorderList count="2">
    <Disorder id="1">
      <OrphaCode>183660</OrphaCode>
      <Name lang="en">Agammaglobulinémie liée à l'X</Name>
      <DisorderType id="21394"><Name lang="en">Disease</Name></DisorderType>
      <DisorderGeneAssociationList count="2">
        <DisorderGeneAssociation>
          <SourceOfValidation>PMID:1</SourceOfValidation>
          <Gene id="1">
            <Name lang="en">Bruton tyrosine kinase</Name>
            <Symbol>BTK</Symbol>
            <OrphaCode>119051</OrphaCode>
            <ExternalReferenceList count="3">
              <ExternalReference id="1"><Source>HGNC</Source><Reference>1133</Reference></ExternalReference>
              <ExternalReference id="2"><Source>OMIM</Source><Reference>300300</Reference></ExternalReference>
              <ExternalReference id="3"><Source>Ensembl</Source><Reference></Reference></ExternalReference>
            </ExternalReferenceList>
          </Gene>
          <DisorderGeneAssociationType id="17949"><Name lang="en">Disease-causing germline mutation(s) in</Name></DisorderGeneAssociationType>
          <DisorderGeneAssociationStatus id="17991"><Name lang="en">Assessed</Name></DisorderGeneAssociationStatus>
        </DisorderGeneAssociation>
        <DisorderGeneAssociation>
          <Gene id="2"><Name lang="en">No symbol gene</Name><Symbol></Symbol></Gene>
          <DisorderGeneAssociationStatus id="17991"><Name lang="en">Not yet assessed</Name></DisorderGeneAssociationStatus>
        </DisorderGeneAssociation>
      </DisorderGeneAssociationList>
    </Disorder>
    <Disorder id="2">
      <OrphaCode> 999999 </OrphaCode>
      <Name lang="en">Made-up &amp; mixed disorder</Name>
      <DisorderGeneAssociationList count="1">
        <DisorderGeneAssociation>
          <Gene id="3"><Symbol>DOCK8</Symbol><OrphaCode>1</OrphaCode></Gene>
          <DisorderGeneAssociationType id="1"><Name lang="en">Role in the phenotype of</Name></DisorderGeneAssociationType>
        </DisorderGeneAssociation>
      </DisorderGeneAssociationList>
    </Disorder>
  </DisorderList>
</JDBOR>"""
(RAW / "orphadata_en_product6.xml").write_bytes(orpha.encode("latin-1") + b"\n")

# --- Patient HPO table (the HPO extraction output layout) -------------------------
patients = [
    ["sample_id", "hp_id", "hp_name", "source", "original_entry", "confidence"],
    ["S1", "HP:0004313", "Decreased antibody level in blood", "rag_direct", "IgG low", "0.9"],
    ["S1", "HP:0002090", "Pneumonia", "rag_direct", '"cough\tand fever"', "0.8"],
    ["S1", "HP:0001880", "Eosinophilia", "text_match", "Eos 1.2", "1"],
    ["S1", "HP:0000006", "Autosomal dominant inheritance", "rag_direct", "AD", "1"],
    ["S1", "HP:0002719", "Recurrent infections", "rag_direct", "infections", "1"],
    ["S1", "HP:0002715", "Abnormality of the immune system", "propagated", "", "1"],
    ["S1", "HP:0004313", "Decreased antibody level in blood", "embedded_id", "HP:0004313", "1"],
    ["S2", "HP:0001250", "Seizure", "rag_direct", "seizures", "1"],
    ["S2", "HP:0000964", "Eczema", "rag_direct", '"eczema, ""severe"""', "1"],
    ["S2", "HP:0002069", "Bilateral tonic-clonic seizure", "rag_direct", "GTC", "1"],
    ["S3", "HP:0000006", "Autosomal dominant inheritance", "rag_direct", "AD", "1"],
    ["S.4", "HP:0001880", "Eosinophilia", "rag_direct", "", "1"],
    ["S.4", "HP:9999999", "Unknown", "rag_direct", "x", "1"],
    ["S5", "HP:0011123", 'Inflammatory "abnormality" of the skin', "rag_direct", "skin", "1"],
    ["S5", "HP:0002721", "Immunodeficiency", "rag_direct", "PID", "1"],
    ["S5", "HP:0000999", "obsolete Old term", "rag_direct", "old", "1"],
    ["S6", "HP:0002090", "Pneumonia", "propagated", "", "1"],
]
tsv(HERE / "patient_hpos.tsv", patients)

# --- HPO terms by name (the new Rust input) ----------------------------------------
tsv(HERE / "hpo_terms.tsv", [
    ["sample_id", "hpo_term"],
    ["T1", "Immunodeficiency"],
    ["T1", "hp:0002090"],
    ["T1", "  eczematous DERMATITIS "],
    ["T2", "HP_0001250"],
])
print("fixtures written")
