# acmg_auto
The repository storing all the scripts and data needed for running acmg based variant prioritization pipeline

Breakdown of annotation resources used for variant deleterious effect evaluation
1. splicing disruption effect: SpliceAI (VEP_plugin + standalone client) || Pangolin (standalone client, deprecated since Pangolin is not easy to deploy)
2. 5-UTR uORF disruption effect: UTRannotator (VEP_plugin)
3. General Deleterious effect prediction: CADD (standalone client) || PrimateAI (VEP_plugin) || AlphaMissense (VEP_plugin)
4. Haplo-insufficiency: LOEUF (VEP_plugin)
5. Transcript disruption effect prediction: NMD (VEP_plugin)
6. Clinical variants: ClinVar (bcftools annotate)
7. Population-wise allele frequency + number of homozygous carriers: gnomADv4 (bcftools annotate)
8. Conservation: Conservation (VEP_plugin)
9. Codon based evaluation: SameCodon (VEP_plugin, only used with internet connection in database mode, needs separate running)
