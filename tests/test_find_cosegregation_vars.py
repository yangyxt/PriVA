from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from find_cosegregation_vars import (  # noqa: E402
    check_variant_presence,
    find_cosegregating_variants,
    parse_gt,
)


def test_pysam_integer_genotypes_preserve_explicit_x_allele_states() -> None:
    assert parse_gt((0, 1)) == [0, 1]
    assert parse_gt((None, None)) is None
    assert check_variant_presence([1], 1, "chrX") == "Hemizygous"
    assert check_variant_presence([1, 1], 1, "chrX") == "Hemizygous"
    assert check_variant_presence([0, 1], 2, "chrX") == "Heterozygous"
    assert check_variant_presence([1, 1], 2, "chrX") == "Homozygous"
    assert check_variant_presence([0, 1], 3, "chrX") is None


def test_cosegregation_counts_known_x_sexes_and_excludes_unknown(
    tmp_path: Path,
) -> None:
    pedigree = tmp_path / "family.ped"
    pedigree.write_text(
        "F1\tMALE\t0\t0\t1\t2\n"
        "F1\tFEMALE\t0\t0\t2\t2\n"
        "F1\tUNKNOWN\t0\t0\t3\t2\n",
        encoding="utf-8",
    )
    vcf = tmp_path / "family.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=X>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        "MALE\tFEMALE\tUNKNOWN\n"
        "X\t100\t.\tA\tG\t.\tPASS\t.\tGT\t1\t0/1\t0/1\n",
        encoding="utf-8",
    )

    result = find_cosegregating_variants(vcf, pedigree, mode="both")

    assert result[("X", 100, "A", "G")] == {
        "Affected_segregated_inds": 2,
        "Unaffected_segregated_inds": 0,
        "Affected_homozygous_inds": 0,
        "Unaffected_nonhomo_inds": 0,
        "Male_affected_segregated_inds": 1,
        "Male_unaffected_segregated_inds": 0,
        "Female_affected_segregated_inds": 1,
        "Female_unaffected_segregated_inds": 0,
    }
