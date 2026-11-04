from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from acmg_consequence import get_variant_type  # noqa: E402


@pytest.mark.parametrize(
    "hgvsp, expected",
    [
        ("ENSP00000439902.1:p.Asn3264LeufsTer12", "frameshift"),
        ("p.Arg123GlyfsTer5", "frameshift"),
        ("p.Arg123Glyfs*5", "frameshift"),
        ("p.R123Gfs*5", "frameshift"),
        ("p.Gly123fs", "frameshift"),
        ("p.Arg123fsTer5", "frameshift"),
        ("p.Arg123Ter", "nonsense"),
        ("p.Arg123*", "nonsense"),
        ("p.Ala123delinsGly", "delins"),
        ("p.Ala123del", "deletion"),
        ("p.Ala123dup", "insertion_type"),
        ("p.Ala123_Gly124insTrp", "insertion_type"),
        ("p.Ala123Gly", "missense"),
        ("p.Met1?", None),
        ("p.Ala123=", None),
        (float("nan"), None),
    ],
)
def test_variant_type(hgvsp, expected):
    assert get_variant_type(hgvsp) == expected
