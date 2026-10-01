import random
import sys
from pathlib import Path

import pytest

TESTS = Path(__file__).parent
sys.path.insert(0, str(TESTS))


@pytest.fixture(scope="session")
def fake_refs(tmp_path_factory) -> Path:
    """Random stand-ins for the CARP rDNA library and the organelle genomes."""
    d = tmp_path_factory.mktemp("refs")
    r = random.Random(9)

    def s(n):
        return "".join(r.choices("ACGT", k=n))

    lib = [("18S", 1800), ("18S", 1700), ("5.8S", 160), ("25S", 3400), ("5S", 120)]
    (d / "rdna_library.fasta").write_text(
        "".join(f">u{i}#rDNA/45S_rDNA/{c}\n{s(n)}\n" for i, (c, n) in enumerate(lib)))
    (d / "plastid_NC_014057.1.fa").write_text(">NC_014057.1\n" + s(122_169) + "\n")
    (d / "mito_PP555264.1.fa").write_text(">PP555264.1\n" + s(50_000) + "\n")
    return d


@pytest.fixture(scope="session")
def testdata_noreads(fake_refs, tmp_path_factory) -> Path:
    """Synthetic data set (SPEC §11.1) without reads."""
    import make_testdata

    out = tmp_path_factory.mktemp("testdata")
    make_testdata.main([str(out), "--refs", str(fake_refs), "--no-reads"])
    return out
