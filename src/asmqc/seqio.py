"""Minimal sequence-file helpers shared by validation, prep and modules."""

import gzip
import hashlib
from collections.abc import Iterator
from pathlib import Path
from typing import BinaryIO

GZIP_MAGIC = b"\x1f\x8b"


def is_gzip(path: Path) -> bool:
    with path.open("rb") as fh:
        return fh.read(2) == GZIP_MAGIC


def open_bytes(path: Path) -> BinaryIO:
    """Open a plain or gzip-compressed file for binary reading."""
    return gzip.open(path, "rb") if is_gzip(path) else path.open("rb")


def fasta_names(path: Path) -> Iterator[str]:
    """Yield sequence names (first word after '>') without loading sequences."""
    with open_bytes(path) as fh:
        for line in fh:
            if line.startswith(b">"):
                yield line[1:].split(None, 1)[0].decode() if line[1:].strip() else ""


def file_md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()
