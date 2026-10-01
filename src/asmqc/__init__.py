"""asmqc: comparable QC of chromosome-level pea genome assemblies."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("asmqc")
except PackageNotFoundError:
    __version__ = "0+unknown"
