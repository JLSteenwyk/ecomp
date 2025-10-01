"""Project-wide configuration constants for the eComp tool."""

from __future__ import annotations

FORMAT_VERSION = "0.2.0"
FORMAT_VERSION_TUPLE = (0, 2, 0)
HEADER_MAGIC = b"ECOMP001"
LEGACY_HEADER_STRUCT = ">8sBBBQ"
HEADER_STRUCT = ">8sBBBQQ"
METADATA_SUFFIX = ".json"
DEFAULT_OUTPUT_FORMAT = "fasta"
SUPPORTED_INPUT_FORMATS = {
    "fasta": {"extensions": {".fasta", ".fa", ".faa", ".fna"}},
    "phylip": {"extensions": {".phy", ".phylip"}},
}


def detect_format_from_suffix(path: str) -> str | None:
    """Return an alignment format inferred from *path* or ``None`` if unknown."""

    import pathlib

    suffix = pathlib.Path(path).suffix.lower()
    if not suffix:
        return None
    for fmt, data in SUPPORTED_INPUT_FORMATS.items():
        if suffix in data["extensions"]:
            return fmt
    return None
