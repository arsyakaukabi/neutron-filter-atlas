"""SHA-256 of files (data manifest, frozen-input checks)."""
import hashlib


def sha256(path):
    """SHA-256 checksum of a file, as used in `data/MANIFEST.csv`.

    Parameters:
        path: file path.
    Returns:
        hexadecimal digest (str).
    Example:
        >>> from filteratlas.common.hashing import sha256
        >>> sha256("data/MANIFEST.csv")[:12]
    """
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
