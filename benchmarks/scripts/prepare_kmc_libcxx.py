"""Small recorded portability patch for upstream KMC with Apple libc++.

Replaces GNU's copy_n header/namespace with the C++11 standard equivalent.
No counting or filtering algorithm is changed. Apply only to an isolated clone.
"""
import argparse
from pathlib import Path
import subprocess


def apply(source: Path, patch_path: Path) -> None:
    names = ["kmc_core/defs.h", "kmc_api/kmer_defs.h"]
    contents = [(source / name).read_text() for name in names]
    if not all("#include <ext/algorithm>" in value for value in contents):
        raise ValueError("Expected upstream GNU headers not found; inspect source before patching")
    for name, value in zip(names, contents):
        (source / name).write_text(value.replace("#include <ext/algorithm>", "#include <algorithm>")
                                 .replace("using __gnu_cxx::copy_n;", "using std::copy_n;"))
    patch_path.write_bytes(subprocess.check_output(["git", "diff", "--", *names], cwd=source))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--patch", type=Path, required=True)
    arguments = parser.parse_args()
    apply(arguments.source, arguments.patch)
