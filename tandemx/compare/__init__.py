"""Read-vs-assembly comparison interfaces for TandemX."""

from tandemx.compare.mvp import (
    AssemblyReadComparison,
    CompareConfig,
    classify_assembly_read_ratio,
    compare_assembly_to_reads,
    compare_toy_abundance,
    read_arrays_bed,
    read_copy_number,
    write_comparisons,
)

__all__ = [
    "AssemblyReadComparison",
    "CompareConfig",
    "classify_assembly_read_ratio",
    "compare_assembly_to_reads",
    "compare_toy_abundance",
    "read_arrays_bed",
    "read_copy_number",
    "write_comparisons",
]
