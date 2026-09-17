"""Check that structural controls preserve length and change tile identity."""

from __future__ import annotations

import pytest

from benchmarks.scripts.build_native_equal_length_challenge import OPERATIONS, transform


def test_equal_length_challenge_operations_are_non_silent() -> None:
    tiles = [f"ACG{'ACGT'[index % 4]}TT{'ACGT'[(index // 4) % 4]}{'ACGT'[index // 16]}"
             for index in range(20)]
    for operation in OPERATIONS:
        result = transform(tiles, operation)
        assert len(result) == len(tiles)
        assert all(len(value) == 8 for value in result)
        assert (result == tiles) is (operation == "intact")
    assert transform(tiles, "swap_adjacent_05_06")[5:7] == [tiles[6], tiles[5]]
    assert transform(tiles, "replace_tile_05_with_06")[5] == tiles[6]
    assert tiles[5] == "ACGCTTCA"  # Source list was not edited in place.
    with pytest.raises(ValueError, match="Unknown"):
        transform(tiles, "undeclared")
