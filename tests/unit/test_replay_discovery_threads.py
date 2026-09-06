import pytest

from benchmarks.scripts.replay_discovery import command_with_threads


def test_replay_discovery_thread_override_preserves_input_command() -> None:
    command = ["python", "-m", "tandemx.cli", "discover", "--threads", "1", "--reads", "reads.fa"]
    updated, baseline, replay = command_with_threads(command, 4)
    assert command[command.index("--threads") + 1] == "1"
    assert updated[updated.index("--threads") + 1] == "4"
    assert (baseline, replay) == (1, 4)


def test_replay_discovery_thread_override_rejects_ambiguous_or_invalid_values() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        command_with_threads(["discover"], 1)
    with pytest.raises(ValueError, match="Threads"):
        command_with_threads(["discover", "--threads", "1"], 0)
    with pytest.raises(ValueError, match="Baseline has an invalid"):
        command_with_threads(["discover", "--threads", "zero"], None)
