from benchmarks.scripts.run_unitfinder_interface_smoke import build_stage_manifest


def test_unitfinder_smoke_stage_manifest_keeps_network_off_and_input_read_only(
    tmp_path,
) -> None:
    config = {
        "image_tag": "tandemx/unitfinder:test",
        "docker": {"platform": "linux/amd64", "network": "none"},
        "command_arguments": ["--seq", "/input/test.fa", "--name", "test", "--flag", "1"],
    }
    input_dir = tmp_path / "input"
    run_dir = tmp_path / "run"
    input_dir.mkdir()
    run_dir.mkdir()
    manifest = build_stage_manifest(config, input_dir, run_dir, "docker")
    assert len(manifest["stages"]) == 2
    for stage in manifest["stages"]:
        command = stage["command"]
        assert command[command.index("--network") + 1] == "none"
        assert f"{input_dir.resolve()}:/input:ro" in command
    smoke = manifest["stages"][1]["command"]
    assert "/opt/unitFinder/bin/unitFinder.py" in smoke
    assert smoke[-6:] == config["command_arguments"]
