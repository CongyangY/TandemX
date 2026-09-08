from pathlib import Path

from benchmarks.unitfinder.source_audit import scan_python_sources


def test_unitfinder_static_scan_retains_syntax_and_shell_failures(tmp_path: Path) -> None:
    (tmp_path / "good.py").write_text(
        "import numpy\nfrom Bio import SeqIO\nimport os\nos.system('trf input.fa')\n"
    )
    (tmp_path / "bad.py").write_text("value = (\n")
    result = scan_python_sources(tmp_path, ["good.py", "bad.py"])
    assert result["python_file_count"] == 2
    assert result["syntax_failure_count"] == 1
    assert result["syntax_failures"][0]["file"] == "bad.py"
    assert result["os_system_call_count"] == 1
    assert {"numpy", "Bio", "os"} <= set(result["import_roots"])
