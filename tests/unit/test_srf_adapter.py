import pytest

from benchmarks.scripts.run_srf_pilot import parse_bed, workflow


def test_native_bed_filtering_keeps_coordinates_and_motif_length(tmp_path):
    path = tmp_path / "srf.bed"
    path.write_text("r\t10\t90\tsrf#circ1-20\t0.9\t100\t20\t1\n"
                    "r\t90\t100\tsrf#circ1-20\t0.9\t100\t20\t0\n")
    motifs = {"srf#circ1-20": "ACGT" * 5}
    result = parse_bed(path, motifs)
    assert [(r.start, r.end, r.period) for r in result] == [(10, 90, 20)]
    assert len(parse_bed(path, motifs, keep_only=False)) == 2
    with pytest.raises(ValueError, match="length"):
        parse_bed(path, {"srf#circ1-20": "ACGT"})
    path.write_text("")
    assert parse_bed(path, motifs) == []
    path.write_text("unknown\tformat\n")
    with pytest.raises(ValueError, match="Malformed"):
        parse_bed(path, motifs)


def test_empty_successful_kmer_dump_skips_native_graph_without_fake_fasta(tmp_path, monkeypatch):
    from pathlib import Path
    import benchmarks.scripts.run_srf_pilot as runner

    called = []

    def fake_process(command, stdout, stderr, timeout):
        called.append(command[0])
        if command[0] == "dump":
            Path(command[-1]).write_text("")
        elif command[0] != "kmc":
            pytest.fail("Empty dump must not be sent to native graph assembly")
        return dict(exit_code=0, runtime_seconds=0.1, peak_rss_mib=1, timed_out=False)

    monkeypatch.setattr(runner, "run_process", fake_process)
    reads = tmp_path / "reads.fa"
    reads.write_text(">test\nACGT\n")
    tools = {key: Path(key) for key in ("kmc", "dump", "srf", "k8", "minimap2", "utils")}
    outdir = tmp_path / "output"
    record = workflow(reads, outdir, tools, 100, 5)
    assert called == ["kmc", "dump"]
    assert record["status"] == "no_eligible_kmers"
    assert record["skipped_stages"][0] == "assemble"
    assert not (outdir / "srf.fa").exists()
