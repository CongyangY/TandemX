"""Static, source-backed figures for the family report.

The renderer deliberately treats abundance and assembly representation as
estimates.  Missing values are omitted from quantitative marks and retained in
the source tables; no threshold or biological state is inferred here.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import textwrap
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch

_BLUE = "#2f638f"
_GOLD = "#c18d32"
_ORANGE = "#c46d3b"
_INK = "#263238"
_MUTED = "#64727a"
_GRID = "#d9e0e3"
_UNRESOLVED = "#a5afb4"
_FIG_NAMES = ("summary", "family_abundance_vs_assembly", "top_underrepresented", "family_landscape", "family_hierarchy")


def _num(row: dict[str, Any], key: str) -> float | None:
    value = row.get(key)
    if value in (None, "", "NA", "NaN", "--", "unresolved"):
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _families(data: dict[str, Any]) -> list[dict[str, Any]]:
    rows = data.get("families", [])
    return [dict(row) for row in rows if isinstance(row, dict)]


def _label(value: float | int, unit: str = "bp") -> str:
    value = float(value)
    if abs(value) >= 1e9:
        return f"{value / 1e9:.3g} Gb"
    if abs(value) >= 1e6:
        return f"{value / 1e6:.3g} Mb"
    if abs(value) >= 1e3:
        return f"{value / 1e3:.3g} kb"
    return f"{value:.0f} {unit}"


def _style(ax: plt.Axes, title: str, subtitle: str = "") -> None:
    ax.set_title(title, loc="left", fontsize=11, color=_INK, fontweight="bold", pad=12)
    if subtitle:
        ax.text(0, 1.01, subtitle, transform=ax.transAxes, fontsize=8, color=_MUTED, va="bottom")
    ax.grid(axis="both", color=_GRID, linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#89979d")
    ax.tick_params(colors=_MUTED, labelsize=8)


def _empty(ax: plt.Axes, message: str = "No quantitative data available") -> None:
    ax.text(0.5, 0.5, message, transform=ax.transAxes, ha="center", va="center", color=_MUTED, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])


def _write_source(path: Path, rows: Iterable[dict[str, Any]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join("" if row.get(col) is None else str(row.get(col)) for col in columns) + "\n")


def _save(fig: plt.Figure, name: str, out: Path, source_rows: list[dict[str, Any]], columns: list[str]) -> dict[str, str]:
    out.mkdir(parents=True, exist_ok=True)
    svg, pdf, png = out / f"{name}.svg", out / f"{name}.pdf", out / f"{name}.png"
    fig.savefig(svg, format="svg", facecolor="white", bbox_inches="tight")
    fig.savefig(pdf, format="pdf", facecolor="white", bbox_inches="tight")
    fig.savefig(png, format="png", dpi=180, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    source = out / f"{name}.tsv"
    _write_source(source, source_rows, columns)
    svg_text = svg.read_text(encoding="utf-8", errors="replace")
    receipt = out / f"{name}.receipt.json"
    files = {"svg": svg, "pdf": pdf, "png": png, "source_tsv": source}
    payload: dict[str, Any] = {"figure": name, "source_tsv": source.name, "files": {}}
    for kind, path in files.items():
        payload["files"][kind] = {"path": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size}
    payload["svg_node_counts"] = {"text": len(re.findall(r"<text(?: |>)", svg_text)), "raster": len(re.findall(r"<image(?: |>)", svg_text))}
    receipt.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"svg": str(svg), "pdf": str(pdf), "png": str(png), "source_tsv": str(source), "receipt": str(receipt)}


def _scatter(rows: list[dict[str, Any]], out: Path) -> dict[str, str]:
    fig, ax = plt.subplots(figsize=(6.8, 5.1), constrained_layout=True)
    _style(ax, "Read abundance and assembly representation", "Estimated family-level bp; equality line is a visual reference")
    points = []
    for row in rows:
        x, y = _num(row, "estimated_abundance_bp"), _num(row, "assembly_representation_bp")
        if x is not None and y is not None and x >= 0 and y >= 0:
            points.append((row, x, y))
    if not points:
        _empty(ax)
    else:
        maxv = max(max(x, y) for _, x, y in points) or 1
        positive = [v for _, x, y in points for v in (x, y) if v > 0]
        linthresh = max(min(positive) * 0.5, maxv * 1e-5) if positive else 1
        ax.set_xscale("symlog", linthresh=linthresh); ax.set_yscale("symlog", linthresh=linthresh)
        for row, x, y in points:
            confidence = str(row.get("confidence") or "")
            color = _UNRESOLVED if confidence.lower() in {"unresolved", "unknown", ""} else _BLUE
            marker = "o" if color == _BLUE else "o"
            ax.scatter(x, y, s=32, c=color, edgecolors="white", linewidths=0.7, alpha=0.9, marker=marker, zorder=3)
        bound = maxv * 1.15
        ax.plot([0, bound], [0, bound], color=_INK, linestyle="--", linewidth=0.9, label="equality")
        ax.set_xlim(left=0); ax.set_ylim(bottom=0)
        ax.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=_BLUE, markersize=6, label="confidence reported"), Line2D([0], [0], marker="o", color="w", markerfacecolor=_UNRESOLVED, markersize=6, label="unresolved / absent"), Line2D([0], [0], color=_INK, linestyle="--", label="equality")], frameon=False, fontsize=8, loc="upper left")
        ax.set_xlabel("Estimated read abundance (bp)"); ax.set_ylabel("Assembly representation (bp)")
    source = [{"family_id": r.get("family_id"), "estimated_abundance_bp": _num(r, "estimated_abundance_bp"), "assembly_representation_bp": _num(r, "assembly_representation_bp"), "assembly_read_ratio": _num(r, "assembly_read_ratio"), "confidence": r.get("confidence"), "warning": r.get("warning")} for r in rows]
    return _save(fig, "family_abundance_vs_assembly", out, source, list(source[0]) if source else ["family_id", "estimated_abundance_bp", "assembly_representation_bp", "assembly_read_ratio", "confidence", "warning"])


def _top_deficit(rows: list[dict[str, Any]], out: Path) -> dict[str, str]:
    fig, ax = plt.subplots(figsize=(7.4, 5.2), constrained_layout=True)
    _style(ax, "Families with the largest abundance deficit", "Sorted by estimated read abundance minus assembly representation; values are estimates")
    valid = [(r, _num(r, "abundance_deficit_bp")) for r in rows]
    valid = [(r, d) for r, d in valid if d is not None and d > 0]
    valid.sort(key=lambda x: x[1], reverse=True); valid = valid[:10]
    if not valid: _empty(ax)
    else:
        labels = [str(r.get("family_id", "")) for r, _ in valid][::-1]
        deficits = [d for _, d in valid][::-1]
        read = [(_num(r, "estimated_abundance_bp") or 0) for r, _ in valid][::-1]
        assembly = [(_num(r, "assembly_representation_bp") or 0) for r, _ in valid][::-1]
        y = list(range(len(labels))); height = .34
        ax.barh([v - height / 2 for v in y], read, height=height, color=_BLUE, label="read abundance")
        ax.barh([v + height / 2 for v in y], assembly, height=height, color=_GOLD, label="assembly representation")
        ax.axvline(0, color=_INK, linewidth=0.8)
        ax.set_yticks(y, labels); ax.set_xlabel("Estimated sequence amount (bp)")
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: _label(x, "bp")))
        ax.legend(frameon=False, fontsize=8, loc="lower right")
        maxv = max(read + assembly + [1])
        for ypos, (r, d) in enumerate(valid[::-1]):
            conf = str(r.get("confidence") or "unresolved")
            ax.text(maxv * 1.01, ypos, f"deficit {_label(d)} · {conf}", va="center", fontsize=7, color=_MUTED, clip_on=False)
        ax.set_xlim(0, maxv * 1.32)
    source = [{"family_id": r.get("family_id"), "estimated_abundance_bp": _num(r, "estimated_abundance_bp"), "assembly_representation_bp": _num(r, "assembly_representation_bp"), "abundance_deficit_bp": d, "confidence": r.get("confidence"), "warning": r.get("warning")} for r, d in valid]
    return _save(fig, "top_underrepresented", out, source, list(source[0]) if source else ["family_id", "estimated_abundance_bp", "assembly_representation_bp", "abundance_deficit_bp", "confidence", "warning"])


def _landscape(rows: list[dict[str, Any]], out: Path) -> dict[str, str]:
    fig, ax = plt.subplots(figsize=(7.0, 5.2), constrained_layout=True)
    use_abundance = any(_num(r, "estimated_abundance_bp") is not None for r in rows)
    ylabel = "Estimated abundance (bp)" if use_abundance else "Supporting reads (count)"
    _style(ax, "Family landscape", f"Monomer length versus {ylabel.lower()}; colour encodes GC fraction")
    value_key = "estimated_abundance_bp" if use_abundance else "support_read_count"
    valid = [(r, _num(r, "monomer_length_bp"), _num(r, value_key)) for r in rows]
    valid = [(r, l, a) for r, l, a in valid if l is not None and a is not None and l > 0 and a >= 0]
    if not valid: _empty(ax)
    else:
        gc = [_num(r, "gc_fraction") for r, _, _ in valid]; c = [v if v is not None else None for v in gc]
        cmap = matplotlib.colormaps["cividis"]
        colors = [cmap(v) if v is not None and 0 <= v <= 1 else "#a5afb4" for v in c]
        ax.scatter([l for _, l, _ in valid], [a for _, _, a in valid], color=colors, s=34, edgecolors="white", linewidths=.6)
        ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=max(min([a for _, _, a in valid if a > 0], default=1) * .5, 1))
        ax.set_xlabel("Monomer length (bp)"); ax.set_ylabel(ylabel)
        # Discrete legend keeps the exported SVG vector-only (Matplotlib colorbars
        # can embed a raster gradient) while preserving the GC encoding.
        ax.legend(handles=[Patch(facecolor=cmap(v), edgecolor="none", label=lab) for v, lab in ((.15, "low GC"), (.5, "mid GC"), (.85, "high GC"))], title="GC fraction", frameon=False, fontsize=7, title_fontsize=8, loc="best")
    source = [{"family_id": r.get("family_id"), "monomer_length_bp": l, "landscape_measure": value_key, "landscape_value": a, "estimated_abundance_bp": _num(r, "estimated_abundance_bp"), "gc_fraction": _num(r, "gc_fraction"), "support_read_count": _num(r, "support_read_count"), "warning": r.get("warning")} for r, l, a in valid]
    return _save(fig, "family_landscape", out, source, list(source[0]) if source else ["family_id", "monomer_length_bp", "estimated_abundance_bp", "gc_fraction", "support_read_count", "warning"])


def _hierarchy(data: dict[str, Any], rows: list[dict[str, Any]], out: Path) -> dict[str, str]:
    fig, ax = plt.subplots(figsize=(8.0, 5.2), constrained_layout=False)
    _style(ax, "Candidate family architecture relationships", "Sequence-supported heuristic links; candidate period multiples are not validated HOR structures")
    edges = [e for e in data.get("architecture_edges", []) if isinstance(e, dict)]
    ranked = {str(r.get("family_id")): i for i, r in enumerate(sorted(rows, key=lambda r: (_num(r, "estimated_abundance_bp") or 0), reverse=True))}
    edges.sort(key=lambda e: (ranked.get(str(e.get("shorter_family_id")), 10**9), ranked.get(str(e.get("longer_family_id")), 10**9)))
    total_edges = len(edges)
    edges = edges[:6]
    nodes = sorted({str(e.get(k)) for e in edges for k in ("shorter_family_id", "longer_family_id") if e.get(k)}, key=lambda n: ranked.get(n, 10**9))[:8]
    edges = [e for e in edges if str(e.get("shorter_family_id")) in nodes and str(e.get("longer_family_id")) in nodes]
    ax.set_title(f"Candidate family architecture relationships ({len(edges)}/{total_edges} edges shown)", loc="left", fontsize=11, color=_INK, fontweight="bold", pad=12)
    if not edges or not nodes: _empty(ax, "No candidate relationship edges available")
    else:
        lengths = {str(r.get("family_id")): (_num(r, "monomer_length_bp") or 0) for r in rows}
        layers = sorted(nodes, key=lambda n: (lengths.get(n, 0), n))
        pos = {n: (0.15 if i < len(layers) / 2 else 0.85, .12 + .76 * (i % max(1, math.ceil(len(layers) / 2))) / max(1, math.ceil(len(layers) / 2) - 1)) for i, n in enumerate(layers)}
        for e in edges:
            a, b = str(e.get("shorter_family_id")), str(e.get("longer_family_id")); x1, y1 = pos[a]; x2, y2 = pos[b]
            color = _GOLD if str(e.get("edge_type")) == "putative_period_multiple" else _UNRESOLVED
            ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops={"arrowstyle": "->", "color": color, "lw": 1.3, "alpha": .8})
            ratio = _num(e, "length_ratio"); ident = _num(e, "local_identity")
            edge_label = f"r={ratio:.2f}" if ratio is not None else "r=NA"
            edge_label += f", id={ident:.2f}" if ident is not None else ", id=NA"
            ax.text((x1 + x2) / 2, (y1 + y2) / 2 + .035, edge_label, fontsize=6.5, color=color, ha="center", va="center")
        for n, (x, y) in pos.items():
            ax.scatter([x], [y], s=480, facecolor="white", edgecolor=_BLUE, linewidth=1.2, zorder=3)
            ax.text(x, y, n, ha="center", va="center", fontsize=7, color=_INK, zorder=4)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
        ax.legend(handles=[Line2D([0], [0], color=_GOLD, lw=1.5, label="putative period multiple"), Line2D([0], [0], color=_UNRESOLVED, lw=1.5, label="unresolved related / partial")], frameon=False, fontsize=8, loc="lower left")
    fig.subplots_adjust(left=.02, right=.98, top=.84, bottom=.10)
    source = [dict(e) for e in edges]
    cols = list(source[0]) if source else ["shorter_family_id", "longer_family_id", "edge_type", "status", "warning"]
    return _save(fig, "family_hierarchy", out, source, cols)


def _evidence_cards(rows: list[dict[str, Any]], out: Path) -> dict[str, str]:
    """Render a compact, editable evidence card for each of the six largest families."""
    ranked = sorted(rows, key=lambda r: (_num(r, "abundance_deficit_bp") or 0, _num(r, "estimated_abundance_bp") or 0), reverse=True)[:6]
    # A fixed one-column grid prevents long metadata strings from collapsing a
    # constrained-layout two-column grid on large real catalogues.
    fig, axes = plt.subplots(max(1, len(ranked)), 1, figsize=(10, max(3.3, 2.05 * max(1, len(ranked)))), squeeze=False, constrained_layout=False)
    fig.suptitle("Family-level evidence cards", x=.03, ha="left", fontsize=14, fontweight="bold", color=_INK)
    fig.subplots_adjust(left=.035, right=.965, top=.92, bottom=.04, hspace=.28)
    for index, ax in enumerate(axes[:, 0]):
        ax.axis("off")
        if index >= len(ranked):
            continue
        row = ranked[index]
        card = FancyBboxPatch((.01, .04), .98, .9, boxstyle="round,pad=.012,rounding_size=.018", facecolor="#f7fafb", edgecolor="#cbd7dc", linewidth=1)
        ax.add_patch(card)
        family = str(row.get("family_id") or "unresolved family")
        confidence = str(row.get("confidence") or "unresolved")
        ax.text(.05, .82, family, fontsize=11, fontweight="bold", color=_INK)
        ax.text(.95, .82, confidence, fontsize=8, color=_BLUE if confidence.lower() != "unresolved" else _MUTED, ha="right")
        abundance, assembly, ratio = (_num(row, key) for key in ("estimated_abundance_bp", "assembly_representation_bp", "assembly_read_ratio"))
        deficit = _num(row, "abundance_deficit_bp")
        lines = [f"read abundance: {_label(abundance) if abundance is not None else 'unresolved'}", f"assembly representation: {_label(assembly) if assembly is not None else 'unresolved'}", f"assembly/read ratio: {ratio:.4g}" if ratio is not None else "assembly/read ratio: unresolved", f"abundance deficit: {_label(deficit) if deficit is not None else 'unresolved'}"]
        cross = row.get("cross_k_values", row.get("cross_k"))
        if cross not in (None, "", [], {}):
            text = json.dumps(cross, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
            lines.append("cross-k evidence: " + (text[:100] + "…" if len(text) > 101 else text))
        warning = str(row.get("warning") or "")
        if warning:
            warning_codes = [code for code in warning.split(";") if code]
            readable = {
                "cluster_identity": "Cluster similarity threshold",
                "genome_background_uniqueness_not_verified": "Genome uniqueness not assessed",
                "haploid_depth_estimated_from_total_read_bases_and_genome_size": "Depth estimated from total reads",
                "assembly_missing": "Assembly evidence unavailable",
                "zero_estimate": "Zero estimate reported",
            }
            shown_codes = [readable.get(code, code.replace("_", " ").capitalize()) for code in warning_codes[:3]]
            shown = ";".join(shown_codes)
            remaining = len(warning_codes) - min(3, len(warning_codes))
            suffix = f" ({remaining} further warnings in table)" if remaining else ""
            lines.extend(textwrap.wrap("Notes: " + shown + suffix, width=92, subsequent_indent="       "))
        ax.text(.06, .68, "\n".join(lines), va="top", fontsize=8.2, color=_INK, linespacing=1.28)
    source = []
    for row in ranked:
        source.append({"family_id": row.get("family_id"), "estimated_abundance_bp": _num(row, "estimated_abundance_bp"), "assembly_representation_bp": _num(row, "assembly_representation_bp"), "assembly_read_ratio": _num(row, "assembly_read_ratio"), "abundance_deficit_bp": _num(row, "abundance_deficit_bp"), "confidence": row.get("confidence"), "warning": row.get("warning"), "cross_k_values": json.dumps(row.get("cross_k_values", row.get("cross_k")), ensure_ascii=False, sort_keys=True) if row.get("cross_k_values", row.get("cross_k")) not in (None, "", [], {}) else None})
    columns = list(source[0]) if source else ["family_id", "estimated_abundance_bp", "assembly_representation_bp", "assembly_read_ratio", "abundance_deficit_bp", "confidence", "warning", "cross_k_values"]
    return _save(fig, "family_evidence_cards", out, source, columns)


def _summary(rows: list[dict[str, Any]], out: Path, data: dict[str, Any]) -> dict[str, str]:
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.2), constrained_layout=True)
    fig.suptitle("TandemX family evidence summary", x=.03, ha="left", fontsize=14, fontweight="bold", color=_INK)
    ax = axes[0, 0]; _style(ax, "Family evidence by confidence")
    counts: dict[str, int] = {}
    for r in rows: counts[str(r.get("confidence") or "unresolved")] = counts.get(str(r.get("confidence") or "unresolved"), 0) + 1
    if counts: ax.bar(list(counts), list(counts.values()), color=[_BLUE if k.lower() != "unresolved" else _UNRESOLVED for k in counts]); ax.set_ylabel("Families"); ax.tick_params(axis="x", rotation=25)
    else: _empty(ax)
    ax = axes[0, 1]; _style(ax, "Abundance versus assembly representation")
    valid = [(r, _num(r, "estimated_abundance_bp"), _num(r, "assembly_representation_bp")) for r in rows]; valid = [(r,x,y) for r,x,y in valid if x is not None and y is not None and x >= 0 and y >= 0]
    if valid:
        ax.scatter([x for _,x,_ in valid], [y for _,_,y in valid], s=18, color=_BLUE, alpha=.75); ax.set_xscale("symlog"); ax.set_yscale("symlog"); ax.set_xlabel("Read abundance (bp)"); ax.set_ylabel("Assembly (bp)")
    else: _empty(ax)
    ax = axes[1, 0]; _style(ax, "Largest abundance deficits")
    deficits = sorted([(str(r.get("family_id")), _num(r, "abundance_deficit_bp")) for r in rows if _num(r, "abundance_deficit_bp") is not None], key=lambda x:x[1], reverse=True)[:6]
    if deficits: ax.barh([x[0] for x in deficits][::-1], [x[1] for x in deficits][::-1], color=_ORANGE); ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x,_: _label(x)))
    else: _empty(ax)
    ax = axes[1, 1]; _style(ax, "Monomer length distribution")
    lengths = [_num(r, "monomer_length_bp") for r in rows]; lengths = [x for x in lengths if x is not None and x > 0]
    if lengths: ax.hist(lengths, bins=min(12, max(3, len(set(lengths)))), color=_GOLD, edgecolor="white"); ax.set_xlabel("Monomer length (bp)"); ax.set_ylabel("Families")
    else: _empty(ax)
    source = [{"family_id": r.get("family_id"), "confidence": r.get("confidence"), "estimated_abundance_bp": _num(r,"estimated_abundance_bp"), "assembly_representation_bp": _num(r,"assembly_representation_bp"), "abundance_deficit_bp": _num(r,"abundance_deficit_bp"), "monomer_length_bp": _num(r,"monomer_length_bp"), "warning": r.get("warning")} for r in rows]
    return _save(fig, "summary", out, source, list(source[0]) if source else ["family_id","confidence","estimated_abundance_bp","assembly_representation_bp","abundance_deficit_bp","monomer_length_bp","warning"])


def write_report_figures(outdir: Path, data: dict[str, Any]) -> dict[str, str]:
    """Write five static figure families and auditable source/receipt files.

    Returns a flat mapping of ``<figure>_<format>`` to absolute paths.  The
    ``*_svg`` keys are suitable for inlining into the companion HTML report.
    """
    outdir = Path(outdir)
    out = outdir / "figures"
    rows = _families(data)
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "pdf.fonttype": 42, "axes.labelcolor": _INK, "figure.facecolor": "white"})
    results = {}
    generated = {"summary": _summary(rows, out, data), "family_abundance_vs_assembly": _scatter(rows, out), "top_underrepresented": _top_deficit(rows, out), "family_landscape": _landscape(rows, out), "family_hierarchy": _hierarchy(data, rows, out), "family_evidence_cards": _evidence_cards(rows, out)}
    for name, files in generated.items():
        for kind, path in files.items(): results[f"{name}_{kind}"] = path
    manifest = out / "figure_manifest.json"
    manifest.write_text(json.dumps({"schema_version": data.get("schema_version"), "figures": generated}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    results["figure_manifest"] = str(manifest)
    return results
