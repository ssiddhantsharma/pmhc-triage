"""Rich terminal rendering of scored targets (presentation only).

Kept separate so the package imports without rich; the CLI falls back to plain text if
rich is absent. render: one panel per target. render_rank: one table across all
(target, population) sorted by effective-N, with missing rows last (never dropped).
"""

from __future__ import annotations

from rich.console import Console, Group
from rich.panel import Panel
from rich.table import Table
from rich.text import Text


def _fmt_n(v) -> str:
    return f"{v:,.0f}" if isinstance(v, (int, float)) else "n/a"


def _fmt_frac(v) -> str:
    return f"{v:.3f}" if isinstance(v, (int, float)) else "n/a"


def _fmt_ci(ci) -> str:
    if not ci or ci[0] is None:
        return "-"
    return f"{ci[0]:,.0f}-{ci[1]:,.0f}"


def _effn_text(row) -> Text:
    """Effective-N cell: green when computed, dim 'n/a' when a factor is missing."""
    n = row.get("effective_n")
    if n is None:
        return Text("n/a", style="dim red")
    return Text(_fmt_n(n), style="bold green")


def _target_panel(ts) -> Panel:
    t = Table(box=None, pad_edge=False, expand=False)
    for col, just in [("population", "left"), ("effective-N", "right"),
                      ("MC 95%", "right"), ("incidence", "right"),
                      ("antigen", "right"), ("coverage", "right"), ("note", "left")]:
        t.add_column(col, justify=just)
    for row in ts.rows():
        mc = row.get("effective_n_mc_ci95")
        note = "" if row.get("effective_n") is not None else "; ".join(row.get("reasons", []))
        t.add_row(
            row["population"],
            _effn_text(row),
            Text(_fmt_ci(mc), style="cyan") if mc else Text("-", style="dim"),
            _fmt_n(row.get("incidence")),
            _fmt_frac(row.get("antigen_fraction")),
            _fmt_frac(row.get("hla_coverage")),
            Text(note[:48], style="dim red") if note else "",
        )
    body = [t]
    for w in ts.warnings:
        body.append(Text(f"! {w}", style="dim yellow"))
    title = f"[bold]{ts.gene} {ts.variant}[/]  ·  {ts.disease}"
    return Panel(Group(*body), title=title, title_align="left", border_style="blue")


def render(scores, console: Console | None = None) -> None:
    """Print one rich panel per target."""
    console = console or Console()
    for ts in scores:
        console.print(_target_panel(ts))


def render_rank(scores, console: Console | None = None) -> None:
    """Print a single table across all (target, population), most-addressable first.

    Rows with a missing effective-N sort to the bottom (kept + shown with the reason),
    never dropped, consistent with the tool's surface-never-hide rule.
    """
    console = console or Console()
    rows = []
    for ts in scores:
        for r in ts.rows():
            rows.append((ts, r))
    # sort: computable first (by effective-N desc), missing last
    rows.sort(key=lambda tr: (tr[1].get("effective_n") is None,
                              -(tr[1].get("effective_n") or 0)))

    t = Table(title="pmhc-triage — ranked by effective addressable-N",
              title_style="bold", title_justify="left")
    for col, just in [("target", "left"), ("population", "left"),
                      ("effective-N", "right"), ("MC 95%", "right"),
                      ("incidence", "right"), ("antigen", "right"),
                      ("coverage", "right"), ("note", "left")]:
        t.add_column(col, justify=just)
    for ts, r in rows:
        mc = r.get("effective_n_mc_ci95")
        note = "" if r.get("effective_n") is not None else "; ".join(r.get("reasons", []))[:40]
        t.add_row(
            f"{ts.gene} {ts.variant}",
            r["population"],
            _effn_text(r),
            Text(_fmt_ci(mc), style="cyan") if mc else Text("-", style="dim"),
            _fmt_n(r.get("incidence")),
            _fmt_frac(r.get("antigen_fraction")),
            _fmt_frac(r.get("hla_coverage")),
            Text(note, style="dim red") if note else "",
        )
    console.print(t)
    console.print(Text("effective-N = incidence x antigen-fraction x HLA-coverage; "
                       "missing rows (a factor unavailable) sort last, never dropped.",
                       style="dim"))


def render_crossreact(peptide, allele, contact, safety, agents, incidents, scanned,
                      *, top=20, console=None):
    console = console or Console()
    L = len(peptide)
    anchors = [i for i in range(L) if i not in contact]
    motif = Text()
    for i, c in enumerate(peptide):
        motif.append(c, style="dim" if i in anchors else "bold cyan")
    body = [
        Text.assemble(("peptide ", "bold"), (peptide, "bold white"), ("  on  ", "dim"),
                      (allele, "bold white"), (f"   ({L}-mer)", "dim")),
        Text.assemble(("recognition motif  ", "bold"), motif,
                      (f"   anchors P{','.join(str(i + 1) for i in anchors)} (buried); "
                       "TCR-contact in bold", "dim")),
    ]
    label, note = safety
    body.append(Text.assemble(("antigen  ", "bold"), (label, "yellow"), (f"  {note}", "dim")))

    if scanned is None:
        body.append(Text("self-scan  run with --proteome <human.fasta> to flag near-self mimics",
                         style="dim"))
    elif scanned.is_missing:
        body.append(Text(f"self-scan unavailable: {'; '.join(scanned.warnings)}", style="red"))
    else:
        hits = scanned.value
        body.append(Text.assemble(
            ("self-scan  ", "bold"),
            (f"{len(hits)} self peptide(s) within threshold", "red" if hits else "green"),
            (f"   [{scanned.provenance.method}]", "dim")))
        if hits:
            t = Table(box=None, pad_edge=False)
            for c in ("protein", "pos", "peptide", "mism"):
                t.add_column(c)
            for h in hits[:top]:
                t.add_row(h.accession, str(h.position), h.peptide, str(h.mismatches))
            body.append(t)
            if len(hits) > top:
                body.append(Text(f"  ... +{len(hits) - top} more", style="dim"))
        body.append(Text("  coarse near-self screen; blind to distant mimicry "
                         "(MAGE-A3/titin differ at 4 contacts) -- no hit is not safety",
                         style="dim italic"))

    for title, items in (("clinical agents", agents),
                         ("cross-reactivity incidents (landmines)", incidents)):
        if items:
            t = Table(box=None, pad_edge=False, title=title, title_style="bold", title_justify="left")
            for c in ("name", "antigen", "peptide", "allele", "modality/status", "note"):
                t.add_column(c)
            for a in items:
                t.add_row(a.name, a.antigen, a.peptide, a.allele,
                          f"{a.modality} / {a.status}", Text(a.note, style="dim"))
            body.append(t)

    console.print(Panel(Group(*body), title="[bold]Cross-reactivity / self-safety[/]",
                        title_align="left", border_style="red"))
