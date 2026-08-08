"""Rich terminal rendering of scored targets.

Kept separate from the core so the package stays importable without rich; the CLI
falls back to plain text (``cli._print_summary``) if rich is unavailable. Nothing
here computes anything -- it only presents a :class:`~pmhc_triage.score.TargetScore`.

Two views:
- :func:`render`         -- one panel per target, per-population factor table.
- :func:`render_rank`    -- one table across all (target, population), sorted by
                            effective addressable-N (the "which target reaches the
                            most patients" view). Missing rows sort to the bottom,
                            never silently dropped.
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
    never dropped -- consistent with the tool's surface-never-hide rule.
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
