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

from .montecarlo import prob_greater

# P(top > row) at or above this counts as separable. 0.95 is the usual convention, stated
# here rather than buried: it is a threshold, not a measurement.
SEPARABLE_P = 0.95


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
    # carry the PopulationScore, not just its dict: the Monte-Carlo draws hang off it and a
    # sorted column of point estimates reads as a ranking whether or not it is one.
    rows = [(ts, ps.to_dict(), ps.mc) for ts in scores for ps in ts.per_population.values()]
    # sort: computable first (by effective-N desc), missing last
    rows.sort(key=lambda tr: (tr[1].get("effective_n") is None,
                              -(tr[1].get("effective_n") or 0)))
    top_mc = next((mc for _, _, mc in rows if mc is not None and mc.draws is not None), None)

    t = Table(title="pmhc-triage — ranked by effective addressable-N",
              title_style="bold", title_justify="left")
    for col, just in [("target", "left"), ("population", "left"),
                      ("effective-N", "right"), ("MC 95%", "right"),
                      ("P(top>)", "right"),
                      ("incidence", "right"), ("antigen", "right"),
                      ("coverage", "right"), ("note", "left")]:
        t.add_column(col, justify=just)
    n_tied = 0
    for ts, r, ps_mc in rows:
        mc = r.get("effective_n_mc_ci95")
        p_top = None if ps_mc is None or top_mc is None else prob_greater(top_mc, ps_mc)
        if p_top is None:
            ptxt = Text("-", style="dim")
        elif ps_mc is top_mc:
            ptxt = Text("top", style="bold green")
        elif p_top >= SEPARABLE_P:
            ptxt = Text(f"{p_top:.2f}", style="green")
        else:
            ptxt = Text(f"{p_top:.2f}", style="yellow")
            n_tied += 1
        note = "" if r.get("effective_n") is not None else "; ".join(r.get("reasons", []))[:40]
        t.add_row(
            f"{ts.gene} {ts.variant}",
            r["population"],
            _effn_text(r),
            Text(_fmt_ci(mc), style="cyan") if mc else Text("-", style="dim"),
            ptxt,
            _fmt_n(r.get("incidence")),
            _fmt_frac(r.get("antigen_fraction")),
            _fmt_frac(r.get("hla_coverage")),
            Text(note, style="dim red") if note else "",
        )
    console.print(t)
    console.print(Text("effective-N = incidence x antigen-fraction x HLA-coverage; "
                       "missing rows (a factor unavailable) sort last, never dropped.",
                       style="dim"))
    if top_mc is not None:
        msg = (f"P(top>) is P(top row > this row) over INDEPENDENT Monte-Carlo draws, not interval "
               f"overlap. {n_tied} row(s) are NOT separable from the top at P>={SEPARABLE_P:.2f}"
               + (": read the order, not the ranking." if n_tied else "."))
        console.print(Text(msg, style="yellow" if n_tied else "dim"))


def _agent_table(items, *, header):
    t = Table(box=None, pad_edge=False, title_style="bold", header_style=f"bold {header}")
    for c in ("name", "antigen", "peptide", "allele", "modality/status", "note"):
        t.add_column(c)
    for a in items:
        t.add_row(Text(a.name, style="bold"), Text(a.antigen, style="cyan"),
                  Text(a.peptide, style="bold magenta"), Text(a.allele, style="green"),
                  f"{a.modality} / {a.status}", Text(a.note, style="dim"))
    return t


def _benign_panel(benign, top):
    """Mimics eluted from benign tissue. Stronger evidence than a proteome match: these were
    measured, not merely translated."""
    if benign.value is None:
        return Panel(Text("; ".join(benign.warnings) or "unavailable", style="yellow"),
                     title="[bold]Benign-tissue ligands[/]", title_align="left",
                     border_style="yellow")
    hits = benign.value
    body = [Text.assemble(
        (f"{len(hits)} benign ligand(s) share the contact face  ",
         "bold red" if hits else "bold green"),
        (f"[{benign.provenance.method}]", "dim"))]
    if hits:
        t = Table(box=None, pad_edge=False, header_style="bold red")
        for c in ("peptide", "mism", "allele call", "tissues"):
            t.add_column(c)
        for h in hits[:top]:
            mc = "bold red" if h.mismatches <= 1 else ("yellow" if h.mismatches == 2 else "dim")
            call = {"s": Text("strong", style="bold red"), "w": Text("weak", style="yellow"),
                    "n": Text("non", style="dim")}.get(h.binding_call, Text("-", style="dim"))
            t.add_row(Text(h.peptide, style="bold magenta"), Text(str(h.mismatches), style=mc),
                      call, Text(", ".join(h.tissues[:4]) or "-", style="cyan"))
        body.append(t)
        if len(hits) > top:
            body.append(Text(f"... +{len(hits) - top} more", style="dim"))
    body.append(Text("the peptide was MEASURED on benign tissue; the allele call is PREDICTED "
                     "from the donor genotype, and tissue cannot be paired to allele here",
                     style="dim italic"))
    return Panel(Group(*body), title="[bold]Benign-tissue ligands (HLA Ligand Atlas)[/]",
                 title_align="left", border_style="red" if hits else "green")


def render_crossreact(peptide, allele, contact, safety, agents, incidents, scanned,
                      *, top=20, console=None, benign=None):
    console = console or Console()
    L = len(peptide)
    anchors = [i for i in range(L) if i not in contact]
    motif = Text()
    for i, c in enumerate(peptide):
        motif.append(c, style="dim yellow" if i in anchors else "bold bright_cyan")
    label, note = safety
    safe_color = "green" if label == "neoantigen" else "yellow"
    head = Group(
        Text.assemble(("peptide  ", "bold"), (peptide, "bold bright_white"), ("  on  ", "dim"),
                      (allele, "bold bright_white"), (f"   ({L}-mer)", "dim")),
        Text.assemble(("motif    ", "bold"), motif,
                      (f"   anchors P{','.join(str(i + 1) for i in anchors)} ", "dim yellow"),
                      ("(buried); ", "dim"), ("TCR-contact", "bold bright_cyan"),
                      (" is what a binder reads", "dim")),
        Text.assemble(("antigen  ", "bold"), (label, f"bold {safe_color}"), (f"  {note}", "dim")),
    )
    console.print(Panel(head, title="[bold]Recognition motif / safety[/]",
                        title_align="left", border_style="cyan"))

    if scanned is None:
        console.print(Panel(
            Text("run with --proteome <human.fasta> to flag near-self mimics", style="dim"),
            title="[bold]Proteome self-scan[/]", title_align="left", border_style="blue"))
    elif scanned.is_missing:
        console.print(Panel(Text('; '.join(scanned.warnings), style="red"),
                            title="[bold]Proteome self-scan[/]", title_align="left",
                            border_style="red"))
    else:
        hits = scanned.value
        body = [Text.assemble(
            (f"{len(hits)} self peptide(s) within threshold  ", "bold red" if hits else "bold green"),
            (f"[{scanned.provenance.method}]", "dim"))]
        if hits:
            annotated = any(h.presented is not None for h in hits)
            t = Table(box=None, pad_edge=False, header_style="bold red")
            for c in ("protein", "pos", "peptide", "mism") + (("presented",) if annotated else ()):
                t.add_column(c)
            for h in hits[:top]:
                mc = "bold red" if h.mismatches <= 1 else ("yellow" if h.mismatches == 2 else "dim")
                cells = [Text(h.accession, style="cyan"), str(h.position),
                         Text(h.peptide, style="bold magenta"), Text(str(h.mismatches), style=mc)]
                if annotated:
                    cells.append({True: Text("yes", style="bold red"),
                                  False: Text("no", style="green"),
                                  None: Text("?", style="dim")}[h.presented])
                t.add_row(*cells)
            body.append(t)
            if annotated:
                body.append(Text("presented = MHCflurry predicts this allele displays the self "
                                 "peptide; 'no' is a likely false positive, not a cleared one "
                                 "(hits are annotated, never dropped)", style="dim"))
            if len(hits) > top:
                body.append(Text(f"... +{len(hits) - top} more", style="dim"))
        body.append(Text("coarse near-self screen; blind to distant mimicry "
                         "(MAGE-A3/titin differ at 4 contacts) -- no hit is not safety",
                         style="dim italic"))
        console.print(Panel(Group(*body), title="[bold]Proteome self-scan[/]",
                            title_align="left", border_style="red" if hits else "green"))

    if benign is not None:
        console.print(_benign_panel(benign, top))

    if agents:
        console.print(Panel(_agent_table(agents, header="green"),
                            title="[bold green]Clinical agents[/]", title_align="left",
                            border_style="green"))
    if incidents:
        console.print(Panel(_agent_table(incidents, header="red"),
                            title="[bold red]Cross-reactivity incidents (landmines)[/]",
                            title_align="left", border_style="red"))
