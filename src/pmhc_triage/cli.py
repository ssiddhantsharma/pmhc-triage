"""Command-line entry point: score (flags), run (YAML config), validate (preflight).

score/run compute effective addressable-N and write results.csv + provenance.json.
validate is the fail-fast preflight -- it checks joins and exits non-zero on issues
*before* any run.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from .caching import cached_client
from .opentargets import associated_targets, resolve_disease
from .pipeline import TargetSpec, preflight, run_target
from .report import write_provenance_json, write_results_csv


def _csv_list(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def _spec_from_args(a: argparse.Namespace) -> TargetSpec:
    return TargetSpec(
        gene=a.gene,
        variant=a.variant,
        disease=a.disease,
        study=a.study,
        alleles=_csv_list(a.alleles),
        populations=_csv_list(a.populations),
        freqs_path=a.freqs,
        burden_path=a.burden,
        uniprot=a.uniprot,
        predict_alleles=getattr(a, "predict_alleles", False),
        presentation_threshold=getattr(a, "presentation_threshold", 2.0),
        studies=_csv_list(a.studies) or None,
        variants=_csv_list(a.variants) or None,
        antigen_mode=getattr(a, "antigen_mode", "mutation"),
        expression_threshold=getattr(a, "expression_threshold", 1.0),
    )


def _print_summary(ts) -> None:
    print(f"\n{ts.gene} {ts.variant} / {ts.disease}")
    for row in ts.rows():
        if row["effective_n"] is not None:
            ci = row.get("effective_n_ci95_antigen_only")
            ci_str = f" [{ci[0]:,.0f}-{ci[1]:,.0f} antigen-CI]" if ci else ""
            print(
                f"  {row['population']:12s} effective_N = {row['effective_n']:,.0f}{ci_str}  "
                f"(incidence {row['incidence']:,.0f} x antigen {row['antigen_fraction']} "
                f"x coverage {row['hla_coverage']})"
            )
        else:
            print(f"  {row['population']:12s} effective_N = n/a  ({'; '.join(row['reasons'])})")
    for w in ts.warnings:
        print(f"  ! {w}")


def _render(scores: list, *, rank: bool = False) -> None:
    """Render scored targets with rich; fall back to plain text if rich is absent."""
    try:
        from rich.console import Console

        from .render import render, render_rank
    except ImportError:  # pragma: no cover - rich is a core dep, this is belt-and-suspenders
        for ts in scores:
            _print_summary(ts)
        return
    console = Console()
    (render_rank if rank else render)(scores, console)


def _add_score_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--gene", required=True)
    p.add_argument("--variant", default="", help="e.g. G12D (required for mutation mode)")
    p.add_argument("--antigen-mode", dest="antigen_mode", choices=["mutation", "expression"],
                   default="mutation", help="mutation (variant) or expression (RNA-seq, e.g. CTAs)")
    p.add_argument("--expression-threshold", dest="expression_threshold", type=float, default=1.0,
                   help="expression z-score cutoff for --antigen-mode expression")
    p.add_argument("--disease", required=True)
    p.add_argument("--study", required=True, help="cBioPortal study id")
    p.add_argument("--studies", default="", help="comma-separated study ids to POOL antigen fraction over")
    p.add_argument("--variants", default="", help="comma-separated variants to pool (e.g. 'G12D,G12V,G12C')")
    p.add_argument("--alleles", default="", help="comma-separated, e.g. 'A*11:01,A*03:01' (manual path)")
    p.add_argument("--predict-alleles", dest="predict_alleles", action="store_true",
                   help="predict presenting alleles via MHCflurry (needs --uniprot + [presentation] extra)")
    p.add_argument("--presentation-threshold", dest="presentation_threshold", type=float, default=2.0,
                   help="MHCflurry presentation percentile threshold (default 2.0)")
    p.add_argument("--populations", required=True, help="comma-separated")
    p.add_argument("--freqs", help="AFND-format frequency file")
    p.add_argument("--burden", help="burden CSV (disease,population,incidence)")
    p.add_argument("--uniprot", help="UniProt accession (enables variant WT validation)")
    p.add_argument("--out", default="results.csv")
    p.add_argument("--provenance-out", dest="provenance_out", default=None)
    p.add_argument("--cache", default=None, help="cache dir for reproducible/offline re-runs")
    p.add_argument("--rank", action="store_true",
                   help="render one comparison table sorted by effective addressable-N")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pmhc-triage",
        description="HLA-coverage-adjusted effective addressable-population "
        "estimates for pMHC / T-cell immunotherapy targets.",
    )
    sub = p.add_subparsers(dest="command", metavar="{score,discover,validate,run,fetch-freqs}")

    s = sub.add_parser("score", help="score one target from flags")
    _add_score_args(s)

    v = sub.add_parser("validate", help="preflight one target (fail-fast, exits non-zero on issues)")
    _add_score_args(v)

    dsc = sub.add_parser("discover", help="disease -> candidate targets (Open Targets association)")
    dsc.add_argument("--disease", required=True)
    dsc.add_argument("--top", type=int, default=20)
    dsc.add_argument("--out", default="discover.csv")
    dsc.add_argument("--cache", default=None, help="cache dir for reproducible/offline re-runs")

    r = sub.add_parser("run", help="score a batch of targets from a YAML config")
    r.add_argument("--config", required=True)
    r.add_argument("--out", default="results.csv")
    r.add_argument("--provenance-out", dest="provenance_out", default=None)
    r.add_argument("--cache", default=None, help="cache dir for reproducible/offline re-runs")
    r.add_argument("--rank", action="store_true",
                   help="render one comparison table sorted by effective addressable-N")

    ff = sub.add_parser(
        "fetch-freqs",
        help="fetch AFND allele frequencies to a local --freqs file (NOT bundled; "
        "your data-use with AFND under its CC BY-NC terms)",
    )
    ff.add_argument("--alleles", default="",
                    help="comma-separated alleles (default: a common class-I panel)")
    ff.add_argument("--out", default="freqs.tsv")

    return p


def _provenance_path(out: str, provided: str | None) -> Path:
    if provided:
        return Path(provided)
    return Path(out).with_suffix(".provenance.json")


def _client(a: argparse.Namespace):
    return cached_client(a.cache) if getattr(a, "cache", None) else None


def _cmd_score(a: argparse.Namespace) -> int:
    spec = _spec_from_args(a)
    client = _client(a)
    try:
        ts = run_target(spec, client=client)
    finally:
        if client is not None:
            client.close()
    _render([ts], rank=getattr(a, "rank", False))
    csv_path = write_results_csv([ts], a.out)
    prov_path = write_provenance_json([ts], _provenance_path(a.out, a.provenance_out))
    print(f"\nwrote {csv_path} and {prov_path}")
    return 0


def _cmd_validate(a: argparse.Namespace) -> int:
    spec = _spec_from_args(a)
    client = _client(a)
    try:
        issues = preflight(spec, client=client)
    finally:
        if client is not None:
            client.close()
    if not issues:
        print("preflight OK -- inputs join cleanly.")
        return 0
    print("preflight found issues:", file=sys.stderr)
    for i in issues:
        print(f"  - {i}", file=sys.stderr)
    return 1


def _cmd_discover(a: argparse.Namespace) -> int:
    client = _client(a)
    try:
        did = resolve_disease(a.disease, client=client)
        if did.is_missing:
            print(f"could not resolve disease {a.disease!r}: {did.warnings}", file=sys.stderr)
            return 1
        rows = associated_targets(did.value, top=a.top, client=client)
    finally:
        if client is not None:
            client.close()
    if rows.is_missing:
        print(f"no associated targets: {rows.warnings}", file=sys.stderr)
        return 1
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["symbol", "ensembl_id", "association_score"])
        w.writeheader()
        w.writerows(rows.value)
    print(f"{a.disease} ({did.value}): top {len(rows.value)} associated targets -> {a.out}")
    for r in rows.value[:10]:
        print(f"  {r['symbol']:10s} {r['association_score']:.3f}")
    print("\nNote: Open Targets association is disease-agnostic; feed these into `score` "
          "with a variant + study to get pMHC addressable-N.")
    return 0


def _cmd_fetch_freqs(a: argparse.Namespace) -> int:
    from .fetch import COMMON_CLASS_I, DEFAULT_POPULATIONS, fetch_frequencies

    alleles = _csv_list(a.alleles) or list(COMMON_CLASS_I)
    pops = DEFAULT_POPULATIONS
    print(f"fetching {len(alleles)} alleles x {len(pops)} cohorts from AFND "
          "(your data-use under AFND CC BY-NC; not redistributed)...")
    result = fetch_frequencies(alleles, pops)
    n = result.write_tsv(a.out)
    print(f"wrote {n} (allele,population) rows -> {a.out}  [cohorts: "
          + ", ".join(f"{p.label}={p.afnd_population}" for p in pops) + "]")
    if result.warnings:
        print(f"{len(result.warnings)} surfaced (not written as 0):")
        for w in result.warnings[:12]:
            print(f"  ! {w}")
    return 0


def _cmd_run(a: argparse.Namespace) -> int:
    import yaml

    config = yaml.safe_load(Path(a.config).read_text())
    specs = [TargetSpec.from_dict(t) for t in config.get("targets", [])]
    if not specs:
        print("no targets in config", file=sys.stderr)
        return 2
    client = _client(a)
    scores = []
    try:
        for spec in specs:
            scores.append(run_target(spec, client=client))
    finally:
        if client is not None:
            client.close()
    _render(scores, rank=getattr(a, "rank", False))
    csv_path = write_results_csv(scores, a.out)
    prov_path = write_provenance_json(scores, _provenance_path(a.out, a.provenance_out))
    print(f"\nwrote {csv_path} and {prov_path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.command:
        build_parser().print_help()
        return 0
    handlers = {
        "score": _cmd_score,
        "validate": _cmd_validate,
        "discover": _cmd_discover,
        "run": _cmd_run,
        "fetch-freqs": _cmd_fetch_freqs,
    }
    return handlers[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
