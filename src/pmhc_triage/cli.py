"""Command-line entry point: score (flags), run (YAML config), validate (preflight).

score/run compute effective addressable-N and write results.csv + provenance.json.
validate is the fail-fast preflight -- it checks joins and exits non-zero on issues
*before* any run.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

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
    )


def _print_summary(spec: TargetSpec, ts) -> None:
    print(f"\n{spec.gene} {spec.variant} / {spec.disease}")
    for row in ts.rows():
        if row["effective_n"] is not None:
            print(
                f"  {row['population']:12s} effective_N = {row['effective_n']:,.0f}  "
                f"(incidence {row['incidence']:,.0f} x antigen {row['antigen_fraction']} "
                f"x coverage {row['hla_coverage']})"
            )
        else:
            print(f"  {row['population']:12s} effective_N = n/a  ({'; '.join(row['reasons'])})")
    for w in ts.warnings:
        print(f"  ! {w}")


def _add_score_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--gene", required=True)
    p.add_argument("--variant", required=True, help="e.g. G12D")
    p.add_argument("--disease", required=True)
    p.add_argument("--study", required=True, help="cBioPortal study id")
    p.add_argument("--alleles", required=True, help="comma-separated, e.g. 'A*11:01,A*03:01'")
    p.add_argument("--populations", required=True, help="comma-separated")
    p.add_argument("--freqs", help="AFND-format frequency file")
    p.add_argument("--burden", help="burden CSV (disease,population,incidence)")
    p.add_argument("--uniprot", help="UniProt accession (enables variant WT validation)")
    p.add_argument("--out", default="results.csv")
    p.add_argument("--provenance-out", dest="provenance_out", default=None)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pmhc-triage",
        description="HLA-coverage-adjusted effective addressable-population "
        "estimates for pMHC / T-cell immunotherapy targets.",
    )
    sub = p.add_subparsers(dest="command", metavar="{score,run,validate}")

    s = sub.add_parser("score", help="score one target from flags")
    _add_score_args(s)

    v = sub.add_parser("validate", help="preflight one target (fail-fast, exits non-zero on issues)")
    _add_score_args(v)

    r = sub.add_parser("run", help="score a batch of targets from a YAML config")
    r.add_argument("--config", required=True)
    r.add_argument("--out", default="results.csv")
    r.add_argument("--provenance-out", dest="provenance_out", default=None)

    return p


def _provenance_path(out: str, provided: str | None) -> Path:
    if provided:
        return Path(provided)
    return Path(out).with_suffix(".provenance.json")


def _cmd_score(a: argparse.Namespace) -> int:
    spec = _spec_from_args(a)
    ts = run_target(spec)
    _print_summary(spec, ts)
    csv_path = write_results_csv([ts], a.out)
    prov_path = write_provenance_json([ts], _provenance_path(a.out, a.provenance_out))
    print(f"\nwrote {csv_path} and {prov_path}")
    return 0


def _cmd_validate(a: argparse.Namespace) -> int:
    spec = _spec_from_args(a)
    issues = preflight(spec)
    if not issues:
        print("preflight OK -- inputs join cleanly.")
        return 0
    print("preflight found issues:", file=sys.stderr)
    for i in issues:
        print(f"  - {i}", file=sys.stderr)
    return 1


def _cmd_run(a: argparse.Namespace) -> int:
    import yaml

    config = yaml.safe_load(Path(a.config).read_text())
    specs = [TargetSpec.from_dict(t) for t in config.get("targets", [])]
    if not specs:
        print("no targets in config", file=sys.stderr)
        return 2
    scores = []
    for spec in specs:
        ts = run_target(spec)
        _print_summary(spec, ts)
        scores.append(ts)
    csv_path = write_results_csv(scores, a.out)
    prov_path = write_provenance_json(scores, _provenance_path(a.out, a.provenance_out))
    print(f"\nwrote {csv_path} and {prov_path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.command:
        build_parser().print_help()
        return 0
    return {"score": _cmd_score, "validate": _cmd_validate, "run": _cmd_run}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
