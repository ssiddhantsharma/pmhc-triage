"""Command-line entry point.

Subcommands are declared here so the interface is visible from day one; each is
wired up as its module lands. The scaffold intentionally errors clearly on
not-yet-built commands rather than pretending to work.
"""

from __future__ import annotations

import argparse
import sys


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pmhc-triage",
        description="HLA-coverage-adjusted effective addressable-population "
        "estimates for pMHC / T-cell immunotherapy targets.",
    )
    sub = p.add_subparsers(dest="command", metavar="{score,discover,run}")

    s = sub.add_parser("score", help="score one explicit target (gene + variant + disease)")
    s.add_argument("--gene")
    s.add_argument("--variant")
    s.add_argument("--disease")
    s.add_argument("--populations", help="comma-separated population names")
    s.add_argument("--out")

    d = sub.add_parser("discover", help="disease -> candidate pMHC targets via Open Targets")
    d.add_argument("--disease")
    d.add_argument("--top", type=int, default=20)
    d.add_argument("--out")

    r = sub.add_parser("run", help="batch score targets from a YAML config")
    r.add_argument("--config")
    r.add_argument("--report", choices=["md"], default=None)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.command:
        build_parser().print_help()
        return 0
    print(
        f"'{args.command}' is not implemented yet in this scaffold. "
        "Built so far: provenance primitive + hla population-coverage core.",
        file=sys.stderr,
    )
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
