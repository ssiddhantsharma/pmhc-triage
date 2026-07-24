# pmhc-triage

**HLA-coverage-adjusted effective addressable-population estimates for pMHC / T-cell immunotherapy targets.**

## Statement of need

Target prioritization tools (e.g. Open Targets) are disease-agnostic and ignore the one variable unique to pMHC-directed therapeutics: **HLA restriction**. A binder against *peptide-X on HLA-A\*02:01* is useless in an A\*02:01-negative patient. Raw disease incidence therefore overstates who a pMHC therapeutic can actually reach.

`pmhc-triage` estimates the **effective addressable patient population** for a pMHC target:

```
effective addressable N = disease incidence
                        × antigen-positive fraction
                        × HLA population coverage (of the alleles that present the peptide)
```

…and attaches Open Targets tractability as context. **Every output number carries its source, URL, query date, and derivation method** (provenance-first): the point is sourced, reproducible numbers — exactly what a language model cannot be trusted to give.

## Scope

| In scope | Out of scope |
|---|---|
| pMHC / T-cell targets (HLA-restricted) | Non-pMHC modalities — the formula (`× HLA coverage × presentation`) is undefined without an MHC |
| Disease → candidate-target discovery (chains Open Targets, optional) | Computed TAM / pricing forecast — that's judgment; faking it discredits the provenance premise |
| Analog pricing *anchors* (sourced reference, never multiplied) | Disease-agnostic discovery as a product (that's Open Targets' job) |

**Principle:** compute anything with a full sourced provenance chain; refuse anything that requires a forecast assumption.

## Modules

| Module | Role | Source |
|---|---|---|
| `provenance` | every value is a `Sourced` datum with its origin | — |
| `hla` | diploid population coverage (**load-bearing core**) | AFND, Bui et al. 2006 method |
| `burden` | disease incidence | curated GLOBOCAN/GBD bundle *(planned)* |
| `antigen` | antigen-positive fraction | cBioPortal REST *(planned)* |
| `presentation` | peptide → presenting alleles | MHCflurry *(planned)* |
| `tractability` | druggability context (not a driver of the score) | Open Targets GraphQL *(planned)* |
| `score` | combine + rank + provenance log | — *(planned)* |

## Status

Early scaffold. Built: the **provenance primitive** and the **`hla` diploid-coverage core** (validated against analytically-derived values; to be cross-checked against the IEDB Population Coverage tool). API-calling modules are next.

## Install

```bash
pip install -e ".[dev]"          # core + test tooling
pip install -e ".[presentation]" # adds MHCflurry (pulls TensorFlow)
```

## License

Apache-2.0.
