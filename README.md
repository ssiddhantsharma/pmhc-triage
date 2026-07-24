# pmhc-triage

**HLA-coverage-adjusted effective addressable-population estimates for pMHC / T-cell immunotherapy targets** — with a source on every number.

## Statement of need

Target-prioritization tools (e.g. Open Targets) are disease-agnostic and ignore the one variable unique to pMHC-directed therapeutics: **HLA restriction**. A binder against *peptide-X on HLA-A\*02:01* is useless in an A\*02:01-negative patient, so raw disease incidence overstates who a pMHC therapeutic can reach.

`pmhc-triage` estimates the **effective addressable patient population**:

```
effective addressable N = disease incidence
                        × antigen-positive fraction
                        × HLA population coverage (of the alleles presenting the peptide)
```

…attaches Open Targets tractability as context, and carries a **provenance record on every value** (source, URL, query date, method) — because sourced, reproducible numbers are exactly what a language model can't be trusted to give.

## Install

```bash
pip install -e ".[dev]"            # core + test tooling
pip install -e ".[presentation]"   # optional: MHCflurry allele prediction (torch backend)
mhcflurry-downloads fetch models_class1_presentation   # once, if using presentation
```

Core has no heavy deps (numpy, pandas, httpx, pyyaml). MHCflurry is optional — the manual-allele path runs without it.

## Commands

```bash
# 1) validate — fail-fast preflight (study resolves? populations in freq file? variant WT matches UniProt?)
pmhc-triage validate --gene KRAS --variant G12D --disease PDAC \
  --study paad_tcga_pan_can_atlas_2018 --alleles "A*11:01,A*03:01" \
  --populations Europe,EastAsia --freqs afnd.tsv --burden burden.csv --uniprot P01116

# 2) score — one target from flags -> results.csv + results.provenance.json
pmhc-triage score --gene KRAS --variant G12D --disease PDAC \
  --study paad_tcga_pan_can_atlas_2018 --alleles "A*11:01,A*03:01" \
  --populations Europe,EastAsia --freqs afnd.tsv --burden burden.csv \
  --cache .cache --out kras.csv

# 3) discover — disease -> candidate targets (Open Targets association; feed into `score`)
pmhc-triage discover --disease "pancreatic cancer" --top 20 --out targets.csv

# 4) run — batch from a YAML config
pmhc-triage run --config targets.yaml --out batch.csv
```

`--cache DIR` makes any run reproducible/offline: responses are stored with the real fetch datetime and replayed byte-identically (provenance reports the original fetch time, not the re-run).

### Advanced modes

```bash
# predict presenting alleles instead of supplying them (needs [presentation] extra)
pmhc-triage score --gene KRAS --variant G12D --disease PDAC --study paad_tcga_pan_can_atlas_2018 \
  --predict-alleles --uniprot P01116 --presentation-threshold 2 --populations Europe --freqs afnd.tsv --burden b.csv

# pool the antigen fraction across studies and/or variants (per-study N recorded)
pmhc-triage score ... --studies "paad_tcga_pan_can_atlas_2018,paad_cptac_2021" --variants "G12D,G12V,G12C"

# expression antigen (cancer-testis antigens like PRAME): RNA-seq z-score instead of a mutation
pmhc-triage score --gene PRAME --antigen-mode expression --expression-threshold 1.0 --disease ... --study ...

# HLA class II works too (DRB1 is a valid proxy; DQ/DP carry an alpha-beta heterodimer caveat)
```

### YAML config (`run`)

```yaml
targets:
  - gene: KRAS
    variant: G12D
    disease: PDAC
    study: paad_tcga_pan_can_atlas_2018
    alleles: [A*11:01, A*03:01]
    populations: [Europe, EastAsia]
    freqs: afnd.tsv
    burden: {Europe: 100000, EastAsia: 150000}   # inline, or a path to a CSV
    uniprot: P01116                                # optional, enables WT validation
```

## Inputs you supply (not bundled)

- **`--freqs`** — an AFND-format allele-frequency table you export (AFND is CC BY-NC; we never redistribute it).
- **`--burden`** — an incidence CSV (`disease,population,incidence`), or use the shipped cited GLOBOCAN-2022 starter bundle via `pmhc_triage.load_bundled(...)` (World-level; refine per-region).
- **`--alleles`** — presenting alleles (manual), or predict them with the optional MHCflurry path.

## Modules

| Module | Role | Source |
|---|---|---|
| `provenance` | every value is a `Sourced` datum (value/source/url/date/method + `extra`) | — |
| `hla` | diploid population coverage (load-bearing core) | AFND freqs, Bui et al. 2006 method |
| `sequences` | canonical protein sequence | UniProt REST (CC BY 4.0) |
| `peptides` | variant → spanning k-mers (refuses on WT mismatch) | — |
| `presentation` | peptides → presenting alleles | manual, or MHCflurry (Apache-2.0) |
| `antigen` | antigen-positive fraction (+ N, Wilson CI) | cBioPortal REST (ODbL) |
| `burden` | incidence (manual / CSV / cited bundle) | GLOBOCAN 2022 (bundle) |
| `opentargets` | tractability (context) + associated targets (discovery) | Open Targets (CC0) |
| `identity` | conservative population/disease aliasing + join hints | — |
| `caching` | on-disk, query-date-stamped, offline replay | — |
| `score` | the exit gate: join-guard, missing→None (never 0), provenance log | — |
| `pipeline` | orchestrator (`run_target`) + `preflight` | — |

## Design invariants

- **Provenance-first.** Every emitted number carries its source; a missing input is surfaced (`None` + reason), never silently zeroed.
- **Join guard.** Incidence × coverage are combined only per identical population label — never cross-population.
- **Uncertainty visible.** Fractions carry sample size N and a Wilson 95% CI.
- **MHC-only, on purpose.** The HLA multiplier only exists for pMHC/T-cell targets; other modalities are out of scope (that's Open Targets' / other tools' turf).

## Limitations — read before trusting a number

These are real and mostly **not** fixable in code; the tool surfaces them rather than hiding them.

1. **The estimate is unvalidated against ground truth.** The test suite verifies *plumbing* (that each factor is fetched and combined correctly), not that `effective_N` matches any real-world addressable population. Treat outputs as *sourced, reproducible estimates*, not measured truth.
2. **The multiplicative model assumes independence.** `incidence × antigen-fraction × HLA-coverage` treats the three as independent. They aren't: HLA-driven immunoediting can select against presented neoantigens, so antigen-positivity and HLA type correlate. The point estimate is biased by an unknown amount/direction. (The score output carries this warning.)
3. **Uncertainty is under-reported.** `effective_n_ci95_antigen_only` propagates *only* the antigen-fraction sampling CI. Incidence and HLA-coverage uncertainty are **not** modeled, so the true interval is wider than shown. And cBioPortal cohorts are convenience samples (referral/sequencing bias), so even the antigen Wilson CI is optimistically narrow.
4. **Curated burden is World-level.** The shipped bundle can't legitimately join with a population-specific HLA coverage (only with a `World` population). Supply per-region incidence for population-specific estimates; three of the five bundle figures are approximate (flagged at runtime).
5. **MHCflurry threshold is a lever, not a fact.** Which alleles count as "presenting" depends on `--presentation-threshold` (default 2%); a looser threshold admits more alleles and inflates coverage. Sensitivity is on you.

## Data-source licenses (verified)

UniProt **CC BY 4.0** · Open Targets data **CC0** / code **Apache-2.0** · MHCflurry **Apache-2.0** · cBioPortal **ODbL** (some studies restrict commercial use) · AFND **CC BY-NC** (not bundled — you supply it). This package is **Apache-2.0** and ships no non-commercially-licensed data.

## Tests

```bash
pytest -q        # hermetic (mocked APIs); a gated live MHCflurry test runs if installed
```

## License

Apache-2.0.
