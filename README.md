# pmhc-triage

**HLA-coverage-adjusted estimates of how many patients a pMHC / T-cell immunotherapy target can actually reach — with a primary source on every number.**

## Why

Target-prioritization tools (Open Targets, etc.) are HLA-blind. But a pMHC-directed therapy — a TCR-T, a TCR-mimic antibody, a bispecific — only works if the patient's HLA can present the target peptide. A binder against *peptide-X on HLA-A\*02:01* is useless in an A\*02:01-negative patient. So **raw disease incidence overstates who a pMHC target can reach**, and by an amount that varies wildly across populations (A\*02:01 is common in Europe, A\*11:01 in East Asia).

`pmhc-triage` corrects for this:

```
effective addressable N = disease incidence
                        × antigen-positive fraction
                        × HLA population coverage (of the alleles presenting the peptide)
```

Every value carries its source, URL, query date, and method — because a sourced, reproducible number is exactly what a language model will otherwise invent.

## Install

```bash
pip install -e ".[dev]"                 # core + tests
pip install -e ".[presentation]"        # optional: MHCflurry allele prediction
```

## Quickstart (out of the box)

```bash
# 1) fetch allele frequencies from AFND (license-clean, on your machine)
pmhc-triage fetch-freqs --out freqs.tsv

# 2) score a target — incidence auto-filled from the shipped cited bundle
pmhc-triage score --gene KRAS --variant G12D --disease "pancreatic cancer" \
  --study paad_tcga_pan_can_atlas_2018 --alleles "A*03:01,A*11:01" \
  --populations Europe,EastAsia,SouthAsia --freqs freqs.tsv
```

→ per-population effective-N, a Monte-Carlo 95% interval, and a provenance record for every factor. Other commands: `discover` (disease → candidate targets), `run` (batch YAML), `validate` (fail-fast preflight).

## Design principles

- **Provenance-first** — every number carries its source; a missing input is surfaced (`None` + reason), never silently zeroed.
- **Join guard** — incidence × coverage combine only within the *same* population label, never across.
- **Uncertainty visible** — factors carry sample size + CI; `effective_N` carries a seeded Monte-Carlo interval.
- **Never bundle licensed data** — AFND (CC BY-NC) is fetched on your machine, not redistributed.

## Does the number hold up?

See [`validation/`](validation/) for the evidence, honestly reported: the tool recovers known KRAS G12D restriction (A\*03:01/A\*11:01) from sequence; HLA-adjustment reorders targets across restriction families and populations but is a flat multiplier within a shared HLA superfamily; the multiplicative model's independence bias is bounded from Marty et al. 2017. It is **not** validated against a ground-truth addressable population — none exists — so treat outputs as *sourced, reproducible estimates*, not measured truth.

## Data & license

Code: **Apache-2.0**. Data attribution and terms in [`NOTICE`](NOTICE) (GCO/GLOBOCAN incidence, cBioPortal ODbL, AFND CC BY-NC, UniProt, Open Targets, MHCflurry).
