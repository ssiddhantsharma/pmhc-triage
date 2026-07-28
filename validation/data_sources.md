# Data-source reconnaissance (probed 2026-07-24/25)

Exact endpoints + results, so the next session doesn't re-probe.

## Fetchable now (primary, reproducible)

- **cBioPortal REST** — `https://www.cbioportal.org/api` — ✅ JSON, live. The tool's
  `variant_frequency` / `expression_positive_fraction` use it. PDAC sequenced N=179.
- **MHCflurry 2.x** — local (`[presentation]` extra). Restriction from sequence.
- **AFND** — `http://www.allelefrequencies.net/hla6006a.asp?hla_locus_type=Classical&hla_allele1=A*NN:NN`
  returns a ~325 KB HTML table. Data-row column order (0-indexed `<td>`):
  `[1]=allele  [3]=population  [4]=%individuals  [5]=allele_frequency  [7]=sample_size`.
  e.g. `A*11:01 | American Samoa | | 0.1600 | | 51`. Aggregates many studies/pop —
  choose reference cohorts deliberately (NMDP broad groups recommended).

## BLOCKED (behind SPA or controlled access)

- **GCO / GLOBOCAN incidence** — `gco.iarc.who.int` serves a 3095-byte SPA shell for
  `gco-api/v1/*`, `today/api/*`, and all `factsheets/.../*.pdf`. Old `gco.iarc.fr`
  hard-301s to it. No JSON reachable. → manual CSV export from the UI.
- **Per-patient TCGA HLA / neoantigens** — GDC PanImmune page
  (`gdc.cancer.gov/about-data/publications/panimmune`) is a real static file server
  (`https://api.gdc.cancer.gov/data/<uuid>`), but the HLA-bearing files are
  **CONTROLLED** (dbGaP):
  - `0d3ee0a7-0557-447b-9ada-bc7838d1effb` — TCGA_pMHC_SNV_sampleSummary_MC3_v0.2.8.CONTROLLED_170404.tsv
  - `5ecff039-a050-4ad3-bff0-c112cbd8d2ff` — TCGA_PCA.mc3...sample_neoantigens_10062017.tsv (CONTROLLED)
  - `1c8cfe5f-e52d-41ba-94da-f15ea1337efc` — mc3.v0.2.8.PUBLIC.maf.gz (**public** mutations, open)
  - `1a7d7be8-675d-4e60-a105-19d4121bdebf` — merged_sample_quality_annotations.tsv (open)

## Literature anchor (verify numbers from the paper, not memory)

- Marty et al. 2017, *Cell* 171:1272 — MHC-I genotype restricts the oncogenic
  mutational landscape (the immunoediting effect the tool's independence assumption ignores).
