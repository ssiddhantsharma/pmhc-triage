# pmhc-triage — validation log

Goal (from the plan): move the tool from *"plumbing verified"* to *"the number is
believable"* via two findings — **(1)** does HLA-adjustment reorder pMHC target
priorities vs raw incidence, and **(2)** how large is the bias from the tool's
independence assumption (immunoediting). Every number below is fetched from a
**primary source**; anything not cleanly fetchable from this environment is marked
BLOCKED with the exact one-step human unblock, **not** faked.

Reproduced with the tool's own functions on 2026-07-24/25.

> Note: `kras_freqs_nmdp.tsv` (raw AFND frequencies, CC BY-NC) is **not committed** —
> regenerate with `pmhc-triage fetch-freqs --alleles "A*02:01,A*03:01,A*11:01,C*03:04,C*04:01,C*06:02,C*07:01,C*07:02"`.

---

## Status at a glance

| Leg | Source | Fetchable here? | State |
|---|---|---|---|
| Antigen-positive fraction | cBioPortal REST | ✅ live | **VERIFIED** (`antigen_fractions_cbioportal.json`) |
| Neoantigen HLA restriction | MHCflurry 2.x (local) | ✅ local | **VERIFIED** (`kras_variant_restriction_mhcflurry.json`) |
| HLA allele frequency + N | AFND | ✅ parseable (HTML) | column map solved; reference-population choice pending |
| Incidence — World | in-repo cited bundle | ✅ | **VERIFIED** (GLOBOCAN 2022, Bray 2024) |
| Incidence — per region | GCO / GLOBOCAN | ❌ SPA/redirect | **BLOCKED** → 1 manual export |
| Per-patient TCGA HLA type | GDC PanImmune / TCIA | ❌ controlled / SPA | **BLOCKED** → dbGaP-auth download |

---

## VERIFIED result #1 — the tool recovers known KRAS G12D biology

`antigen_fractions_cbioportal.json` (live cBioPortal, TCGA PanCanAtlas):

| Target | Cancer | fraction | n | Wilson 95% |
|---|---|---|---|---|
| KRAS **G12D** | PDAC | **0.274** | 49/179 | 0.214–0.343 |
| KRAS G12V | PDAC | 0.184 | 33/179 | 0.134–0.248 |
| KRAS G12R | PDAC | 0.140 | 25/179 | 0.096–0.198 |
| KRAS G12D | CRC | 0.109 | 58/534 | 0.085–0.138 |
| KRAS G12C | LUAD | 0.124 | 70/566 | 0.099–0.153 |
| TP53 R175H | BRCA | 0.020 | 21/1066 | 0.013–0.030 |

G12D-in-PDAC 27.4% matches the canonical literature figure (~27%). ✅

**Restriction (MHCflurry, strict 0.5% presentation percentile):**
KRAS G12D → **{A\*03:01, A\*11:01, C\*03:04, C\*04:01}**. A\*03:01 / A\*11:01 is the
textbook A3-superfamily restriction of KRAS G12D (Wang & Rosenberg). The tool
recovers it from sequence alone. ✅ *(C\*03:04/C\*04:01 are additional predicted
restrictions, less characterized — reported, not asserted.)*

## VERIFIED result #2 — threshold sensitivity is large and real (dogfoods the new sweep)

Same G12D peptides, my new `sweep_presenting_alleles` over real MHCflurry output:

| threshold | # presenting alleles |
|---|---|
| ≤ 0.5% | 4 (the defensible restriction) |
| ≤ 1.0% | 9 |
| ≤ 2.0% | 13 (the old default — **3×** the strict set) |
| ≤ 5.0% | 17 (nearly the whole panel) |

The presentation threshold is *the* dominant lever on coverage, exactly as the
`sensitivity` module warns. **Any coverage/effective-N number must state its threshold;
0.5% is the defensible default for a specific neoantigen, not 2%.**

---

## VERIFIED result #3 — ranking-flip within PDAC: NULL (and that's the finding)

Completed end-to-end on fully open data (cBioPortal antigen + MHCflurry restriction +
AFND frequencies). Reference cohorts (large, real, spanning the A\*11:01 gradient):
**Europe = Germany DKMS (N=3.46M)**, **East Asia = Hong Kong Chinese BMDR (N=7,595)**,
**South Asia = India South UCBB (N=11,446)**. Frequencies+N in `kras_freqs_nmdp.tsv`;
result in `ranking_flip_result.json`; figure `ranking_flip_figure.png`.

**HLA-blind order G12D > G12V > G12R > G12C is preserved in all three populations — no flip.**
Reason (real biology): every KRAS codon-12 variant shares A3-superfamily restriction
(A\*03:01 + A\*11:01), and those alleles are **anti-correlated** across populations
(A\*03:01 ≈0.15 EUR / 0.01 EAS; A\*11:01 ≈0.05 EUR / 0.30 EAS), so each variant stays
robustly addressable everywhere (coverage 0.53–0.83). KRAS neoantigens are effectively
**pan-population** via the A3 superfamily — a genuinely useful observation.

**But HLA-adjustment is not inert.** In Europe it **compresses the G12D-vs-G12V gap from
49% → 6%** (0.274 vs 0.184 antigen → 0.162 vs 0.152 addressable), because G12V additionally
gains A\*02:01 (28% in EUR). G12D and G12V become a statistical tie (antigen CIs overlap).

**Refined thesis (honest):** HLA-adjustment reorders targets **across** HLA-restriction
families and drives **cross-population absolute-N** differences — **not** within a shared-
restriction family, where it acts as a near-common multiplier that can still compress/expand
apparent gaps. A true rank-flip demo needs differently-restricted targets (→ cross-indication,
which needs the blocked per-region incidence) — see Finding #1.

## VERIFIED result #5 — capstone: end-to-end on 100% primary-sourced data

The full pipeline run for KRAS G12D in pancreas with **every input traced to a primary
source**: incidence = GCO Cancer Today 2024 (per-country UI exports, `burden_per_region.csv`),
antigen = cBioPortal live (0.274), coverage = AFND (`kras_freqs_nmdp.tsv`). No placeholders.
Result in `capstone_per_region.json`, figure `capstone_figure.png`.

**Effective addressable-N / year (KRAS G12D pancreas, TCR-mimic-addressable):**

| region (cohort) | pancreas incidence | × antigen | × coverage | = effective-N | MC 95% |
|---|---|---|---|---|---|
| Europe (Germany) | 22,941 | 0.274 | 0.591 | **3,708** | 2,877–4,646 |
| EastAsia (China) | 122,597 | 0.274 | 0.662 | **22,229** | 17,248–27,835 |
| SouthAsia (India) | 20,477 | 0.274 | 0.532 | **2,983** | 2,316–3,742 |

China dominates (6× Europe) — and honestly, for KRAS this is **incidence-driven**, not HLA:
China's pancreas incidence is 5.3× Germany's, while KRAS's A3-superfamily coverage is nearly
flat (0.53–0.66) across all three. *(Antigen fraction is applied population-agnostically from
a mostly-European TCGA cohort — a documented caveat.)*

**The HLA lever, isolated (same pancreas incidence × 0.27, restriction varied):**

| target restriction | Europe | EastAsia | SouthAsia | best region |
|---|---|---|---|---|
| A3-superfamily (real KRAS) | 2,267 | 17,209 | 2,013 | EastAsia |
| A\*02:01-restricted | 3,018 | 11,335 | 504 | EastAsia |
| A\*11:01-restricted | 637 | 16,733 | 1,490 | EastAsia |

This is the decision-relevant HLA signal, now on real incidence: holding the German pancreas
incidence fixed, the addressable-N swings **4.7× (637 → 3,018) purely by which allele restricts
the target**, and the *best-suited restriction* reorders by region — **A\*02:01 targets favour
Europe, A\*11:01/A3 targets favour East/South Asia**. That is exactly the cross-restriction /
skewed-population value the tool exists to quantify, here on 100% sourced data.

## VERIFIED result #4 — independence-assumption bias, bounded from Marty 2017

The tool multiplies `antigen × coverage` as if independent. Marty et al. 2017 (*Cell*
171:1272) show they are not (HLA-driven immunoediting). Their supplement is **openly
downloadable** from the Elsevier CDN (`ars.els-cdn.com/.../S0092867417311443-mmc4/5.xlsx`);
I read the numbers directly (not from memory / not via summary). Extracted to
`marty2017_independence_bias.json`.

**Effect is real but modest and target-specific:**
- Recurrent drivers concentrate in poorly-presented residues: worst-presented quartile
  averages **2.34×** the TCGA recurrence of the best-presented quartile (mean count 13.5
  vs 5.8, n=1000). *(Global Spearman ≈ 0 — the signal lives in the recurrent tail, not
  across all singletons.)*
- **Patient-level** signal is broad: within-patient enrichment for poorly-presented
  mutations is FDR-significant in **18/30 tissues**, incl. **PAAD OR=1.42 (p=3.8e-9)**.
- **Per-mutation** signal is weak: **0/30 tissues** FDR-significant on the within-mutation
  test (median OR 1.06). So for any *single* target the depletion is usually not individually
  significant.

**For the tool's flagship target it's ~a non-issue:** KRAS **G12D** is the *least*-presented
common KRAS variant — **0%** of the population presents it at PHBR<1, **28.7%** at PHBR<4.
Immunoediting had almost no presentable substrate to act on, so `antigen ⊥ coverage` is
approximately safe for G12D. (Ranking of KRAS G12 by presentability, PHBR<4:
G12D 0.29 < G12S 0.50 < G12C 0.54 < G12R 0.60 < G12A 0.69 < G12V 0.82 — anti-correlated
with recurrence, the immunoediting signature.)

**Bound on the tool's bias:** where editing acted, the independent product *overestimates*
the present-competent addressable fraction; magnitude is bounded by the ORs (≲1.4× at the
population level in PAAD, negligible for poorly-presented targets like G12D).

**Coverage calibration vs Marty PHBR (corrected).** I earlier claimed the tool's coverage
is "~2× optimistic" — that was from G12D alone and is **retracted**. Comparing all four
KRAS variants (`coverage_calibration.json`): tool coverage / Marty patient-level PHBR<4 =
**mean 1.30×, range 0.86×–2.06×** (G12D 2.06, G12V 1.00, G12R 1.27, G12C 0.86 — the tool
sometimes *under*states). So it is **not** a clean systematic bias: the two metrics
(≥1 allele presents any spanning peptide at MHCflurry 0.5% vs harmonic-mean-best-rank over
the diploid genotype) measure related-but-different quantities. Honest guidance: read
`coverage` as order-of-magnitude presentability, not a calibrated per-patient probability;
it agrees with PHBR to ~1.3× on average but can differ ~2× for a given variant.

**Why the real-target cross-restriction flip is biologically hard (a finding, not a gap).**
Attempting KRAS G12C (A3) vs TP53 R175H (A\*02:01) in LUAD: R175H has **no** predicted
presenter at strict 0.5% and is rare in LUAD (3/566). That is the Marty immunoediting effect
directly — R175H persists *because* it is poorly presented — so common, well-presented,
differently-restricted neoantigen *pairs* in one cancer are rare by construction. The
cross-restriction value is therefore demonstrated with the real-frequency A\*02-vs-A\*11
coverage reordering above (result #3 / the CLI dogfood), not with a forced real-pair flip.

## Finding #1 — ranking-flip (does HLA-adjustment reorder targets?)

**Cleanest incidence-free design (fully open data):** within *one* indication (PDAC),
rank the KRAS codon-12 variants by antigen-fraction alone vs antigen-fraction × HLA
coverage. Same disease ⇒ incidence cancels, so only cBioPortal (✅) + AFND (✅) +
MHCflurry restriction (✅) are needed. A flip = HLA-adjustment materially reorders
targets = the finding.

**Remaining step (methodological, yours to set):** which AFND cohorts define each
population label. AFND aggregates many small studies per allele; a defensible choice is
the large NMDP reference cohorts (European / East-Asian / South-Asian, N ~10⁴–10⁶). I
will not pick these unsupervised — wrong reference cohorts silently corrupt coverage
(the #7 trap). Column map is solved: AFND row = `[1]=allele, [3]=population,
[5]=allele_freq, [7]=sample_size`; sample_size feeds the new Monte-Carlo coverage CI.

**Cross-indication / absolute-N ranking** additionally needs per-region incidence →
BLOCKED (see below).

## Finding #2 — independence-assumption bias (immunoediting)

The tool multiplies `incidence × antigen × coverage` as if independent. They are not:
HLA-driven immunoediting selects against presentable neoantigens (Marty et al. 2017,
*Cell* 171:1272 — established effect). In-cohort test: among carriers of a G12D-presenting
allele (A\*03:01/A\*11:01), is G12D frequency lower than in non-carriers?

- Mutations: cBioPortal / mc3 PUBLIC MAF — ✅ open.
- **Per-patient HLA genotype: BLOCKED.** The TCGA HLA/neoantigen files
  (`TCGA_pMHC_SNV_sampleSummary_MC3_v0.2.8.CONTROLLED_170404.tsv`,
  `..._sample_neoantigens_10062017.tsv`) on the GDC PanImmune page are **dbGaP
  controlled-access**. TCIA (tcia.at) is a SPA.
- **Fallback with zero blocked data:** cite Marty 2017's published magnitude as the
  literature bound on the bias (pull the exact numbers from the paper, not memory).

---

## Blockers → exact one-step unblocks

1. **Per-region incidence.** GCO (`gco.iarc.who.int`) and all factsheet PDFs return a
   3095-byte SPA shell; every `gco-api`/`today/api` endpoint 301-redirects to it.
   → **You** open GCO Cancer Today, filter to the cancer/region, click *Download CSV*
   (a human click through JS). Drop the CSV in `validation/`; the tool's
   `load_burden_table` parses it directly. This is the "world-burden → per-region" task.
2. **Per-patient TCGA HLA.** → **You** download the CONTROLLED file above with your
   dbGaP/TCGA credentials (UUID list captured in `data_sources.md`). Then the
   immunoediting test runs on fully real data.

## Reproduce

```bash
cd ~/pmhc-triage && source .venv/bin/activate
# antigen leg (live):        python validation/scripts... (see git log / commands in this session)
# restriction + sweep:       needs the [presentation] extra (mhcflurry)
```

Data captured: `antigen_fractions_cbioportal.json`, `kras_variant_restriction_mhcflurry.json`.
