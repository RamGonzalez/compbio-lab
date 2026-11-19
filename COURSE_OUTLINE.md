# Population Genomics: from the Terminal to PCA and ADMIXTURE

**Format:** 7 sessions × ~90 min, hands-on, each student on their own laptop (Windows + WSL2/Ubuntu).
**Audience:** Physical anthropology, biology and biotechnology students. No assumed programming or command-line experience.
**Hardware budget:** assume **8 GB RAM and 4 cores** per laptop. Every command in the course must run comfortably within that (use `--threads 4`, `plink --memory 3000`, `admixture -j4`).
**Final goal:** each student can run PCA and ADMIXTURE on a small dataset and explain what the results do and do not say.

## Design principles

- **Terminal first, but always in service of data.** From session 2 students touch real genetic files, not toy `ls` exercises.
- **Concept, then run, then interpret.** Every analysis session has: intuition (no maths first), a run, and an interpretation exercise.
- **One small dataset all course long.** Students get to know it, so by session 7 they can predict what PCA/ADMIXTURE will show and be surprised only by the interesting parts.
- **Everything is a script.** By session 4 students save commands in a `.sh` file; by session 7 they can rerun the whole analysis with one command.
- **Instructor prepares the data.** Students never download or subset reference data themselves.

## Software

One conda environment on WSL, created by the students in session 1 (installed with Miniforge, environment named `compbio-lab`):

- `plink` (1.9) and `plink2`
- `bcftools` (and `htslib` for bgzip/tabix)
- `admixture` (ADMIXTURE)
- `r-base` with `ggplot2`, `tidyverse` (plotting is done in **R** throughout)

Keep everything within 8 GB / 4 cores: small SNP sets, few runs in parallel, no memory-hungry tools.

## Data (to be prepared later)

Not needed yet. Plan when we get to it: a small 1000 Genomes subset, with **MXL** as the admixed group plus continental anchors (e.g. YRI, CEU, CHB), roughly 40 individuals per population, a few chromosomes only, so every step runs in seconds to minutes within the laptop budget. Provided as a VCF, a QC'd PLINK fileset and a `populations.tsv`. Optionally one or two "trap" variants for session 7 (relatives left in, unbalanced sample).

---

## Session 1: Welcome, setup and course overview

**Learning objectives**
- Everyone knows each other and what each student hopes to get out of the course.
- Everyone can open WSL and has a working conda environment with the course tools.
- Students know what the course covers, what they will be able to do at the end, and how sessions are run.

**Plan (90 min)**
- 0–20: Introductions (name, background, what they know about genetics/computing, one thing they want to learn). Brief course presentation: goals, the 7-session outline, expectations, how to ask for help.
- 20–35: Open WSL, first look at the prompt, `pwd`, `ls`, `cd`, working in the Linux home (`~/`) rather than `/mnt/c`.
- 35–60: Install Miniforge/conda. Explain what a package manager and an environment are, and why we use one.
- 60–80: Create the `compbio-lab` env and install `plink`, `plink2`, `bcftools`, `admixture`, R. Run each tool with `--version`/`--help`.
- 80–90: Run the `check_setup.sh` script that verifies every tool. Note who still has problems; fix them before session 2 (office hours or chat).

**Deliverable/homework:** a screenshot or the output of `check_setup.sh` showing everything installed. Reserve a lot of slack: install problems always overrun, and downloads on slow connections take time.

## Session 2: The terminal and your first genomic file

**Learning objectives**
- Navigate the filesystem, find Windows files from Linux (`/mnt/c/...`).
- Use `pwd ls cd mkdir cp mv rm less head tail wc cat`, tab completion, history.
- Use `grep`, `cut`, `sort`, `uniq`, pipes and redirection.
- Understand paths, file extensions, and compressed files.

**Plan**
- 0–10: Quick check that the env activates for everyone.
- 10–40: Guided navigation and file commands on a practice folder.
- 40–70: Pipes on real data: count samples per population from the metadata table with `cut | sort | uniq -c`; peek at a VCF with `zcat | head`.
- 70–85: First shell script that prints sample counts per population.
- 85–90: Cheat sheet handout.

**Deliverable/homework:** `session2.sh` that outputs the counts table.

## Session 3: Human genetics and genetic variation, concepts

**Learning objectives**
- Genome, chromosome, locus, allele, genotype, SNP; homozygous vs heterozygous.
- Allele frequency; why most variants are rare and a few are common.
- Hardy-Weinberg expectation (the idea and p², 2pq, q²).
- Human variation is mostly *within* populations; structure arises from drift, migration, and demographic history (Out of Africa, serial founder effects, admixture).
- Ancestry is a statistical description of genetic similarity, not a biological "race" category. Build in discussion time; the anthropology students will have a lot to add.

**Plan**
- 0–30: Interactive lecture with drawings (drift demo on paper or a tiny simulator).
- 30–60: Hands-on: pull genotypes of ~5 SNPs for a few individuals with `bcftools query`; compute an allele frequency by hand and compare with `plink --freq`.
- 60–80: Frequencies of well-known differentiated SNPs per population (e.g. *SLC24A5*, *LCT*, Duffy) and why they differ.
- 80–90: Students write predictions: how many clusters, where will MXL fall (revisited in session 7).

**Deliverable/homework:** one-page prediction sheet.

## Session 4: Genotype data in practice (VCF, PLINK, QC, LD pruning)

**Learning objectives**
- Read a VCF: header, fixed columns, genotype fields (`0|1`, `./.`).
- Know the PLINK `.bed/.bim/.fam` files and convert VCF → PLINK.
- Filter: biallelic SNPs, missingness, MAF, and why each matters.
- Understand linkage disequilibrium and why PCA/ADMIXTURE need pruned SNPs.
- Relatedness conceptually (and `plink2 --king-cutoff`).

**Plan**
- 0–20: Anatomy of a VCF (annotated printout).
- 20–50: `plink --vcf ... --make-bed`; `--missing`, `--freq`, `--hardy`; inspect outputs; one histogram in R.
- 50–70: LD concept, `--indep-pairwise 50 5 0.2`, compare SNP counts before/after.
- 70–90: Put QC commands into `01_qc.sh` producing the clean, pruned fileset used from here on.

**Deliverable/homework:** `01_qc.sh` runs from scratch; table of SNP counts after each step.

## Session 5: PCA, concepts and run

**Learning objectives**
- Intuition: PCA rotates the data to axes of greatest variance; each individual gets coordinates (PCs).
- Eigenvalues / % variance explained, loadings (briefly).
- Run `plink2 --pca` and plot PC1 vs PC2 coloured by population in R.
- Interpret clusters, gradients (admixed individuals between clusters), outliers.

**Plan**
- 0–25: Intuition without matrix algebra, using an analogy biologists will find natural, then the genotype matrix (individuals × SNPs) → PCs.
- 25–45: Run PCA, inspect `.eigenvec` / `.eigenval`.
- 45–75: R plots: PC1–PC2, PC1–PC3, scree plot, coloured by population.
- 75–90: Interpretation in pairs.

**Interpretation questions**
- What does PC1 separate? PC2? Which populations lie between others?
- Where do MXL individuals fall, and why the spread?
- What happens if you drop one population and rerun? (PCs change meaning.)

**Deliverable/homework:** `02_pca.sh` + `plot_pca.R`; a paragraph interpreting the plot.

## Session 6: ADMIXTURE, concepts and run

**Learning objectives**
- Model intuition: each genome is a mixture of K ancestral populations; outputs are the Q matrix (proportions) and P matrix (allele frequencies).
- K is a modelling choice; the "ancestral populations" are statistical constructs, not real historical groups.
- Run `admixture --cv` for K = 2…5; choose K using cross-validation error *and* biological sense.
- Stacked barplots; cluster labels/colours are arbitrary between runs.

**Plan**
- 0–25: Concept (mixture model, K; relation to PCA: same structure, different representation).
- 25–55: Run K=2..5 in a bash `for` loop with `-j` set to 4 at most; the data are small enough that runs take seconds to a few minutes.
- 55–80: Extract CV error with `grep`, plot CV vs K, barplot of Q sorted by population.
- 80–90: Discussion: what does K=3 show for MXL? What do K=2 and K=4 add?

**Deliverable/homework:** `03_admixture.sh` + `plot_admixture.R`; barplots for K=2–5.

## Session 7: Integration, pitfalls and mini project

**Learning objectives**
- Combine PCA and ADMIXTURE: agreements, differences, limits.
- Recognise pitfalls: unbalanced sampling, relatives, too few SNPs, no LD pruning, over-interpreting K, reading admixture as history or dates.
- Reproduce the whole pipeline with one script.

**Plan**
- 0–10: Revisit session 3 predictions.
- 10–50: Mini project in pairs on a "trap" variant of the data: rerun PCA and ADMIXTURE, document what changed and why.
- 50–80: 5-minute presentations per pair: one figure, three sentences of interpretation, one caveat.
- 80–90: Wrap-up: what you can and cannot conclude; next steps (local ancestry, f-statistics, reference panels, ethics of ancestry inference).

**Deliverable:** final 1–2 page report with PCA and ADMIXTURE figures and an interpretation.

---

## Cross-cutting notes

- **Pacing risks:** session 1 installs and session 6 runtime are the most likely overruns. Keep a precomputed-results folder so a stuck student can continue.
- **WSL tips:** work in the Linux home (`~/`), not `/mnt/c`, for speed; copy figures to Windows to view them (`explorer.exe .` opens the current folder).
- **Assessment (light):** small deliverable per session + final report; participation in interpretation discussions.
- **Ethics thread (2–3 short moments, not a separate session):** ancestry vs race, sampling and consent in genomic datasets, limits of reference panels.

## Materials to build (suggested order)

1. Session 1 slides and `setup/` (install steps and `check_setup.sh`).
2. `data/`: small dataset and README documenting provenance (later).
3. `sessions/01…07/`: notes, exercises, solution scripts, cheat sheets.
4. `results_reference/`: instructor's precomputed PCA/ADMIXTURE outputs.
