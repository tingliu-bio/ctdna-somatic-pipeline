# ctDNA Somatic Variant Calling Pipeline

Somatic variant calling pipeline for circulating tumor DNA (ctDNA) from targeted panel sequencing, implemented with Snakemake.

Tumor-normal paired design: cell-free DNA (cfDNA) from plasma as tumor input, matched buffy coat as germline normal. Designed for ultra-sensitive detection of low-frequency somatic variants for minimal residual disease (MRD) monitoring.

---

## Pipeline Overview

```
FASTQ_merged_Tumours/ ──→ ConsensusCruncher ──→ consensus_output_Tumours/
                                                   ├── sscs_sc/*.sorted.sscs.sc.sorted.bam
                                                   └── dcs_sc/*.sorted.dcs.sc.sorted.bam

FASTQ_merged_Normals/ ──→ ConsensusCruncher ──→ consensus_output_Normals/
                                                   ├── sscs_sc/*_GL.sorted.sscs.sc.sorted.bam
                                                   └── dcs_sc/*_GL.sorted.dcs.sc.sorted.bam

consensus_output_Tumours ──┐
                           ├─→ Mutect2 ──→ GetPileupSummaries ──→ CalculateContamination
consensus_output_Tumours  ─┘           ──→ FilterMutectCalls ──→ SelectVariants (PASS,tumor-only)

Filtered VCF ──→ VEP
             ──→ vt normalize/decompose ──→ ANNOVAR
```

Run separately for `sscs_sc` and `dcs_sc` BAM types.

## Steps

1. **ConsensusCruncher** — UMI-based consensus for tumor and normal FASTQs separately (fastq2bam → consensus); outputs sscs_sc and dcs_sc BAMs
2. **BedToIntervalList** — Convert amplicon and target BED files to interval list format
3. **CollectTargetedPcrMetrics** — Per-sample coverage metrics
4. **Mutect2** — Tumor-normal somatic variant calling (run for sscs_sc and dcs_sc independently)
5. **Contamination filtering** — GetPileupSummaries → CalculateContamination → FilterMutectCalls
6. **SelectVariants** — Retain PASS variants; extract tumor sample only
7. **VEP** — Functional annotation (variant class, SIFT, PolyPhen, regulatory)
8. **vt + ANNOVAR** — Normalize, decompose, and annotate with population frequency and clinical databases

---

## Requirements

- Snakemake ≥ 5.20
- GATK 4.x
- ConsensusCruncher
- samtools, BWA
- vt
- VEP (offline cache, GRCh38)
- ANNOVAR (humandb: ensGene, gnomad211_genome, cosmic70, avsnp147, clinvar)

---

## Usage

**1. Configure paths**

Edit `config/config.yaml` to set paths to the reference genome, target BED files, tool executables, and annotation databases.

**2. Prepare sample table**

Edit `config/samples.tsv`:

| sample | fq1 | fq2 | normal_bam | tumor_sm | normal_sm |
|--------|-----|-----|------------|----------|-----------|
| SAMPLE1 | SAMPLE1_R1.fastq.gz | SAMPLE1_R2.fastq.gz | SAMPLE1_GL.bam | SAMPLE1_T | SAMPLE1_GL |

**3. Run**

```bash
# Dry run
snakemake -s workflow/Snakefile -n

# Submit to Slurm cluster
sbatch run_pipeline.sh
```

---

## Output

```
results/
├── consensus/          # ConsensusCruncher BAMs
├── consensus_output/   # DCS consensus BAMs
├── coverage/           # Interval lists and PCR metrics
├── vcfs/               # Raw, filtered, and final VCFs
├── vep/                # VEP annotation tables
└── annovar/            # Normalized, decomposed, ANNOVAR-annotated VCFs
```

---

## Notes

- Reference genome: GRCh38 (no-alt analysis set)
- ANNOVAR databases: ensGene, gnomAD 2.1.1, COSMIC 70, dbSNP 147, ClinVar
- Designed for HPC environments with Slurm
