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
consensus_output_Normals  ─┘           ──→ FilterMutectCalls ──→ SelectVariants (PASS, tumor-only)

Filtered VCF ──→ VEP
             ──→ vt normalize/decompose ──→ ANNOVAR
```

Run separately for `sscs_sc` and `dcs_sc` BAM types.

---

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

| sample | tumor_fq1 | tumor_fq2 | normal_fq1 | normal_fq2 | normal_prefix |
|--------|-----------|-----------|------------|------------|---------------|
| SAMPLE001_T | SAMPLE001_T_R1.fastq.gz | SAMPLE001_T_R2.fastq.gz | SAMPLE001_GL_R1.fastq.gz | SAMPLE001_GL_R2.fastq.gz | SAMPLE001_GL |

**3. Run**

```bash
# Dry run
snakemake -s workflow/ConsensusCruncher_tumor_normal_hg38.smk -n

# Run locally
snakemake -s workflow/ConsensusCruncher_tumor_normal_hg38.smk --cores 8

# Submit to Slurm cluster
sbatch run_pipeline.sh
```

---

## Output

```
{outdir}/
├── consensus_Tumours/          # Tumor fastq2bam BAMs
├── consensus_output_Tumours/   # Tumor SSCS/DCS consensus BAMs
│   ├── sscs_sc/
│   └── dcs_sc/
├── consensus_Normals/          # Normal fastq2bam BAMs
├── consensus_output_Normals/   # Normal SSCS/DCS consensus BAMs
│   ├── sscs_sc/
│   └── dcs_sc/
├── coverage/                   # Interval lists and PCR metrics
├── vcfs/{bam_type}/            # Raw, filtered, and final VCFs
├── vep/{bam_type}/             # VEP annotation tables
└── annovar/{bam_type}/         # Normalized, decomposed, ANNOVAR-annotated VCFs
```

---

## Notes

- Reference genome: GRCh38 (no-alt analysis set)
- ANNOVAR databases: ensGene, gnomAD 2.1.1, COSMIC 70, dbSNP 147, ClinVar
- Designed for HPC environments with Slurm
