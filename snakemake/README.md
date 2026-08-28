# ctDNA Somatic Variant Calling Pipeline

Somatic variant calling pipeline for circulating tumor DNA (ctDNA) from targeted panel sequencing, implemented with Snakemake.

Tumor-normal paired design: cell-free DNA (cfDNA) from plasma as tumor input, matched buffy coat as germline normal. Designed for ultra-sensitive detection of low-frequency somatic variants for minimal residual disease (MRD) monitoring.

---

## Pipeline Overview

```
FASTQ_merged_Tumours/  (*_T[1-5]_R1/R2.fastq.gz)
  ──→ ConsensusCruncher ──→ consensus_output_Tumours/{sample}/sscs_sc/*.sorted.sscs.sc.sorted.bam
                                                              dcs_sc/*.sorted.dcs.sc.sorted.bam

FASTQ_merged_Normals/  (*_GL_R1/R2.fastq.gz)
  ──→ ConsensusCruncher ──→ consensus_output_Normals/{normal}/sscs_sc/*_GL.sorted.sscs.sc.sorted.bam
                                                              dcs_sc/*_GL.sorted.dcs.sc.sorted.bam

tumour_normal_SNV_calling/{bam_type}/
  ├── tumor_samples/       (symlinks to consensus_output_Tumours)
  └── germline_controls/   (symlinks to consensus_output_Normals)
         │
         ▼
      Mutect2 ──→ GetPileupSummaries ──→ CalculateContamination
              ──→ FilterMutectCalls  ──→ SelectVariants (PASS, tumor-only)

Filtered VCF ──→ VEP
             ──→ vt normalize/decompose ──→ ANNOVAR
```

Run separately for `sscs_sc` and `dcs_sc` BAM types.

---

## Steps

1. **ConsensusCruncher** — UMI-based consensus for tumor and normal FASTQs separately (fastq2bam → consensus); outputs sscs_sc and dcs_sc BAMs
2. **BedToIntervalList** — Convert amplicon and target BED files to interval list format
3. **CollectTargetedPcrMetrics** — Per-sample coverage metrics
4. **Symlinks** — Organise consensus BAMs into `tumour_normal_SNV_calling/{bam_type}/tumor_samples/` and `germline_controls/` for Mutect2 input
5. **Mutect2** — Tumor-normal somatic variant calling (run for sscs_sc and dcs_sc independently)
6. **Contamination filtering** — GetPileupSummaries → CalculateContamination → FilterMutectCalls
7. **SelectVariants** — Retain PASS variants; extract tumor sample only
8. **VEP** — Functional annotation (variant class, SIFT, PolyPhen, regulatory)
9. **vt + ANNOVAR** — Normalize, decompose, and annotate with population frequency and clinical databases

---

## Requirements

- Snakemake ≥ 5.20
- GATK 4.1.8.1
- ConsensusCruncher
- BWA 0.7.15, samtools 1.20, Python 3.7
- vt 0.577
- VEP 98 (offline cache, GRCh38)
- ANNOVAR (humandb: ensGene, gnomad211_genome, cosmic70, avsnp147, clinvar)

---

## Usage

**1. Configure paths**

Edit `config/config.yaml` to set paths to the reference genome, FASTQ input directories, target BED files, tool executables, and annotation databases.

**2. Prepare FASTQ directories**

Place merged FASTQ files in two directories:

```
FASTQ_merged_Tumours/
  ├── SAMPLE001_T1_R1.fastq.gz
  ├── SAMPLE001_T1_R2.fastq.gz
  ├── SAMPLE001_T2_R1.fastq.gz
  └── SAMPLE001_T2_R2.fastq.gz

FASTQ_merged_Normals/
  ├── SAMPLE001_GL_R1.fastq.gz
  └── SAMPLE001_GL_R2.fastq.gz
```

Tumor files must follow the pattern `*_T[1-5]_R1/R2.fastq.gz`. Normal files must follow `*_GL_R1/R2.fastq.gz`. Pairing is automatic: `SAMPLE001_T1` and `SAMPLE001_T2` are both matched to `SAMPLE001_GL`.

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
├── consensus_Tumours/               # Tumor fastq2bam BAMs
├── consensus_output_Tumours/        # Tumor SSCS/DCS consensus BAMs
│   └── {sample}/
│       ├── sscs_sc/
│       └── dcs_sc/
├── consensus_Normals/               # Normal fastq2bam BAMs
├── consensus_output_Normals/        # Normal SSCS/DCS consensus BAMs
│   └── {normal_prefix}/
│       ├── sscs_sc/
│       └── dcs_sc/
├── tumour_normal_SNV_calling/       # Symlinked BAMs organised for Mutect2
│   ├── sscs_sc/
│   │   ├── tumor_samples/
│   │   └── germline_controls/
│   └── dcs_sc/
│       ├── tumor_samples/
│       └── germline_controls/
├── coverage/                        # Interval lists and PCR metrics
├── vcfs/{bam_type}/                 # Raw, filtered, and final VCFs
├── vep/{bam_type}/                  # VEP annotation tables
└── annovar/{bam_type}/              # Normalized, decomposed, ANNOVAR-annotated VCFs
```

---

## Notes

- Reference genome: GRCh38 (no-alt analysis set)
- ANNOVAR databases: ensGene, gnomAD 2.1.1, COSMIC 70, dbSNP 147, ClinVar
- Designed for HPC environments with Slurm
