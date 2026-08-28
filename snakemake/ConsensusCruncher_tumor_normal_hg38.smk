import re
from os.path import join, basename
from glob import glob

configfile: "config/config.yaml"

# ── Sample discovery from FASTQ directories ───────────────────────────────────
# Tumors: *_T[1-5]_R1.fastq.gz
tumor_fq1s = sorted(glob(join(config["tumor_fastq_dir"], "*_T[1-5]_R1.fastq.gz")))
if not tumor_fq1s:
    raise ValueError(f"No tumor FASTQs (*_T[1-5]_R1.fastq.gz) found in {config['tumor_fastq_dir']}")
TUMOR_SAMPLES = [re.sub(r'_R1\.fastq\.gz$', '', basename(f)) for f in tumor_fq1s]

# Normals: *_GL_R1.fastq.gz
normal_fq1s = sorted(glob(join(config["normal_fastq_dir"], "*_GL_R1.fastq.gz")))
if not normal_fq1s:
    raise ValueError(f"No normal FASTQs (*_GL_R1.fastq.gz) found in {config['normal_fastq_dir']}")
NORMAL_PREFIXES = [re.sub(r'_R1\.fastq\.gz$', '', basename(f)) for f in normal_fq1s]

# Pair each tumor with its matched normal: strip _T[1-5], append _GL
def get_normal_prefix(tumor_sample):
    base = re.sub(r'_T[1-5]$', '', tumor_sample)
    match = f"{base}_GL"
    if match not in NORMAL_PREFIXES:
        raise ValueError(
            f"No matching normal for '{tumor_sample}': "
            f"expected '{match}_R1.fastq.gz' in {config['normal_fastq_dir']}"
        )
    return match

TUMOR_TO_NORMAL = {s: get_normal_prefix(s) for s in TUMOR_SAMPLES}
UNIQUE_NORMALS  = sorted(set(TUMOR_TO_NORMAL.values()))

BAM_TYPES = ["sscs_sc", "dcs_sc"]

wildcard_constraints:
    sample        = r"[\w\-]+",
    bam_type      = r"sscs_sc|dcs_sc",
    normal_prefix = r"[\w\-]+",

def bam_dot(bam_type):
    """sscs_sc → sscs.sc,  dcs_sc → dcs.sc (ConsensusCruncher filename convention)."""
    return bam_type.replace("_", ".")

### Target rule ----------------------------------------------------------------

rule all:
    input:
        expand(join(config["outdir"], "annovar/{bam_type}/{sample}.annot.norm_decomp_filtered.hg38_multianno.vcf"),
               sample=TUMOR_SAMPLES, bam_type=BAM_TYPES),
        expand(join(config["outdir"], "vep/{bam_type}/{sample}.vep_annot.txt"),
               sample=TUMOR_SAMPLES, bam_type=BAM_TYPES),
        expand(join(config["outdir"], "coverage/{sample}.output_pcr_metrics.txt"),
               sample=TUMOR_SAMPLES),

### Step 1a: ConsensusCruncher — Tumor FASTQ to BAM ---------------------------

rule tumor_fastq2bam:
    input:
        fq1 = lambda wc: join(config["tumor_fastq_dir"], f"{wc.sample}_R1.fastq.gz"),
        fq2 = lambda wc: join(config["tumor_fastq_dir"], f"{wc.sample}_R2.fastq.gz"),
    output:
        bam = join(config["outdir"], "consensus_Tumours/{sample}/bamfiles/{sample}.sorted.bam"),
    params:
        outdir   = join(config["outdir"], "consensus_Tumours/{sample}"),
        ref      = config["ref_genome"],
        umi      = config["umi_index"],
        bwa      = config["bwa_path"],
        samtools = config["samtools_path"],
        picard   = config["picard_jar"],
    message: "Tumor fastq2bam: {wildcards.sample}"
    shell:
        """
        module load bwa/0.7.15
        module load samtools/1.20
        module load python3/3.7.2
        python {config[consensus_cruncher]} fastq2bam \
            --fastq1 {input.fq1} \
            --fastq2 {input.fq2} \
            -o {params.outdir} \
            -r {params.ref} \
            -b {params.bwa} \
            -s {params.samtools} \
            -l {params.umi} \
            -g {params.picard}
        """

### Step 1b: ConsensusCruncher — Tumor consensus (UMI deduplication) ----------

rule tumor_consensus:
    input:
        bam = rules.tumor_fastq2bam.output.bam,
    output:
        sscs = join(config["outdir"], "consensus_output_Tumours/{sample}/sscs_sc/{sample}.sorted.sscs.sc.sorted.bam"),
        dcs  = join(config["outdir"], "consensus_output_Tumours/{sample}/dcs_sc/{sample}.sorted.dcs.sc.sorted.bam"),
    params:
        outdir   = join(config["outdir"], "consensus_output_Tumours/{sample}"),
        samtools = config["samtools_path"],
        cytoband = config["cytoband"],
    message: "Tumor consensus: {wildcards.sample}"
    shell:
        """
        module load samtools/1.20
        module load python3/3.7.2
        python {config[consensus_cruncher]} consensus \
            -i {input.bam} \
            -o {params.outdir} \
            -s {params.samtools} \
            -g hg38 \
            -b {params.cytoband} \
            --cleanup True
        """

### Step 2a: ConsensusCruncher — Normal FASTQ to BAM --------------------------
# Each unique normal is processed once, keyed by normal_prefix (e.g. SAMPLE001_GL)

rule normal_fastq2bam:
    input:
        fq1 = lambda wc: join(config["normal_fastq_dir"], f"{wc.normal_prefix}_R1.fastq.gz"),
        fq2 = lambda wc: join(config["normal_fastq_dir"], f"{wc.normal_prefix}_R2.fastq.gz"),
    output:
        bam = join(config["outdir"], "consensus_Normals/{normal_prefix}/bamfiles/{normal_prefix}.sorted.bam"),
    params:
        outdir   = join(config["outdir"], "consensus_Normals/{normal_prefix}"),
        ref      = config["ref_genome"],
        umi      = config["umi_index"],
        bwa      = config["bwa_path"],
        samtools = config["samtools_path"],
        picard   = config["picard_jar"],
    message: "Normal fastq2bam: {wildcards.normal_prefix}"
    shell:
        """
        module load bwa/0.7.15
        module load samtools/1.20
        module load python3/3.7.2
        python {config[consensus_cruncher]} fastq2bam \
            --fastq1 {input.fq1} \
            --fastq2 {input.fq2} \
            -o {params.outdir} \
            -r {params.ref} \
            -b {params.bwa} \
            -s {params.samtools} \
            -l {params.umi} \
            -g {params.picard}
        """

### Step 2b: ConsensusCruncher — Normal consensus ------------------------------

rule normal_consensus:
    input:
        bam = rules.normal_fastq2bam.output.bam,
    output:
        sscs = join(config["outdir"], "consensus_output_Normals/{normal_prefix}/sscs_sc/{normal_prefix}.sorted.sscs.sc.sorted.bam"),
        dcs  = join(config["outdir"], "consensus_output_Normals/{normal_prefix}/dcs_sc/{normal_prefix}.sorted.dcs.sc.sorted.bam"),
    params:
        outdir   = join(config["outdir"], "consensus_output_Normals/{normal_prefix}"),
        samtools = config["samtools_path"],
        cytoband = config["cytoband"],
    message: "Normal consensus: {wildcards.normal_prefix}"
    shell:
        """
        module load samtools/1.20
        module load python3/3.7.2
        python {config[consensus_cruncher]} consensus \
            -i {input.bam} \
            -o {params.outdir} \
            -s {params.samtools} \
            -g hg38 \
            -b {params.cytoband} \
            --cleanup True
        """

### Step 3: Interval lists and coverage metrics --------------------------------

rule bed_to_interval_list:
    input:
        bam = join(config["outdir"], "consensus_output_Tumours/{sample}/sscs_sc/{sample}.sorted.sscs.sc.sorted.bam"),
    output:
        amp = join(config["outdir"], "coverage/{sample}_amplicon.interval_list"),
        tgt = join(config["outdir"], "coverage/{sample}_targets.interval_list"),
    message: "BedToIntervalList: {wildcards.sample}"
    shell:
        """
        gatk BedToIntervalList \
            -I {config[amplicon_bed]} \
            -O {output.amp} \
            -SD {input.bam}

        gatk BedToIntervalList \
            -I {config[targets_bed]} \
            -O {output.tgt} \
            -SD {input.bam}
        """

rule coverage_metrics:
    input:
        bam = join(config["outdir"], "consensus_output_Tumours/{sample}/sscs_sc/{sample}.sorted.sscs.sc.sorted.bam"),
        amp = rules.bed_to_interval_list.output.amp,
        tgt = rules.bed_to_interval_list.output.tgt,
    output:
        metrics = join(config["outdir"], "coverage/{sample}.output_pcr_metrics.txt"),
        per_tgt = join(config["outdir"], "coverage/{sample}.per_target_coverage.txt"),
    message: "CollectTargetedPcrMetrics: {wildcards.sample}"
    shell:
        """
        gatk CollectTargetedPcrMetrics \
            -I {input.bam} \
            -O {output.metrics} \
            -R {config[ref_genome]} \
            --PER_TARGET_COVERAGE {output.per_tgt} \
            --AMPLICON_INTERVALS {input.amp} \
            --TARGET_INTERVALS {input.tgt}
        """

### Step 4: Symlinks — organise BAMs for Mutect2 ------------------------------
# Creates tumour_normal_SNV_calling/{bam_type}/tumor_samples/ and germline_controls/

rule link_tumor_bam:
    input:
        bam = lambda wc: join(config["outdir"],
            f"consensus_output_Tumours/{wc.sample}/{wc.bam_type}"
            f"/{wc.sample}.sorted.{bam_dot(wc.bam_type)}.sorted.bam"),
    output:
        bam = join(config["outdir"],
            "tumour_normal_SNV_calling/{bam_type}/tumor_samples"
            "/{sample}.sorted.{bam_type}.sorted.bam"),
    message: "Linking tumor BAM [{wildcards.bam_type}]: {wildcards.sample}"
    shell:
        """
        ln -sf $(realpath {input.bam}) {output.bam}
        ln -sf $(realpath {input.bam}).bai {output.bam}.bai
        """

rule link_normal_bam:
    input:
        bam = lambda wc: join(config["outdir"],
            f"consensus_output_Normals/{wc.normal_prefix}/{wc.bam_type}"
            f"/{wc.normal_prefix}.sorted.{bam_dot(wc.bam_type)}.sorted.bam"),
    output:
        bam = join(config["outdir"],
            "tumour_normal_SNV_calling/{bam_type}/germline_controls"
            "/{normal_prefix}.sorted.{bam_type}.sorted.bam"),
    message: "Linking normal BAM [{wildcards.bam_type}]: {wildcards.normal_prefix}"
    shell:
        """
        ln -sf $(realpath {input.bam}) {output.bam}
        ln -sf $(realpath {input.bam}).bai {output.bam}.bai
        """

### Step 5: Mutect2 tumor-normal ----------------------------------------------

rule mutect2:
    input:
        tumor  = lambda wc: join(config["outdir"],
            f"tumour_normal_SNV_calling/{wc.bam_type}/tumor_samples"
            f"/{wc.sample}.sorted.{wc.bam_type}.sorted.bam"),
        normal = lambda wc: join(config["outdir"],
            f"tumour_normal_SNV_calling/{wc.bam_type}/germline_controls"
            f"/{TUMOR_TO_NORMAL[wc.sample]}.sorted.{wc.bam_type}.sorted.bam"),
        tgt    = rules.bed_to_interval_list.output.tgt,
    output:
        vcf = join(config["outdir"], "vcfs/{bam_type}/{sample}.raw.vcf.gz"),
    message: "Mutect2 [{wildcards.bam_type}]: {wildcards.sample}"
    shell:
        """
        module load samtools/1.20
        module load gatk/4.1.8.1
        tumor_sm=$(samtools view -H {input.tumor} | grep '^@RG' | head -1 | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g")
        normal_sm=$(samtools view -H {input.normal} | grep '^@RG' | head -1 | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g")

        gatk Mutect2 \
            -R {config[ref_genome]} \
            -I {input.tumor} \
            -I {input.normal} \
            -tumor $tumor_sm \
            -normal $normal_sm \
            -L {input.tgt} \
            -O {output.vcf} \
            --germline-resource {config[gnomad_vcf]}
        """

### Step 6: Contamination estimation and filtering ----------------------------

rule filter_variants:
    input:
        bam  = lambda wc: join(config["outdir"],
            f"tumour_normal_SNV_calling/{wc.bam_type}/tumor_samples"
            f"/{wc.sample}.sorted.{wc.bam_type}.sorted.bam"),
        vcf  = rules.mutect2.output.vcf,
        tgt  = rules.bed_to_interval_list.output.tgt,
    output:
        pileup = join(config["outdir"], "vcfs/{bam_type}/{sample}.pileups.table"),
        contam = join(config["outdir"], "vcfs/{bam_type}/{sample}.contamination.table"),
        unfilt = join(config["outdir"], "vcfs/{bam_type}/{sample}.unfiltered.vcf.gz"),
        filt   = join(config["outdir"], "vcfs/{bam_type}/{sample}.filtered.vcf.gz"),
        final  = join(config["outdir"], "vcfs/{bam_type}/{sample}.filtered.case_only.vcf.gz"),
    message: "Filtering [{wildcards.bam_type}]: {wildcards.sample}"
    shell:
        """
        module load samtools/1.20
        module load gatk/4.1.8.1
        tumor_sm=$(samtools view -H {input.bam} | grep '^@RG' | head -1 | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g")

        gatk GetPileupSummaries \
            -I {input.bam} \
            -V {config[gnomad_vcf]} \
            -L {input.tgt} \
            -O {output.pileup}

        gatk CalculateContamination \
            -I {output.pileup} \
            -O {output.contam}

        gatk FilterMutectCalls \
            -R {config[ref_genome]} \
            -V {input.vcf} \
            --contamination-table {output.contam} \
            -O {output.unfilt}

        gatk SelectVariants \
            -R {config[ref_genome]} \
            -V {output.unfilt} \
            --exclude-filtered \
            -O {output.filt}

        gatk SelectVariants \
            -R {config[ref_genome]} \
            -V {output.filt} \
            -sn $tumor_sm \
            --exclude-non-variants true \
            -O {output.final}
        """

### Step 7: VEP annotation ----------------------------------------------------

rule vep_annotate:
    input:
        vcf = rules.filter_variants.output.final,
    output:
        txt = join(config["outdir"], "vep/{bam_type}/{sample}.vep_annot.txt"),
    message: "VEP [{wildcards.bam_type}]: {wildcards.sample}"
    shell:
        """
        module load vep/98
        vep --offline \
            --dir_cache {config[vep_cache]} \
            -i {input.vcf} \
            --force_overwrite \
            --variant_class --sift b --polyphen b \
            --nearest symbol --gene_phenotype --regulatory --show_ref_allele \
            --tab \
            -o {output.txt}
        """

### Step 8: Normalize, decompose, ANNOVAR -------------------------------------

rule annovar_annotate:
    input:
        vcf = rules.filter_variants.output.final,
    output:
        norm   = join(config["outdir"], "annovar/{bam_type}/{sample}.filtered_norm.vcf.gz"),
        decomp = join(config["outdir"], "annovar/{bam_type}/{sample}.filtered_norm_decomp.vcf.gz"),
        annot  = join(config["outdir"], "annovar/{bam_type}/{sample}.annot.norm_decomp_filtered.hg38_multianno.vcf"),
    params:
        outprefix = join(config["outdir"], "annovar/{bam_type}/{sample}.annot.norm_decomp_filtered"),
    message: "ANNOVAR [{wildcards.bam_type}]: {wildcards.sample}"
    shell:
        """
        module load vt/0.577
        module load annovar/20180416
        module load tabix
        vt normalize {input.vcf} \
            -r {config[ref_genome]} \
            -o {output.norm}

        vt decompose -s {output.norm} \
            -o {output.decomp}

        table_annovar.pl \
            --buildver hg38 \
            {output.decomp} \
            {config[annovar_db]} \
            --protocol ensGene,gnomad211_genome,cosmic70,avsnp147,clinvar_20190305 \
            --operation g,f,f,f,f \
            --outfile {params.outprefix} \
            --onetranscript \
            --vcfinput
        """
