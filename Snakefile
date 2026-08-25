import pandas as pd
from os.path import join

configfile: "config/config.yaml"

# Sample table columns:
#   sample        - sample ID (tumor prefix)
#   tumor_fq1     - tumor R1 FASTQ
#   tumor_fq2     - tumor R2 FASTQ
#   normal_fq1    - normal (GL) R1 FASTQ
#   normal_fq2    - normal (GL) R2 FASTQ
#   normal_prefix - normal BAM prefix (used to locate GL BAM after ConsensusCruncher)
samples = pd.read_csv(config["sample_table"], sep="\t", index_col="sample")
SAMPLES = samples.index.tolist()

BAM_TYPES = ["sscs_sc", "dcs_sc"]

wildcard_constraints:
    sample    = r"[\w\-]+",
    bam_type  = r"sscs_sc|dcs_sc",

### Target rule ----------------------------------------------------------------

rule all:
    input:
        expand(join(config["outdir"], "annovar/{bam_type}/{sample}.annot.norm_decomp_filtered.hg38_multianno.vcf"),
               sample=SAMPLES, bam_type=BAM_TYPES),
        expand(join(config["outdir"], "vep/{bam_type}/{sample}.vep_annot.txt"),
               sample=SAMPLES, bam_type=BAM_TYPES),
        expand(join(config["outdir"], "coverage/{sample}.output_pcr_metrics.txt"),
               sample=SAMPLES),

### Step 1a: ConsensusCruncher — Tumor FASTQ to BAM ---------------------------

rule tumor_fastq2bam:
    input:
        fq1 = lambda wc: samples.loc[wc.sample, "tumor_fq1"],
        fq2 = lambda wc: samples.loc[wc.sample, "tumor_fq2"],
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

### Step 1b: ConsensusCruncher — Tumor consensus (UMI deduplication) -----------

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
        python {config[consensus_cruncher]} consensus \
            -i {input.bam} \
            -o {params.outdir} \
            -s {params.samtools} \
            -g hg38 \
            -b {params.cytoband} \
            --cleanup True
        """

### Step 2a: ConsensusCruncher — Normal FASTQ to BAM --------------------------

rule normal_fastq2bam:
    input:
        fq1 = lambda wc: samples.loc[wc.sample, "normal_fq1"],
        fq2 = lambda wc: samples.loc[wc.sample, "normal_fq2"],
    output:
        bam = lambda wc: join(config["outdir"],
                              f"consensus_Normals/{wc.sample}/bamfiles/{samples.loc[wc.sample, 'normal_prefix']}.sorted.bam"),
    params:
        outdir   = join(config["outdir"], "consensus_Normals/{sample}"),
        ref      = config["ref_genome"],
        umi      = config["umi_index"],
        bwa      = config["bwa_path"],
        samtools = config["samtools_path"],
        picard   = config["picard_jar"],
    message: "Normal fastq2bam: {wildcards.sample}"
    shell:
        """
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
        sscs = lambda wc: join(config["outdir"],
                               f"consensus_output_Normals/{wc.sample}/sscs_sc/{samples.loc[wc.sample, 'normal_prefix']}.sorted.sscs.sc.sorted.bam"),
        dcs  = lambda wc: join(config["outdir"],
                               f"consensus_output_Normals/{wc.sample}/dcs_sc/{samples.loc[wc.sample, 'normal_prefix']}.sorted.dcs.sc.sorted.bam"),
    params:
        outdir   = join(config["outdir"], "consensus_output_Normals/{sample}"),
        samtools = config["samtools_path"],
        cytoband = config["cytoband"],
    message: "Normal consensus: {wildcards.sample}"
    shell:
        """
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

### Step 4: Mutect2 tumor-normal (sscs_sc and dcs_sc) -------------------------

def get_tumor_bam(wc):
    if wc.bam_type == "sscs_sc":
        return join(config["outdir"], f"consensus_output_Tumours/{wc.sample}/sscs_sc/{wc.sample}.sorted.sscs.sc.sorted.bam")
    else:
        return join(config["outdir"], f"consensus_output_Tumours/{wc.sample}/dcs_sc/{wc.sample}.sorted.dcs.sc.sorted.bam")

def get_normal_bam(wc):
    np = samples.loc[wc.sample, "normal_prefix"]
    if wc.bam_type == "sscs_sc":
        return join(config["outdir"], f"consensus_output_Normals/{wc.sample}/sscs_sc/{np}.sorted.sscs.sc.sorted.bam")
    else:
        return join(config["outdir"], f"consensus_output_Normals/{wc.sample}/dcs_sc/{np}.sorted.dcs.sc.sorted.bam")

rule mutect2:
    input:
        tumor  = get_tumor_bam,
        normal = get_normal_bam,
        tgt    = rules.bed_to_interval_list.output.tgt,
    output:
        vcf = join(config["outdir"], "vcfs/{bam_type}/{sample}.raw.vcf.gz"),
    message: "Mutect2 [{wildcards.bam_type}]: {wildcards.sample}"
    shell:
        """
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

### Step 5: Contamination estimation and filtering -----------------------------

rule filter_variants:
    input:
        bam  = get_tumor_bam,
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

### Step 6: VEP annotation -----------------------------------------------------

rule vep_annotate:
    input:
        vcf = rules.filter_variants.output.final,
    output:
        txt = join(config["outdir"], "vep/{bam_type}/{sample}.vep_annot.txt"),
    message: "VEP [{wildcards.bam_type}]: {wildcards.sample}"
    shell:
        """
        vep --offline \
            --dir_cache {config[vep_cache]} \
            -i {input.vcf} \
            --force_overwrite \
            --variant_class --sift b --polyphen b \
            --nearest symbol --gene_phenotype --regulatory --show_ref_allele \
            --tab \
            -o {output.txt}
        """

### Step 7: Normalize, decompose, ANNOVAR --------------------------------------

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
