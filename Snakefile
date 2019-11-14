import pandas as pd


####################### shell setup #################

shell.prefix("set -eo pipefail; echo BEGIN at $(date);  ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

####################### external configs ############
configfile: "config.yaml"
workdir: "res/"
# for a local test, will be deleted for cloud mode
localrules: all


####################### global parameters ###########
# input is paired uuid table from GDC api.
df = pd.read_csv(config["samples"], sep = "\t")
PID = (df.cohort + "_" + df.patient).tolist()
df["pid"] = PID
TYPES=["tumor", "normal"]
# the subintervals to scatter on
SUBINTS = ['{:0>4}'.format(i) for i in range(config["pieces"])]
INTERVALS = ["{}/{}-scattered.interval_list".format(config["subint_dir"],int1) for int1 in SUBINTS]


# all the target files
rule all:
    input: 
        expand("{pid}/funco_flag", pid=PID),
        annot_merged_filtered_maf=expand("{pid}/annot_merged_filtered.vcf", pid=PID),
        merged_filtered_vcf=expand("{pid}/merged_filtered.vcf", pid=PID)


# set uo the file structure
rule fs_setup:
    params:
        pid=PID
    output:
        touch("{pid}/fs_flag")
    shell:
        """
        mkdir -p {params.pid}/subvcfs
        mkdir -p {params.pid}/tpile_dir
        mkdir -p {params.pid}/npile_dir
        mkdir -p {params.pid}/f1r2
        touch {params.pid}/fs_flag
        """
######################## preparation #######################

# split intervals for scattering jobs
rule split_intervals:
    params:
        ref=config["ref"],
        pieces=config["pieces"],
        wes=config["wes_interval_list"]
    output:
        directory(config["subint_dir"])
    log:
        "logs/split_interval/log"
    shell:
        """
        gatk SplitIntervals \
            -R {params.ref} -L {params.wes} \
            --scatter-count {params.pieces} \
            -O {output[0]} &> {log}
        """

# get signed url from nci-data-commons
rule get_url:
    input:
        rules.fs_setup.output,
        rules.split_intervals.output
    params:
        tumor=lambda wildcards: df.loc[df["pid"] == wildcards.pid, "tumor"].item(),
        normal=lambda wildcards: df.loc[df["pid"] == wildcards.pid,"normal"].item(),
        credential = config["credential"]
    output:
        expand("{{pid}}/{type}_url.txt", type=TYPES)
    log:
        "logs/get_url/{pid}.log"
    script:
        "scripts/localize_from_nci.py"


# localize the signed urls, annotate with sample name
rule localize_with_name:
    input:
        url="{pid}/{type}_url.txt"
    priority: 5
    params:
        ref=config["ref"],
        gs_bai="{pid}_{type}.bai",
        bucket=config["index_bucket"]
    output:
        bam=temp("{pid}/{type}.bam"),
        bai="{pid}/{type}.bai",
        sample_name="{pid}/{type}_name.txt"
    log:
        "logs/localize/{pid}_{type}.log"
    message:
        "{wildcards.pid} - {wildcards.type} is being localized and annotated with file name"
    shell:
        """
        cat {input.url} &> {log}
	echo --------------
        wget "$(cat {input.url})" -O {output.bam} &>> {log}
        gatk BuildBamIndex -I {output.bam} -O {output.bai} &>> {log}
        
        [ gsutil -q stat gs://{params.bucket}/{params.gs_bai} ] && gsutil cp {params.gs_bai} gs://{params.bucket}/{params.gs_bai} &>> {log}
        
        gatk GetSampleName -R {params.ref} -I {output.bam} -O {output.sample_name} &>> {log}
        cat {output.sample_name} &>> {log}
        """

# run M2 on tumor-normal pairs for 
rule scatter_m2:
    input:
        expand("{{pid}}/{type}.bam", type=TYPES)
    params:
	#pid=lambda wildcards: wildcards.pid ,
        subint_dir=config["subint_dir"],
        ref=config["ref"],
        vfc=config["vfc"],
        gnomad=config["gnomad"],
        default_pon=config["default_pon"]
    log:
        "logs/scatter/{pid}_{chr}.log"
    output:
        out_tpile="{pid}/tpile_dir/{chr}-tpile.table",
        out_npile="{pid}/npile_dir/{chr}-npile.table",
        out_f1r2="{pid}/f1r2/{chr}-f1r2.tar.gz",
        out_vcf="{pid}/subvcfs/{chr}.vcf",
        out_stats="{pid}/subvcfs/{chr}.vcf.stats"
        
    shell:
        """
        gatkm="gatk --java-options -Xmx2g"

        tumor="{wildcards.pid}/tumor.bam"
        normal="{wildcards.pid}/normal.bam"
        tumor_name="{wildcards.pid}/tumor_name.txt"
        normal_name="{wildcards.pid}/normal_name.txt"
        interval_file="{params.subint_dir}/{wildcards.chr}-scattered.interval_list"

        normal_command_line="-I ${{normal}} -normal `cat ${{normal_name}}`"
        tumor_command_line="-I ${{tumor}} -tumor `cat ${{tumor_name}}`"

        $gatkm Mutect2 \
            -R {params.ref} \
            $tumor_command_line \
            $normal_command_line \
            --germline-resource {params.gnomad} \
            -pon {params.default_pon} \
            -L  $interval_file \
            -O {output.out_vcf} \
            --f1r2-tar-gz {output.out_f1r2}
        
        echo finish Mutect2 --------------------

        $gatkm GetPileupSummaries \
            -R {params.ref} -I $tumor \
            --interval-set-rule INTERSECTION \
            -L $interval_file \
            -V {params.vfc} \
            -L {params.vfc} \
            -O {output.out_tpile}
        echo Getpiles tumor --------------------

        $gatkm GetPileupSummaries \
            -R {params.ref} -I $normal \
            --interval-set-rule INTERSECTION \
            -L $interval_file \
            -V {params.vfc} \
            -L {params.vfc} \
            -O {output.out_npile} 
        echo Getpiles normal ------------------

        """

# a checkpoint
rule check_all_intervals:
    input:
        expand("{{pid}}/tpile_dir/{chr}-tpile.table", chr=SUBINTS),
        expand("{{pid}}/npile_dir/{chr}-npile.table", chr=SUBINTS),
        expand("{{pid}}/f1r2/{chr}-f1r2.tar.gz", chr=SUBINTS),
        expand("{{pid}}/subvcfs/{chr}.vcf", chr=SUBINTS),
        expand("{{pid}}/subvcfs/{chr}.vcf.stats", chr=SUBINTS)
    output:
        touch("{pid}/flag")
    shell:
        """
        echo "scatter finished" > {output}
        """

############### merge interval results ###############
rule calculate_contamination:
    input:
        flag="{pid}/flag",
        tumor_pile_table=expand("{{pid}}/tpile_dir/{chr}-tpile.table", chr=SUBINTS),
        normal_pile_table=expand("{{pid}}/npile_dir/{chr}-npile.table", chr=SUBINTS)
    params:
        heap_mem=4,
        all_tumor_piles_input = lambda wildcards, input: " -I ".join(input.tumor_pile_table),
        all_normal_piles_input = lambda wildcards, input: " -I ".join(input.normal_pile_table),
    log:
        "logs/contam/{pid}.log"
    message:
        "{wildcards.pid} is being estimated for contamination level"
    group: "merge"
    output:
        normal_pile_table="{pid}/normal_pile.tsv",
        tumor_pile_table="{pid}/tumor_pile.tsv",
        contamination_table="{pid}/contamination.table",
        segments_table="{pid}/segments.table"
    shell:
        """
        gatkm="gatk --java-options -Xmx{params.heap_mem}g"
        echo "-----------------------gather pile up summaries-----------------------"
        $gatkm GatherPileupSummaries \
            -I {params.all_normal_piles_input} \
            --sequence-dictionary "/demo-mount/refs/Homo_sapiens_assembly19.dict" \
            -O {output.normal_pile_table} &> {log}

        echo "-----------------------gather tumor piles-----------------------"
        $gatkm GatherPileupSummaries \
            -I {params.all_tumor_piles_input} \
            --sequence-dictionary "/demo-mount/refs/Homo_sapiens_assembly19.dict" \
            -O {output.tumor_pile_table} &>> {log}

        echo "-----------------------calc contamination-----------------------"
        $gatkm CalculateContamination \
            -I {output.tumor_pile_table} \
            -O {output.contamination_table} \
            --tumor-segmentation {output.segments_table} \
            -matched {output.normal_pile_table} &>> {log}
        """



rule merge_m2:
    priority: 10
    input:
        vcf=expand("{{pid}}/subvcfs/{chr}.vcf", chr=SUBINTS),
        f1r2=expand("{{pid}}/f1r2/{chr}-f1r2.tar.gz", chr=SUBINTS),
        stats=expand("{{pid}}/subvcfs/{chr}.vcf.stats", chr=SUBINTS),
        contamination_table="{pid}/contamination.table",
        segments_table="{pid}/segments.table"
    params:
        all_vcfs_input= lambda wildcards, input: " -I ".join(input.vcf),
        all_stats_input= lambda wildcards, input: " -stats ".join(input.stats),
        all_f1r2_input = lambda wildcards, input: " -I ".join(input.f1r2),
        ref = config["ref"],
        heap_mem=4
    group: "merge"
    message: 
        "{wildcards.pid} calls are being merged and filtered"
    log:
        "logs/merge_calls/{pid}.log"
    output:
        merged_unfiltered_vcf="{pid}/merged_unfiltered.vcf",
        filtering_stats="{pid}/filtering.stats",
        merged_stats="{pid}/merged_unfiltered.stats",
        artifact_prior="{pid}/artifact_prior.tar.gz",
        merged_filtered_vcf="{pid}/merged_filtered.vcf",
    shell:
        """
        gatkm="gatk --java-options -Xmx{params.heap_mem}g"
        echo "-----------------------learn read orientation model-----------------------"
        $gatkm LearnReadOrientationModel \
            -I {params.all_f1r2_input} \
            -O {output.artifact_prior} &> {log}

        echo "----------------------- merge vcfs----------------------------"
        $gatkm MergeVcfs \
            -I {params.all_vcfs_input} \
            -O {output.merged_unfiltered_vcf} &>> {log}

        echo "------------------------merge mutect stats-------------------------------"
        $gatkm MergeMutectStats \
            -stats {params.all_stats_input} \
            -O {output.merged_stats} &>> {log}

        echo --------------filter mutect stats ----------------
        $gatkm FilterMutectCalls \
            -V {output.merged_unfiltered_vcf} \
            -R {params.ref} \
            --contamination-table {input.contamination_table} \
            --tumor-segmentation {input.segments_table} \
            --ob-priors {output.artifact_prior} \
            -stats {output.merged_stats} \
            --filtering-stats {output.filtering_stats} \
            -O {output.merged_filtered_vcf} &>> {log}
        """



rule funcotate:
    input:
        merged_filtered_vcf="{pid}/merged_filtered.vcf"
    params:
        data_source_folder=config["onco_folder"],
        ref=config["ref"],
        interval_list=config["wes_interval_list"],
        heap_mem=10
    output:
        touch("{pid}/funco_flag"),
        annot_merged_filtered_maf="{pid}/annot_merged_filtered.vcf"
    message:
        "{wildcards.pid} is being annotated ...."
    log:
        "logs/funcotate/{pid}.log"
    group: "merge"
    shell:
        """
        gatkm="gatk --java-options -Xmx{params.heap_mem}g"
        echo "-----------------------annotate-----------------------"
        $gatkm Funcotator \
            --data-sources-path {params.data_source_folder} \
            --ref-version hg19 \
            --output-file-format MAF \
            -R {params.ref} \
            -V {input.merged_filtered_vcf} \
            -O {output.annot_merged_filtered_maf} \
            -L {params.interval_list} \
            --remove-filtered-variants true &> {log}
        touch {wildcards.pid}/funco_flag
        """
