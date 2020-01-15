import pandas as pd


####################### shell setup #################

shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

####################### external configs ############
configfile: "config.yaml"
workdir: "res/"
# for a local test, will be deleted for cloud mode
# localrules: all


####################### global parameters ###########
# input is paired uuid table from GDC api.
df = pd.read_csv(config["samples"], sep = "\t")
df = df.dropna() # remove all possible NAs

TYPES=["tumor", "normal"]
melted_df = pd.melt(df, id_vars=['pid'], value_vars=["tumor","normal"])


# the subintervals to scatter on
SUBINTS = ['{:0>4}'.format(i) for i in range(config["pieces"])]
INTERVALS = ["{}/{}-scattered.interval_list".format(config["subint_dir"],int1) for int1 in SUBINTS]


# all the target files
rule all:
    input: 
        expand("{pid}/funco_flag", pid=df["pid"]),
        annot_merged_filtered_maf=expand("{pid}/annot_merged_filtered.vcf", pid=df["pid"]),
        merged_filtered_vcf=expand("{pid}/merged_filtered.vcf", pid=df["pid"])

######################## preparation #######################

# split intervals for scattering jobs
rule split_intervals:
    priority: 3000
    params:
        ref=config["ref"],
        tumor = df.tumor.values[0],
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
        ../auth.sh {params.tumor}
        """

rule retrieve_gs_path:
    log: "logs/retrieve/{pid}"
    priority: 5000 # vm starts with user auth
    params:
        tumor=lambda wildcards: df.loc[df["pid"] == wildcards.pid, "tumor"].iloc[0],
        normal=lambda wildcards: df.loc[df["pid"] == wildcards.pid,"normal"].iloc[0]
    output:
        expand("{{pid}}/{type}_url.txt", type = TYPES)
    script:
        "drs-translator.py"

rule localize_gs:
    log: "logs/localize_gs/{pid}"
    input: expand("{{pid}}/{type}_url.txt", type = TYPES)
    params: 
        tumor=lambda wildcards: df.loc[df["pid"] == wildcards.pid, "tumor"].iloc[0]
    output:
        tumor_bam = temp("{pid}/tumor.bam"),
        normal_bam = temp("{pid}/normal.bam"),
        tumor_bai = temp("{pid}/tumor.bai"),
        normal_bai = temp("{pid}/normal.bai"),
        check = "{pid}/loc.fin"
    shell:
        """
        sleep $[ ( $RANDOM % 1000 )  + 1 ]s
        fallocate -l 20G BAM_HOLDER
        ../auth.sh {params.tumor}
        gcloud auth activate-service-account --key-file=$HOME/auth.json
        [ -f {output.normal_bam} ] || gsutil cp `cat {wildcards.pid}/normal_url.txt` {output.normal_bam}
        [ -f {output.normal_bai} ] || gsutil cp `cat {wildcards.pid}/normal_url.txt`.bai {output.normal_bai}
        rm BAM_HOLDER
        # 
        [ -f {output.tumor_bam} ] || gsutil cp `cat {wildcards.pid}/tumor_url.txt` {output.tumor_bam}
        [ -f {output.tumor_bai} ] || gsutil cp `cat {wildcards.pid}/tumor_url.txt`.bai {output.tumor_bai}
        touch {wildcards.pid}/loc.fin
        """

rule get_file_name:
    input:
        check = "{pid}/loc.fin",
        tumor_bam = "{pid}/tumor.bam",
        normal_bam = "{pid}/normal.bam",
        tumor_bai = "{pid}/tumor.bai",
        normal_bai = "{pid}/normal.bai"
    log: "logs/get_name/{pid}.log"
    output:
        tumor_name = "{pid}/tumor_name.txt",
        normal_name = "{pid}/normal_name.txt"
    shell:
        """
        gatk GetSampleName -I {input.normal_bam} -O {output.normal_name} 
        gatk GetSampleName -I {input.tumor_bam} -O {output.tumor_name} 
        """



# run M2 on tumor-normal pairs for 
rule scatter_m2:
    input:
        subint_dir = config["subint_dir"],
        name_files = expand("{{pid}}/{type}_name.txt", type=TYPES),
        tumor_bam = "{pid}/tumor.bam",
        normal_bam = "{pid}/normal.bam",
        tumor_bai = "{pid}/tumor.bai",
        normal_bai = "{pid}/normal.bai"
    priority: 100
    params:
        subint_dir=config["subint_dir"],
        ref=config["ref"],
        vfc=config["vfc"],
        gnomad=config["gnomad"],
        default_pon=config["default_pon"],
        user_proj=config["user_project"]
    log:
        "logs/scatter/{pid}_{chr}.log"
    output:
        out_tpile=temp("{pid}/tpile_dir/{chr}-tpile.table"),
        out_npile=temp("{pid}/npile_dir/{chr}-npile.table"),
        out_f1r2=temp("{pid}/f1r2/{chr}-f1r2.tar.gz"),
        out_vcf=temp("{pid}/subvcfs/{chr}.vcf"),
        out_stats=temp("{pid}/subvcfs/{chr}.vcf.stats")
    shell:
        """
        gatkm="gatk --java-options -Xmx4g"
        tumor_name="{wildcards.pid}/tumor_name.txt"
        normal_name="{wildcards.pid}/normal_name.txt"

        interval_file="{params.subint_dir}/{wildcards.chr}-scattered.interval_list"

        normal_command_line="-I {input.normal_bam} -normal `cat ${{normal_name}}`"
        tumor_command_line="-I {input.tumor_bam} -tumor `cat ${{tumor_name}}`"

        $gatkm Mutect2 \
            -R {params.ref} \
            $tumor_command_line \
            $normal_command_line \
            --germline-resource {params.gnomad} \
            -pon {params.default_pon} \
            -L  $interval_file \
            -O {output.out_vcf} \
            --f1r2-tar-gz {output.out_f1r2} \
            --gcs-project-for-requester-pays {params.user_proj}\
            --independent-mates
        
        echo finish Mutect2 --------------------

        $gatkm GetPileupSummaries \
            -R {params.ref} -I {input.tumor_bam} \
            --interval-set-rule INTERSECTION \
            -L $interval_file \
            -V {params.vfc} \
            -L {params.vfc} \
            -O {output.out_tpile} 

        echo Getpiles tumor --------------------

        $gatkm GetPileupSummaries \
            -R {params.ref} -I {input.normal_bam} \
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

############### merge interval results ###############
rule calculate_contamination:
    input:
        flag="{pid}/flag",
        tumor_pile_table=expand("{{pid}}/tpile_dir/{chr}-tpile.table", chr=SUBINTS),
        normal_pile_table=expand("{{pid}}/npile_dir/{chr}-npile.table", chr=SUBINTS)
    priority: 1000
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
        normal_pile_table=temp("{pid}/normal_pile.tsv"),
        tumor_pile_table=temp("{pid}/tumor_pile.tsv"),
        contamination_table=temp("{pid}/contamination.table"),
        segments_table=temp("{pid}/segments.table")
    shell:
        """
        gatkm="gatk --java-options -Xmx{params.heap_mem}g"
        echo "-----------------------gather pile up summaries-----------------------"
        $gatkm GatherPileupSummaries \
            -I {params.all_normal_piles_input} \
            --sequence-dictionary "/demo-mount/M2_refs/Homo_sapiens_assembly19.dict" \
            -O {output.normal_pile_table} &> {log}

        echo "-----------------------gather tumor piles-----------------------"
        $gatkm GatherPileupSummaries \
            -I {params.all_tumor_piles_input} \
            --sequence-dictionary "/demo-mount/M2_refs/Homo_sapiens_assembly19.dict" \
            -O {output.tumor_pile_table} &>> {log}

        echo "-----------------------calc contamination-----------------------"
        $gatkm CalculateContamination \
            -I {output.tumor_pile_table} \
            -O {output.contamination_table} \
            --tumor-segmentation {output.segments_table} \
            -matched {output.normal_pile_table} &>> {log}
        
	"""



rule merge_m2:
    priority: 1000
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
        heap_mem=6
    output:
        touch("{pid}/funco_flag"),
        annot_merged_filtered_maf="{pid}/annot_merged_filtered.vcf"
    message:
        "{wildcards.pid} is being annotated ...."
    log:
        "logs/funcotate/{pid}.log"
    shell:
        """
        gatkm="gatk --java-options -Xmx{params.heap_mem}g"
        echo "-----------------------annotate-----------------------" &> {log}
        $gatkm Funcotator \
            --data-sources-path {params.data_source_folder} \
            --ref-version hg19 \
            --output-file-format MAF \
            -R {params.ref} \
            -V {input.merged_filtered_vcf} \
            -O {output.annot_merged_filtered_maf} \
            -L {params.interval_list} \
            --remove-filtered-variants true \
            --disable-sequence-dictionary-validation &> {log}
        """
