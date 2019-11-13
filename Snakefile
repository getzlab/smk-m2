import pandas as pd
localrules: all


shell.prefix("set -eo pipefail; echo BEGIN at $(date);  ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")


configfile: "config.yaml"
workdir: "res/"

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
        expand("{pid}/aligned_merged_filtered.vcf",pid=PID),
        expand("{pid}/annot_merged_filtered.vcf",pid=PID)


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
            -O {output[0]}

        """



# get signed url from nci-data-commons
rule get_url:
    input:
        rules.fs_setup.output,
        rules.split_intervals.output
    params:
        pid=lambda wildcards: wildcards.pid,
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
    params:
        ref=config["ref"],
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
        cat {input.url}
	echo --------------
        wget "$(cat {input.url})" -O {output.bam}
        gatk BuildBamIndex -I {output.bam} -O {output.bai}
        
        [ gsutil -q stat gs://{params.bucket}/{output.bai} ] && gsutil cp {output.bai} gs://{params.bucket}/{wildcards.pid}_{wildcards.type}.bai
        
        gatk GetSampleName -R {params.ref} -I {output.bam} -O {wildcards.pid}/{wildcards.type}_name.txt
        cat {output.sample_name}
        """

# run M2 on tumor-normal pairs for 
rule scatter:
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
        out_stats="{pid}/subvcfs/{chr}.stats"
        
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
        
        echo finish Mutect2 ----------

        $gatkm GetPileupSummaries \
            -R {params.ref} -I $tumor \
            --interval-set-rule INTERSECTION \
            -L $interval_file \
            -V {params.vfc} \
            -L {params.vfc} \
            -O {output.out_tpile}
        echo Getpiles tumor ---------

        $gatkm GetPileupSummaries \
            -R {params.ref} -I $normal \
            --interval-set-rule INTERSECTION \
            -L $interval_file \
            -V {params.vfc} \
            -L {params.vfc} \
            -O {output.out_npile}
        echo Getpiles normal --------

        """

rule check_all_intervals:
    input:
        #rules.scatter.output
        expand("{{pid}}/tpile_dir/{chr}-tpile.table", chr=SUBINTS),
        expand("{{pid}}/npile_dir/{chr}-npile.table", chr=SUBINTS),
        expand("{{pid}}/f1r2/{chr}-f1r2.tar.gz", chr=SUBINTS),
        expand("{{pid}}/subvcfs/{chr}.vcf", chr=SUBINTS)
    output:
        touch("{pid}/flag")
    log:
        "logs/check_all_chrs/{pid}.log"
    shell:
        """
        echo "scatter finished" > {output}
        """


rule merge:
    input:
        flag="{pid}/flag",
        tumor_pile_table=expand("{{pid}}/tpile_dir/{chr}-tpile.table", chr=SUBINTS),
        normal_pile_table=expand("{{pid}}/npile_dir/{chr}-npile.table", chr=SUBINTS),
        f1r2=expand("{{pid}}/f1r2/{chr}-f1r2.tar.gz", chr=SUBINTS),
        vcf=expand("{{pid}}/subvcfs/{chr}.vcf", chr=SUBINTS),
        stats=expand("{{pid}}/subvcfs/{chr}.stats", chr=SUBINTS)

    params:
        all_tumor_piles_input= lambda wildcards, input: " -I ".join(input.tumor_pile_table),
        all_normal_piles_input=lambda wildcards, input: " -I ".join(input.normal_pile_table),
        all_f1r2_input = lambda wildcards, input: " -I ".join(input.f1r2),
        
        all_vcfs_input= lambda wildcards, input: " -I ".join(input.vcf),
        all_stats_input= lambda wildcards, input: " -stats ".join(input.stats),
        interval_list=config["wes_interval_list"],
        heap_mem=config["heapmem"],
        ref=config["ref"],
        ref_dict=config["ref_dict"],
        data_source_folder=config["onco_folder"]
    log:
        "logs/merge/{pid}.log"
    output:
        merged_unfiltered_vcf="{pid}/merged_unfiltered.vcf",
        merged_stats="{pid}/merged_unfiltered.stats",
        artifact_prior="{pid}/artifact_prior.tar.gz",
        normal_pile_table="{pid}/normal_pile.tsv",
        tumor_pile_table="{pid}/tumor_pile.tsv",
        contamination_table="{pid}/contamination.table",
        segments_table="{pid}/segments.table",
        filtering_stats="{pid}/filtering.stats",
        merged_filtered_vcf="{pid}/merged_filtered.vcf",
        aligned_merged_filtered_vcf="{pid}/aligned_merged_filtered.vcf",
        annot_merged_filtered_maf="{pid}/annot_merged_filtered.vcf"


    shell:
        """
        gatkm="gatk --java-options -Xmx{params.heap_mem}g"
        echo "----------------------- merge vcfs----------------------------"
        $gatkm MergeVcfs \
            {params.all_vcfs_input} \
            -O {output.merged_unfiltered_vcf}

        echo "------------------------merge mutect stats-------------------------------"
        $gatkm MergeMutectStats \
            -stats {params.all_stats_input} \
            -O {output.merged_stats}

        echo "-----------------------learn read orientation model-----------------------"
        $gatkm LearnReadOrientationModel \
            -I {params.all_f1r2_input} \
            -O {output.artifact_prior}

        echo "-----------------------gather pile up summaries-----------------------"
        $gatkm GatherPileupSummaries \
            -I {params.all_normal_piles_input} \
            --sequence-dictionary {params.ref_dict} \
            -O {output.normal_pile_table}

        echo "-----------------------gather tumor piles-----------------------"
        $gatkm GatherPileupSummaries \
            -I {params.all_tumor_piles_input} \
            --sequence-dictionary {params.ref_dict} \
            -O {output.tumor_pile_table}

        echo "-----------------------calc contamination-----------------------"
        $gatkm CalculateContamination \
            -I {input.tumor_pile_table} \
            -O {output.contamination_table} \
            --tumor-segmentation {output.segments_table} \
            -matched {output.normal_pile_table}

        echo "-----------------------annotate-----------------------"
        $gatkm Funcotator \
            --data-sources-path {params.data_source_folder} \
            --ref-version hg19 \
            --output-file-format MAF \
            -R {params.ref} \
            -V {output.aligned_merged_filtered_vcf} \
            -O {output.annot_merged_filtered_maf} \
            -L {params.interval_list} \
            --remove-filtered-variants true


        """
