rule chip_bammkdp:
    input:
        bam = "Result/minimap2/{fastqid}.sortedByPos.bam",
    output:
        bam = "Result/minimap2/{fastqid}.sortedByPos.mkdp.bam",
        metric = "Result/minimap2/{fastqid}.sortedByPos.mkdp.txt",
        tmp = temp(directory("Result/Tmp/{fastqid}"))
    shell:
        "picard MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metric} TMP_DIR={output.tmp};"
        "rm {input.bam}"

rule chip_callpeak:
    input:
        bam_t = "Result/minimap2/{expid}_treatment.sortedByPos.mkdp.bam",
        bam_c = "Result/minimap2/{expid}_control.sortedByPos.mkdp.bam",
    output:
        bam_t = "Result/minimap2/{expid}_treatment.sortedByPos.rmdp.clean.bam",
        bam_c = "Result/minimap2/{expid}_control.sortedByPos.rmdp.clean.bam",	
        peak = "Result/Analysis/{expid}_peaks.narrowPeak",
        xls = "Result/Analysis/{expid}_peaks.xls",
        bdg_t = "Result/Analysis/{expid}_treat_pileup.bdg",
        bw_t = "Result/Analysis/{expid}_treat_pileup.bw",
        bdg_c = "Result/Analysis/{expid}_control_lambda.bdg",
        bw_c = "Result/Analysis/{expid}_control_lambda.bw",
    params:
        name = "{expid}",
        chrombed = config["annotation"]["chromBed"],
     	chromlen = config["annotation"]["chromInfo"],
	macs2_option = "{macs2_genome_option} {macs2_cutoff_option} {macs2_lambda_option} {macs2_format_option}"
    log:
        "Result/Log/{expid}_macs2_peak.log"
    benchmark:
        "Result/Benchmark/{expid}_callpeak.benchmark"
    shell:
        "samtools view --threads {threads} -b -L {params.chrombed} -F 400 -o {output.bam_t} {input.bam_t}; "
        "samtools view --threads {threads} -b -L {params.chrombed} -F 400 -o {output.bam_c} {input.bam_c}; "	
        "macs2 callpeak --outdir Result/Analysis -n {params.name} --keep-dup all -B --SPMR -t {output.bam_t} -c {output.bam_c}; "
        "bdg2bw {output.bdg_t} {params.chromlen}; "
        "bdg2bw {output.bdg_c} {params.chromlen}; "
