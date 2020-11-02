rule atac_bammkdp:
    input:
        bam = "Result/minimap2/{fastqid}.sortedByPos.bam"
    output:
        bam = "Result/minimap2/{fastqid}.sortedByPos.mkdp.bam",
        metric = "Result/minimap2/{fastqid}.sortedByPos.mkdp.txt",
        tmp = temp(directory("Result/Tmp/{fastqid}"))
    shell:
        "picard MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metric} TMP_DIR={output.tmp};"
        "rm {input.bam}"

rule chip_bamrmdp:
    input:
        bam = "{OUT_DIR}/minimap2/{fastqid}.sortedByPos.mkdp.bam",
    output:
        bam = "{OUT_DIR}/minimap2/{fastqid}.sortedByPos.rmdp.clean.bam",
    params:
        chrombed = config["annotation"]["chromBed"],
    shell:
        "samtools view --threads {threads} -b -L {params.chrombed} -F 0x400 -o {output.bam} {input.bam};"
	
rule atac_callpeak:
    input:
        bam = "Result/minimap2/{fastqid}.sortedByPos.rmdp.clean.bam" 
    output:
        peak = "Result/Analysis/{fastqid}_peaks.narrowPeak",
        bdg = "Result/Analysis/{fastqid}_treat_pileup.bdg",
        xls = "Result/Analysis/{fastqid}_peaks.xls",
    params:
        name = "{fastqid}",
    log:
        "Result/Log/{fastqid}_macs2_peak.log"
    benchmark:
        "Result/Benchmark/{fastqid}_callpeak.benchmark"
    shell:
        "macs2 callpeak -g hs --outdir Result/Analysis -n {params.name} --keep-dup all -B -q 0.05 -f BAMPE --SPMR -t {input.bam};"

rule chip_bdg2bw:
    input:
        bdg = "{OUT_DIR}/Analysis/{bdgname}.bdg",
    output:
        bw = "{OUT_DIR}/Analysis/{bdgname}.bw",
    params:
        chromlen = config["annotation"]["chromInfo"],
    shell:
        "bdg2bw {input.bdg} {params.chromlen};"