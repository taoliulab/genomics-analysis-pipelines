rule chip_bammkdp:
    input:
        bam = "{OUT_DIR}/minimap2/{fastqid}.sortedByPos.bam",
    output:
        bam = "{OUT_DIR}/minimap2/{fastqid}.sortedByPos.mkdp.bam",
        metric = "{OUT_DIR}/minimap2/{fastqid}.sortedByPos.mkdp.txt",
        tmp = temp(directory("{OUT_DIR}/Tmp/{fastqid}"))
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

rule chip_callpeak:
    input:
        bam_t = "{OUT_DIR}/minimap2/{expid}_treatment.sortedByPos.rmdp.clean.bam",
        bam_c = "{OUT_DIR}/minimap2/{expid}_control.sortedByPos.rmdp.clean.bam",
    output:
        peak = "{OUT_DIR}/Analysis/{expid}_peaks.narrowPeak",
        xls = "{OUT_DIR}/Analysis/{expid}_peaks.xls",
        bdg_t = "{OUT_DIR}/Analysis/{expid}_treat_pileup.bdg",
        bdg_c = "{OUT_DIR}/Analysis/{expid}_control_lambda.bdg",
    params:
        name = "{expid}",
    log:
        "{OUT_DIR}/Log/{expid}_macs2_peak.log"
    benchmark:
        "{OUT_DIR}/Benchmark/{expid}_callpeak.benchmark"
    shell:
        "macs2 callpeak --outdir {OUT_DIR}/Analysis -n {params.name} --nomodel --extsize 150 --keep-dup all -B --SPMR -t {input.bam_t} -c {input.bam_c} {macs2_option}; "

rule chip_bdg2bw:
    input:
        bdg = "{OUT_DIR}/Analysis/{bdgname}.bdg",
    output:
        bw = "{OUT_DIR}/Analysis/{bdgname}.bw",
    params:
        chromlen = config["annotation"]["chromInfo"],
    shell:
        "bdg2bw {input.bdg} {params.chromlen};"