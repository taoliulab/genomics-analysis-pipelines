rule atac_qcstat:
    input:
        bam = "Result/minimap2/{fastqid}.sortedByPos.rmdp.clean.bam",
        peak = "Result/Analysis/{fastqid}_peaks.narrowPeak",
    output:
        qc_stat = "Result/QC/{fastqid}.stat.txt",
        bam = "Result/minimap2/{fastqid}.sortedByPos.rmdp.clean.unique.bam",
        bed = "Result/minimap2/{fastqid}.sortedByPos.rmdp.clean.unique.bed",
    params:
        promoter = config["annotation"]["promoter"]
    threads:
        config["options"]["cores"]
    benchmark:
        "Result/Benchmark/{fastqid}_BulkQCStat.benchmark"
    shell:
        "echo 'flagstat:' > {output.qc_stat};"
        "samtools flagstat --threads {threads} {input.bam} >> {output.qc_stat};"
        "samtools view -F 2316 -f 0x2 -q 30 -b -o {output.bam} {input.bam};"
        "echo 'mapped Q30 reads:' >> {output.qc_stat};"
        "samtools view {output.bam} -c >> {output.qc_stat};"
        "bedtools bamtobed -i {output.bam} > {output.bed};"
        "echo 'chrM reads:' >> {output.qc_stat};"
        "grep -c 'chrM' {output.bed} >> {output.qc_stat} || true;"
        "echo 'non chrM reads:' >> {output.qc_stat};"
        "grep -v 'chrM' -c {output.bed} >> {output.qc_stat};"
        "echo 'non chrM reads in promoter:' >> {output.qc_stat};"
        "grep -v 'chrM' {output.bed} | bedtools intersect -wa -a - -b {params.promoter} -u | wc -l >> {output.qc_stat} || true;"
        "echo 'non chrM reads in peak:' >> {output.qc_stat};"
        "grep -v 'chrM' {output.bed} | bedtools intersect -wa -a - -b {input.peak} -u | wc -l >> {output.qc_stat} || true ;"

rule atac_peakqc:
    input:
        peak = "Result/Analysis/{fastqid}_peaks.narrowPeak",
        peakxls = "Result/Analysis/{fastqid}_peaks.xls",
    output:
        peak_qc = "Result/QC/{fastqid}.peakstat.txt",
    params:
        promoter = config["annotation"]["promoter"],
        chrMregion = config["annotation"]["MtBed"],
        blacklist = config["annotation"]["blacklist"],
        DHS = config["annotation"]["DHS"],
    threads:
        config["options"]["cores"],
    benchmark:
        "Result/Benchmark/{fastqid}_PeakQCStat.benchmark",
    shell:
        "grep 'total fragments in treatment' {input.peakxls} | perl -pe 's/# //' > {output.peak_qc};"
        #"grep 'fragments after filtering in treatment' {input.peakxls} | perl -pe 's/# //' >> {output.peak_qc};"
        "echo 'total number of peaks:' >> {output.peak_qc};"
        "wc -l {input.peak} | cut -f 1 -d' ' >> {output.peak_qc};"
        "echo 'number of peaks over FC 2:' >> {output.peak_qc};"
        "awk '$7>=2{{print}}' {input.peak} | wc -l >> {output.peak_qc};"
        "echo 'number of peaks in blacklist regions:' >> {output.peak_qc};"
        "bedtools intersect -a {input.peak} -b {params.blacklist} -u | wc -l >> {output.peak_qc};"
        "echo 'number of peaks in chrM:' >> {output.peak_qc};"
        "grep -w chrM {input.peak} | wc -l >> {output.peak_qc};"
        "echo 'number of peaks in promoter regions:' >> {output.peak_qc};"
        "bedtools intersect -a {input.peak} -b {params.promoter} -u | wc -l >> {output.peak_qc};"
        "echo 'number of peaks in DHS regions:' >> {output.peak_qc};"
        "bedtools intersect -a {input.peak} -b {params.DHS} -u | wc -l >> {output.peak_qc};"

rule atac_frag:
    input:
        bam = "Result/minimap2/{fastqid}.sortedByPos.rmdp.clean.unique.bam",
    output:
        insertl = "Result/minimap2/{fastqid}.sortedByPos.rmdp.clean.unique.bam.insertl",
        insertlsum = "Result/minimap2/{fastqid}.sortedByPos.rmdp.clean.unique.bam.insertl.txt",        
    shell:
        "samtools view {input.bam} | cut -f 9 | awk '$1>0{{print}}' > {output.insertl};"
        "perl -e 'while(<>){{chomp;$bin=int($_/10);$count{{$bin}}+=1}}foreach my $key (sort {{$a<=>$b}} keys %count){print $key*10,'\t',$count{{$key}},'\n'}' {output.insertl} > {output.insertlsum};"
