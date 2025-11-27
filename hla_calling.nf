#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {
    // ids = ["ERR2683153","DRX579471", "DRX579468"]
    ids = ["ERS509592"]
    ch_reads = Channel.fromSRA(ids)
    FASTP(ch_reads)
    ch_fastp = FASTP.out.reads
    FASTQC(ch_fastp)
    // OPTITYPE(ch_fastp)
}
// fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz

process FASTP {
    publishDir "${projectDir}/data/fastp"
    cpus 6
    input:
        tuple val(sample), path(reads)
    output:
        // tuple val(sample), path("fastp_${sample}_{1,2}.fastq"), emit: reads
        tuple val(sample), 
            path("fastp_${sample}_1.fastq"), 
            path("fastp_${sample}_2.fastq"), 
            emit:reads
        path("fastp_${sample}.fastp.html")
    script:
    def (r1, r2) = reads
        """
        fastp \\
            -i "${r1}" \\
            -I "${r2}" \\
            -o "fastp_${sample}_1.fastq" \\
            -O "fastp_${sample}_2.fastq" \\
            --html "fastp_${sample}.fastp.html"
            --thread ${task.cpus}
        """
    }

// process FASTQC {
//     publishDir "${projectDir}/data/fastqc"
//     input:
//         tuple val(sample), path(reads)
//     output:
//         path("*.html")
//     script:
//     """
//     fastqc ${reads}
//     """
// }

process FASTQC {
    publishDir "${projectDir}/data/fastqc"
    cpus 4
    input:
        tuple val(sample), path(r1), path(r2)
    output:
        path("${sample}_R1_fastqc.html")
        path("${sample}_R2_fastqc.html")
    script:
        """
        fastqc ${r1} --outdir . --threads ${task.cpus}
        fastqc ${r2} --outdir . --threads ${task.cpus}
        """
}

