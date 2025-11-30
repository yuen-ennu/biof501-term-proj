#!/usr/bin/env nextflow
nextflow.enable.dsl=2


workflow {
    // ids = ["sra1","sra2", "sra3"]

    // Download references
    ch_human_ref = DOWNLOAD_HUMAN_DB()
    ch_viral_ref = DOWNLOAD_VIRAL_DB()

    // SRA_ID(s) of human gut metagenome smaples (WHOLE GENOME SEQUENCING)
    ids = ["ERR14295483"] 
        // gut (3) - wgs (2) and amplicon(1)
        // running into probelm w read 3

    // ids = ["ERR13983914"]  // human gut metagenome



    // Preprocessing and QC
    ch_reads = Channel.fromSRA(ids)
    FASTP(ch_reads)
    ch_fastp = FASTP.out.reads
    FASTQC(ch_fastp)

    // Find non-human-mapped reads
    ch_unmapped = MINIMAP2(ch_fastp, ch_human_ref)

    // Classification of viral reads
    ch_viral_db = MAKE_VIRAL_DB(ch_viral_ref)
    ch_kraken2 = KRAKEN2(ch_unmapped, ch_viral_db)
    ch_bracken = BRACKEN(ch_kraken2, ch_viral_db)

    // Viral taxonomy visualization 
    KRONA(ch_bracken)
}


process DOWNLOAD_HUMAN_DB {
    errorStrategy 'retry'
    maxRetries 2
    tag "Download Reference"
    output:
        path "hg38.fa.gz"
    script:
    """
    echo "DOWNLOADING HUMAN REFERENCE..."
    curl -L https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz -o hg38.fa.gz
    """
}


process FASTP {
    publishDir "${projectDir}/data/fastp"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), 
            path("fastp_${sample}_1.fastq"), 
            path("fastp_${sample}_2.fastq"), 
            emit:reads
        path("fastp_${sample}.fastp.html")
    script:
    def(r1,r2) = reads
    """
    echo "Preprocessing and filtering for paired-read files..."
    fastp \\
        -i "${r1}" \\
        ${r2 ? "-I ${r2}" : ""} \\
        -o "fastp_${sample}_1.fastq" \\
        ${r2 ? "-O fastp_${sample}_2.fastq" : ""} \\
        --html "fastp_${sample}.fastp.html"
    """
    }

process FASTQC {
    publishDir "${projectDir}/data/fastqc"
    input:
        tuple val(sample), path(r1), path(r2)
    output:
        path("*.html")
    script:
    """
    fastqc ${r1} ${r2}
    """
}

process MINIMAP2 {
    publishDir "${projectDir}/data/minimap2"
    input:
        tuple val(sample), path(r1), path(r2)
        path(reference)
    output:
        tuple val(sample), 
            path("${sample}_unmapped_1.fastq"), 
            path("${sample}_unmapped_2.fastq")
   script:
    """
    echo "Finding non-hg38 reads..."
    minimap2 -a -x sr ${reference} ${r1} ${r2} \
        | samtools view -b -f 4 - \
        | samtools fastq \
            -1 ${sample}_unmapped_1.fastq \
            -2 ${sample}_unmapped_2.fastq \
            -0 /dev/null \
            -s /dev/null
    """
}

process DOWNLOAD_VIRAL_DB{
    errorStrategy 'retry'
    maxRetries 2
    tag "Download Kraken2 Viral Reference"
    output:
        path "k2_viral_20251015.tar.gz"
    script:
    """
    echo "Downloading KRAKEN2 viral reference..."
    curl -L https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20251015.tar.gz -o k2_viral_20251015.tar.gz
    """
}
process MAKE_VIRAL_DB {
    input:
        path db_tar
    output:
        path "viral_db"
    script:
    """
    echo "Making KRAKEN2 viral tar database..."
    mkdir -p viral_db
    tar -xvzf ${db_tar} -C viral_db
    """
}

process KRAKEN2 {
    publishDir "${projectDir}/data/kraken2"
    input:
        tuple val(sample), path(r1), path(r2)
        path viral_db_folder
    output:
        tuple val(sample), path("${sample}_kraken2_report.txt")
    script:
    """
    echo "Finding viral taxonomy..."
    kraken2 --db ${viral_db_folder} \\
        --report ${sample}_kraken2_report.txt \\
        --output ${sample}_kraken2_output.txt \\
        --paired "${r1}" "${r2}"
    """
}

process BRACKEN {
    publishDir "${projectDir}/data/bracken"
    input:
        tuple val(sample), path(kraken_report)
        path viral_db_folder
    output:
        tuple val(sample), path("${sample}_bracken_report.txt")
    script:
    """
    echo "Finding viral abundance..."
    bracken -d ${viral_db_folder} \
            -i ${kraken_report} \
            -r 100 -l S -t 10 \
            -o ${sample}_b_output.txt \
            -w ${sample}_bracken_report.txt
    """
}

process KRAKEN2KRONA {
    publishDir "${projectDir}/data/kreport2krona"
    input:
        tuple val(sample), path(bracken_report)
    output:
        tuple val(sample), path("${sample}_krona.txt")
    script:
    """
    python ${projectDir}/KrakenTools/kreport2krona.py -r ${bracken_report} \
        -o ${sample}_krona.txt --no-intermediate-ranks
    """
}

process KRONA {
    publishDir "${projectDir}/data/krona"
    input:
        tuple val(sample), path(krona_result)
    output:
        path("${sample}_krona.html")
    """
    echo "Making taxonomy visual at ${sample}_krona.html..."
    ktImportText ${krona_result} -o ${sample}_krona.html
    """
}