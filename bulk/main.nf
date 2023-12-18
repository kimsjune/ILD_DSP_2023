nextflow.enable.dsl=2

params.fasta = "$projectDir/sym/hg38.fa"
params.gtf = "$projectDir/sym/gencode.v21.annotation.gtf"
params.sample = "$projectDir/reads/*_{R1,R2}.fastq.gz"
params.readLength= 100

fasta_ch = channel.fromFilePairs(params.fasta, size: 1)
gtf_ch = channel.fromPath(params.gtf)
sample_ch = channel.fromFilePairs(params.sample, size: 2)

process FASTQC {
    publishDir "$projectDir/fastqc", mode:"copy"
    cpus 8
    memory 24.GB
    time '1h'
    clusterOptions '--account=def-muram'

    input:
        tuple val(sample_id), path(fastq)
    output:
        path "fastqc_out/*", emit: dir // the * at the end is crucial 
        tuple val(sample_id), path("fastqc_out/*_fastqc.html"), emit: html
        tuple val(sample_id), path("fastqc_out/*_fastqc.zip"), emit: zip
    script:
    """
    mkdir fastqc_out
    fastqc -f fastq -t $task.cpus $fastq -o fastqc_out
    """
}
process STAR_INDEX {
    publishDir "$projectDir/index"
    cpus 8
    memory 64.GB
    time '3h'
    clusterOptions '--account=def-muram'

    input:
        tuple val(build), path(fasta)
            // instead of path fasta from channel.path()
        path gtf
        val read_length
    output: // The script will generate a folder called STARindex and files within. Nextflow will check if such directory was indeed created
        path "$build", emit: index
        // path "version.yml"
    script:
    """
    STAR --runThreadN $task.cpus --runMode genomeGenerate --genomeDir $build --genomeFastaFiles $fasta --sjdbGTFfile $gtf --sjdbOverhang $read_length 
    """
/*         cat << EOF >> version.yml
    $task.process:
        STAR: \$(STAR --version 2>&1)
    EOF */
    // STAR expects the output sub-directory to be there
}

process STAR_ALIGN {
    publishDir "$projectDir/aligned", mode: "copy"
    cpus 8
    memory 64.GB
    time '3h'
    clusterOptions '--account=def-muram'


    input:
        path index //  Exit code 137 caused by lack of memory if there are too many samples...
        tuple val(sample_id), path(sample_path)
        
    
    output: // All of the outputs from STAR must be listed here or else they won't be kept
            // Nextflow checks to see if these output files were indeed generated 
            // A minor problem is that STAR takes in output prefix, and automatically appends the rest 
            // This means that the code must be run separately to know what outputs (and their names) will be generated. 
            
        tuple val(sample_id), path("*Aligned.out.sam"), emit : bam
        //tuple val(sample_id), path("*Log.final.out"), emit: log
        //tuple val(sample_id), path("*Log.out")
        //tuple val(sample_id), path("*Log.progress.out")
        //tuple val(sample_id), path("*SJ.out.tab")
        path "log/*", emit:logDir

        // not really important to check if other files were generated. Makes dag look simpler.





    script:
    """
    STAR --runMode alignReads --runThreadN $task.cpus --genomeDir $index --readFilesCommand zcat --readFilesIn ${sample_path[0]} ${sample_path[1]} --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix $sample_id
    mkdir log 
    mv *Log.final.out log/
    """

}


process SAMTOOLS_VIEW_SORT {
    publishDir "$projectDir/alignedRemoveDup", mode: "copy"
    cpus 24
    memory 96.GB
    time '2h 30m'
    clusterOptions '--account=def-muram'

    input:
        tuple val(sample_id), path(sample_path)
    output: 
        tuple val(sample_id), path("*.sorted.bam"), emit : sortedBam

    script:
    """
    samtools view -F 1024 -b -@ $task.cpus $sample_path | samtools sort -@ $task.cpus - -o ${sample_id}.sorted.bam 
    """
    // for this command, -o must indicate the exact output file name, matching output name in output:
}

process HTSEQ_COUNT {
    publishDir "$projectDir/countTableRemoveDup", mode: "copy"
    cpus 1
    memory 64.GB
    time '5h'
    clusterOptions '--account=def-muram'

    input:
        tuple val(sample_id), path(sample_path)
        path gtf
    output:
        tuple val(sample_id), path("*.csv")
    script:
    // -i gene_name not yet tested
    """
    htseq-count -i gene_name -f bam -s no -q $sample_path $gtf > ${sample_id}.csv 
    """

}

process MULTIQC {
    publishDir "$projectDir/multiqc", mode: "copy"
    memory 64.GB
    time '1h'
    clusterOptions '--account=def-muram'

    input:
        path('*')

    output:
        path "multiqc_report.html"
    script:
    // the default output file name is known to be "multiqc_report.html" so there's no need for >
    """
    multiqc . 
    """
}

workflow {
    FASTQC(sample_ch)
    STAR_INDEX(fasta_ch, gtf_ch, params.readLength)
    STAR_ALIGN(STAR_INDEX.out.index.collect(), sample_ch)
        // .collect() is the KEY. Or else it returns a queue channel which gets used up after one task 
    SAMTOOLS_VIEW_SORT(STAR_ALIGN.out.bam)
        // assigning each workflow item to another channel doesn't work i.e.
        // ch_STARgenomeGenerate = STARgenomeGenerate(...)
        // STARalign(ch_STARgenomeGenerate.out.collect(), sample_ch) fails
        // I guess it's too redundant anyway?
    HTSEQ_COUNT(SAMTOOLS_VIEW_SORT.out.sortedBam, gtf_ch.collect()) //.collect seems to be the magic operator to iterate over multiple files when there are many vs. one inputs
    MULTIQC(FASTQC.out.dir.mix(STAR_ALIGN.out.logDir).collect())

}
