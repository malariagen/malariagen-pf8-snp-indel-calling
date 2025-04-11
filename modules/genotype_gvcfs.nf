process genotype_gvcfs {
    // This process runs GATK GenotypeGVCFs on a GenomicsDB workspace
    // (untarred from input) for a specific interval, producing a genotyped VCF
    // and its index.

    input:
        tuple path(intervals), path(genomic_db_tar)
        path(ref_fasta)
        path(ref_fasta_index)
        path(ref_dict)

    output:
        tuple path(intervals), path("genotyped.vcf.gz"), path("genotyped.vcf.gz.tbi")

    script:
    jvm_Xmx_memory  = task.memory.mega-200
    java_Xms_memory = jvm_Xmx_memory
    jvm_args        = "-Xmx${jvm_Xmx_memory}m -Xms${java_Xms_memory}m"

    genomic_db = genomic_db_tar.getBaseName()
    """
    set -exo pipefail

    tar xf ${genomic_db_tar}

    cp ${intervals} intervals.bed


    gatk \\
        --java-options \"${jvm_args}\" \\
        GenotypeGVCFs \\
        -R ${ref_fasta} \\
        -O genotyped.vcf.gz \\
        --only-output-calls-starting-in-intervals \\
        --use-new-qual-calculator \\
        -V gendb://${genomic_db} \\
        -L intervals.bed \\
        --annotation-group StandardAnnotation \\
        --annotation-group AS_StandardAnnotation

    rm -Rf ${genomic_db}
    """
}
