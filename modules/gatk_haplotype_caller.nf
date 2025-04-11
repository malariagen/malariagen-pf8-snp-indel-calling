process gatk_haplotype_caller {
    // This process runs GATK HaplotypeCaller on a BAM file over
    // a specific genomic interval to produce a GVCF (.vcf.gz),
    // with optional Pf8-tuned Smith-Waterman parameters, and outputs
    // the compressed variant file and its index.

    publishDir "${params.results_dir}/gvcfs/${sample_id}/", mode: 'copy'

    memory { 8.GB * task.attempt + 2.GB } 

    // No point of using more than one thread, parallelism is implemented thru scatter gather
    cpus 1

    input:
        tuple path(interval_file), val(sample_id), path(bam_file), path(index_file)
        path(reference_file)
        path(reference_index_file)
        path(dict_file)

    output:
        tuple path(interval_file), val(sample_id), path(vcf_file_gz), path(vcf_file_gz_index)

    script:

        jvm_memory = task.memory.mega-200
        jvm_args = "-Xmx${jvm_memory}m"

        interval_name = interval_file.getBaseName();

        base_name         = bam_file.getBaseName()
        vcf_file_gz       = base_name + "_" + interval_name + ".vcf.gz"
        vcf_file_gz_index = vcf_file_gz + ".tbi"
        
        contamination = params.gatk_haplotype_caller_contamination

        sw_opts = ""

        if (params.gatk_haplotype_caller_execution_parameters) {
            sw_opts = """--smith-waterman-dangling-end-gap-extend-penalty -6 \\
                --smith-waterman-dangling-end-gap-open-penalty -110 \\
                --smith-waterman-dangling-end-match-value 25 \\
                --smith-waterman-dangling-end-mismatch-penalty -50 \\
                --smith-waterman-haplotype-to-reference-gap-extend-penalty -6 \\
                --smith-waterman-haplotype-to-reference-gap-open-penalty -110 \\
                --smith-waterman-haplotype-to-reference-match-value 25 \\
                --smith-waterman-haplotype-to-reference-mismatch-penalty -50 \\
                --smith-waterman-read-to-haplotype-gap-extend-penalty -5 \\
                --smith-waterman-read-to-haplotype-gap-open-penalty -30 \\
                --smith-waterman-read-to-haplotype-match-value 10 \\
                --smith-waterman-read-to-haplotype-mismatch-penalty -15"""
        }


        """
        cp ${interval_file} intervals.bed

        # GATK4 Haplotyper has issues with OpenJDK 11

        # Run haplotype caller
        gatk --java-options ${jvm_args} \\
            HaplotypeCaller \\
            -R ${reference_file} \\
            -I ${bam_file} \\
            -O ${vcf_file_gz} \\
            -L intervals.bed \\
            -contamination ${contamination} \\
            --create-output-variant-index true \\
            -ERC ${params.gatk_haplotype_caller_out_mode} \\
            ${sw_opts}
        """
}
