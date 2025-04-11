process variant_recalibration_apply_model {
    // This process applies VQSR models to a VCF using GATK ApplyVQSR for SNP and 
    // INDEL recalibration. It first extracts SNP and INDEL recalibration bundles, 
    // applies recalibration in two steps (INDELs first, then SNPs), and outputs 
    // the recalibrated VCF and its index.
    
    memory "18G"

    input:
        tuple path(intervals), path(input_vcf), path(input_vcf_index)
        each path(snp_recal_bundle)
        each path(indel_recal_bundle)
        val(snp_filter_level)
        val(indel_filter_level)

    output:
        tuple path(intervals), path(output_vcf), path(output_vcf_index)

    script:
        jvm_Xmx_memory = task.memory.mega-200
        jvm_Xms_memory = jvm_Xmx_memory
        jvm_args       = "-Xmx${jvm_Xmx_memory}m -Xms${jvm_Xms_memory}m"

        output_vcf       = input_vcf.getName().replaceAll(/(\.bcftools)?\.vcf(\.gz)?$/,"") + ".recal.vcf.gz"
        output_vcf_index = output_vcf + ".tbi"

        """
        mkdir -p snp_recal indel_recal
        tar xzf ${snp_recal_bundle} -C snp_recal
        tar xzf ${indel_recal_bundle} -C indel_recal

        gatk --java-options "${jvm_args}" \\
                ApplyVQSR \\
                -O tmp.indel.recalibrated.vcf.gz \\
                -V ${input_vcf} \\
                --recal-file indel_recal/*.recal \\
                --tranches-file indel_recal/*.tranches \\
                --truth-sensitivity-filter-level ${indel_filter_level} \\
                --create-output-variant-index true \\
                -mode INDEL 

        gatk --java-options "${jvm_args}" \\
                ApplyVQSR \\
                -O ${output_vcf} \\
                -V tmp.indel.recalibrated.vcf.gz \\
                --recal-file snp_recal/*.recal \\
                --tranches-file snp_recal/*.tranches \\
                --truth-sensitivity-filter-level ${snp_filter_level} \\
                --create-output-variant-index true \\
                -mode SNP 

        rm -Rf snp_recal indel_recal tmp.*
        """
}
