process variant_filtration {
    // This process applies variant filtration to a VCF file using GATK's
    // VariantFiltration tool. It filters variants based on several criteria,
    // including VQSLOD value, region type (e.g., centromere, subtelomeric,
    // mitochondrion), and the presence of VQSLOD. It generates a filtered VCF
    // and its index, with multiple filtering steps, including invalidating
    // previous filters before applying new ones.
    
    memory { "" + (2 + 8 * task.attempt) + "G" }

    publishDir "${params.results_dir}/filtered_vcfs/", mode: 'copy'   

    input:
        tuple path(in_vcf), path(in_vcf_idx)
        each path(reference_fasta)
        each path(reference_fasta_fai)
        each path(reference_fasta_dict)

    output:
        tuple path(out_vcf), path(out_vcf_idx)

    script:
        jvm_Xmx_memory = task.memory.mega-200
        jvm_Xms_memory = jvm_Xmx_memory
        jvm_args       = "-Xmx${jvm_Xmx_memory}m -Xms${jvm_Xms_memory}m"


        vqslod_threshold = params.variant_filtration_vqslod_threshold
        base_name        = in_vcf.getName().replaceAll(/(\.ann)?\.vcf(\.gz)?$/,"")
        out_vcf          = base_name + ".filt.vcf.gz"
        out_vcf_idx      = out_vcf + ".tbi"

        """
        gatk --java-options "${jvm_args}" \\
            VariantFiltration \\
            -V ${in_vcf} \\
            -O filter_cleared.vcf.gz \\
            --invalidate-previous-filters || exit 101

        gatk --java-options "${jvm_args}" \\
            VariantFiltration \\
            -V filter_cleared.vcf.gz \\
            -O ${out_vcf} \\
            --reference ${reference_fasta} \\
            --filter-name "MissingVQSLOD" --filter-expression "!vc.hasAttribute('VQSLOD')" \\
            --filter-name "Low_VQSLOD" --filter-expression "VQSLOD <= ${vqslod_threshold}" \\
            --filter-name "Centromere" --filter-expression "RegionType == 'Centromere'" \\
            --filter-name "InternalHypervariable" --filter-expression "RegionType == 'InternalHypervariable'" \\
            --filter-name "SubtelomericHypervariable" --filter-expression "RegionType == 'SubtelomericHypervariable'" \\
            --filter-name "SubtelomericRepeat" --filter-expression "RegionType == 'SubtelomericRepeat'" \\
            --filter-name "Apicoplast" --filter-expression "RegionType == 'Apicoplast'" \\
            --filter-name "Mitochondrion" --filter-expression "RegionType == 'Mitochondrion'" || exit 102
        """
}
