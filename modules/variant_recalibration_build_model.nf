process variant_recalibration_build_model {
    // This process builds a VQSR model using GATK VariantRecalibrator, based on 
    // input VCFs and parameters, and outputs a tarball containing the recalibration
    // model, tranches, R plots, and related files for later variant quality score
    // recalibration.

    tag "$mode"


    input:
        path(all_input_vcfs)
        val(basename)
        val(mode)
        path(parameters_file)

    output:
        path("${output_bundle}")

    script:
        jvm_Xmx_memory = task.memory.mega-200
        jvm_Xms_memory = jvm_Xmx_memory
        jvm_args       = "-Xmx${jvm_Xmx_memory}m -Xms${jvm_Xms_memory}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"

        output_prefix      = "${basename}_${mode}"
        output_bundle      = "${output_prefix}.recal.tar.gz"
        output_recal       = "${output_prefix}.recal"
        output_recal_index = "${output_prefix}.recal.idx"
        output_tranches    = "${output_prefix}.tranches"
        output_model       = "${output_prefix}.model"
        output_rplot       = "${output_prefix}.rplot"

        """
        set -x
        echo ${mode}
        
        cp ${all_input_vcfs} input_list.list
        cp ${parameters_file} parameters.args

        gatk --java-options "${jvm_args}" \\
            VariantRecalibrator \\
            -V input_list.list \\
            --arguments_file parameters.args \\
            -O ${output_recal} \\
            --tranches-file ${output_tranches} \\
            --trust-all-polymorphic \\
            -mode ${mode} \\
            --output-model ${output_model} \\
            --rscript-file ${output_rplot}
        
        echo "${output_prefix}" > .prefix

        FILES="${output_recal} ${output_recal_index} ${output_model} ${output_tranches} ${output_rplot} .prefix"
        which zip
        tar czf "${output_bundle}" \$FILES
        #rm \$FILES 

        """
}
