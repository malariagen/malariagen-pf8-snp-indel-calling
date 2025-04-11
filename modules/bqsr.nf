process bqsr {
    // This process performs Base Quality Score Recalibration (BQSR) using GATK. It 
    // first generates a recalibration model using the BaseRecalibrator with known sites
    // and BAM files, then applies this model to the BAM file with ApplyBQSR. If the
    // recalibration fails due to a bug, a workaround is applied by copying the original BAM 
    // file. The process outputs a recalibrated BAM file, its index, a recalibration report, 
    // and some validation checks with samtools.

    memory '4G'

    input:
        path(reference_file)
        path(reference_index_file)
        path(reference_dict)
        tuple path(known_sites), path(known_sites_index)
        tuple val(sample_id), path(bam_file), path(bam_file_index)

    output:
        tuple val(sample_id), path("$recalibrated_bam_file"), path("$recalibrated_bam_file_index"), emit: recalibrated
        path("$gatk_recalibration_report"), emit: recal_table

    script:
        base_name = bam_file.getBaseName()
        gatk_recalibration_report   = "${base_name}.bqsr_table"
        recalibrated_bam_file       = "${base_name}.recalibrated.bam"
        recalibrated_bam_file_index = "${base_name}.recalibrated.bai"


        jvm_memory = task.memory.mega-200
        jvm_args = "-Xmx${jvm_memory}m"
        
        """
        set -euxo pipefail
        # Recal model/table inference:
        gatk --java-options "${jvm_args}" BaseRecalibrator \\
            -R ${reference_file} \\
            -I ${bam_file} \\
            --known-sites ${known_sites} \\
            -O ${gatk_recalibration_report} || exit 101

        # Apply the model
        gatk --java-options "${jvm_args}" ApplyBQSR \\
            -R ${reference_file} \\
            -I ${bam_file} \\
            -O ${recalibrated_bam_file} \\
            --create-output-bam-index \\
            --bqsr-recal-file ${gatk_recalibration_report} 2> applyBQSR.err && exit 0


        # Ok apparently for some "bad" samples a current bug makes ApplyBQSR to fail with
        # an MPE... it seems that this is due to lack of good reads. 
        # here we added it for now.
        if grep -Fxq "NullPointerException" applyBQSR.err
        then
            echo "${sample_id} > nasty-fix-applied.txt
            cp ${bam_file} ${recalibrated_bam_file}
            cp ${bam_file_index} ${recalibrated_bam_file_index}
        else
            exit 102
        fi
            ## paranoia checks
            samtools quickcheck ${recalibrated_bam_file}
            samtools view -bS ${recalibrated_bam_file} > /dev/null
            """
}