process mark_duplicates {
    // This process marks duplicate reads in a BAM file using GATK's MarkDuplicates. 
    // It generates a new BAM file with duplicates marked, along with its index. 
    // The process sets memory usage, defines temporary directories, and adjusts parameters for
    // handling large datasets. If the index is not created by GATK, it falls back to using
    // samtools for indexing. It also includes validation checks using samtools to ensure
    // the integrity of the resulting BAM file.

    tag "$base_name"

    input:
        tuple val(sample_id), path(bam_file), path(bam_file_index)

    output:
        tuple val(sample_id), path("$markdup_bam_file"), path("$markdup_bam_index")

    script:
        base_name             = bam_file.getBaseName()
        markdup_bam_file      = base_name + ".markdup.bam"
        markdup_bam_index     = base_name + ".markdup.bam.bai"
        markdup_alt_bam_index = base_name + ".markdup.bai"
        jvm_memory            = task.memory.mega-200
        jvm_args              = "-Xmx${jvm_memory}m"

        
        // got an informal rule of thumb estimate of 250_000 for each 1G for reads 100bp long.
        // assuming that reads are usually at most 150bp long then...: 
        max_records_in_ram = jvm_memory.intdiv(1500) * 250000
        
        """
        set -euxo pipefail
        mkdir -p tmp
        gatk --java-options "${jvm_args}" MarkDuplicates \\
            -I ${bam_file} \\
            -O ${markdup_bam_file} \\
            --ASSUME_SORT_ORDER coordinate \\
            --METRICS_FILE /dev/null \\
            --VALIDATION_STRINGENCY SILENT \\
            --MAX_RECORDS_IN_RAM ${max_records_in_ram} \\
            --TMP_DIR tmp \\
            --CREATE_INDEX true 
        rm -Rf tmp
        ### changing from X.bai to X.bam.bai for the index to keep it consistent with other tools.
        if [[ -f "${markdup_alt_bam_index}" ]]; then
            mv "${markdup_alt_bam_index}" "${markdup_bam_index}"
        ### fall-back to relay on samtools if the index files wasn't produced
        elif [[ ! -f "${markdup_bam_index}" ]]; then
            samtools index ${markdup_bam_file}
        fi
        ### 
        samtools quickcheck ${markdup_bam_file}
        samtools view -bS ${markdup_bam_file} > /dev/null
        """
}