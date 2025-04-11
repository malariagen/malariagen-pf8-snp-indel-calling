process import_gvcfs {
    // This process runs GATK GenomicsDBImport to merge multiple GVCFs
    // for a specific interval into a GenomicsDB workspace, then tars the
    // workspace and outputs it for downstream joint genotyping.

    memory "4G"

    input:
        tuple path(intervals), path(sample_name_map_file)

    output:
        tuple path(intervals), path("${genomic_db}.tar")

    script:

        jvm_Xmx_memory  = task.memory.mega-200
        java_Xms_memory = jvm_Xmx_memory
        jvm_args        = "-Xmx${jvm_Xmx_memory}m -Xms${java_Xms_memory}m"

        genomic_db       = intervals.getName() + ".gdb"
        batch_size       = params.import_gvcfs_batch_size
        interval_padding = params.import_gvcfs_interval_padding

        threads = task.cpus

        """
        cp ${intervals} intervals.bed
        set -exo pipefail

        gatk --java-options \"${jvm_args}\" \\
            GenomicsDBImport \\
            --genomicsdb-workspace-path ${genomic_db} \\
            --batch-size ${batch_size} \\
            -L intervals.bed \\
            --sample-name-map ${sample_name_map_file} \\
            --reader-threads ${threads} \\
            -ip ${interval_padding} 


        # Writes the tar into the current working directory
        tar -cf ${genomic_db}.tar ${genomic_db}
        rm -Rf ${genomic_db}
        """
}
