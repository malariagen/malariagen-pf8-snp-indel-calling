process samtools_stats {
    // This process generates statistics for a BAM file using samtools stats. 
    // It runs the samtools stats command on the input BAM file, using the available 
    // CPU cores to speed up the process. The resulting statistics are saved to a file 
    // named ${base_name}.bamstats. The output file is published to a specified 
    // directory. The process allows for optional configuration of samtools options 
    // (though in this case, the options variable is empty).

    cpus 4

    publishDir "./${params.results_dir}/bams/samtools-stats-files/", mode:'copy'

    input:
        path(bam_file)

    output:
        path("${stats_file}")

    script:
        options               = "" // ???
        base_name             = bam_file.getBaseName()
        stats_file            = "${base_name}.bamstats"	
        samtools_threads_args = task.cpus <= 1 ? "" : "-@ " + (task.cpus - 1)
        
        """
        samtools stats ${samtools_threads_args} ${options} ${bam_file} > ${stats_file}
        """
}
