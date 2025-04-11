process remap_bam_file {
    // This process remaps a BAM file using bwa and samtools. It first extracts the 
    // necessary reference index from a tarball and prepares read group (RG) information. 
    // The input BAM file is then sorted, converted to FASTQ, aligned with bwa mem, and 
    // processed through several samtools steps to fix mates, sort, and recalibrate the 
    // reads. The resulting remapped BAM file is indexed and validated. Temporary files 
    // and directories are cleaned up at the end. It also handles multi-threading for 
    // the operations based on available CPUs.
    
    cpus 8

    tag "$sample_id"

    input:
        path(reference_fasta)
        path(bwa_index_tarball)
        tuple val(sample_id), path(bam_file)

    output:
        tuple val(sample_id), path("$remapped_bam_file"), path("$index_file")

    script:
        base_name        = bam_file.getBaseName()
        remapped_bam_file= "${base_name}.remapped.bam"
        index_file       = "${remapped_bam_file}.bai"
        
        samtools_threads_args = task.cpus <= 1.1 ? "" : " -@ " + (task.cpus - 1)
        bwa_threads_args      = task.cpus <= 1.2 ? "" : " -t " + task.cpus

        """
        set -eux
        mkdir -p bwa_index tmp1 tmp2 ref_cache ref_path
        tar xf ${bwa_index_tarball} -C bwa_index
        if [[ ! -f bwa_index/prefix.txt ]]; then
            echo "Missing prefix.txt file in bwa index tarball; did you use \\"process_reference\\" to generate this index";
            exit 104
        fi
        index_prefix=bwa_index/`cat bwa_index/prefix.txt`
        
        # Get the RG info from the input:
        # Take the oportunity to fix some of the annotations like the sample id.:
        samtools view -H ${bam_file} | grep -a ^@RG | \\
            sed 's/\\tSM:[^\\t]*/\\tSM:${sample_id}/' | \\
            sed 's/\\tPG:[^\\t]*//' > rg_lines.txt
        number_of_rg_lines=`wc -l rg_lines.txt | awk '{print \$1}'`  
    
        if [[ \$number_of_rg_lines == 0 ]]; then
            echo "No RG information in input bam/cram ${bam_file}" && exit 105
        elif [[ \$number_of_rg_lines > 1 ]]; then
            echo "More than one RG lines in input bam/cram ${bam_file}" && exit 106
        fi
        # we need to escape the tabs in order to be able to inline it.
        rg_line="\$(sed 's/\\t/\\\\t/g' < rg_lines.txt)"

        samtools sort -n -O BAM -T tmp1 -l 0 ${samtools_threads_args} ${bam_file} > name_sorted.bam
        samtools bam2fq ${samtools_threads_args} name_sorted.bam | \\
        bwa mem -K 100000000 -R "\$rg_line" -M -p \$index_prefix ${bwa_threads_args} - | \\
        samtools fixmate ${samtools_threads_args} -m - - | \\
        samtools sort ${samtools_threads_args} -l 0 -O BAM -T tmp2 - | \\
        samtools calmd ${samtools_threads_args} -b - ${reference_fasta} | \\
        tee ${remapped_bam_file} | \\
        samtools index ${samtools_threads_args} - ${index_file}
        
        ## These two last samtools command is just to verify that the outfile is a valid .bam 
        samtools quickcheck ${remapped_bam_file}
        samtools view -bS ${samtools_threads_args} ${remapped_bam_file} > /dev/null
    
        rm -Rf tmp1 tmp2 ref_cache ref_path bwa_index rg_lines.txt name_sorted.bam
        """
}