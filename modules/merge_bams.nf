process merge_bams {
	// This process merges multiple BAM files for a given sample using samtools merge.
	// If only one BAM file is provided, it simply creates symbolic links to the BAM
	// file and its index. If multiple files are provided, it merges them, creates an
	// index, and performs integrity checks using samtools. The merged BAM file and
	// its index are then outputted. The process also handles parallelization for
	// samtools operations based on available CPU resources.
	
	cpus 8

	publishDir "${params.results_dir}/bams/${sample_id}/", mode: 'copy'

	input:
		tuple val(sample_id), path(bam_files), path(bam_file_indexes)

	output:
		tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

	script:
		samtools_thread_args = task.cpus <= 1 ? "" : "-@ " + (task.cpus - 1)
		
		"""
		set -euxo pipefail
		count=0
		for bam_file in ${bam_files.collect { bam -> '"' + bam + '"' }.join(" ")}
		do
			echo \$bam_file >> file.list
			(( ++count ))
		done
		if [ \$count -eq 1 ]; then
			ln -s ${bam_files[0]} ${sample_id}.bam
			ln -s ${bam_file_indexes[0]} ${sample_id}.bam.bai
		else
			samtools merge ${samtools_thread_args} ${sample_id}.bam `cat file.list | sort`
			samtools index ${samtools_thread_args} ${sample_id}.bam
			
			### paranoia checks:
			samtools quickcheck ${sample_id}.bam
			samtools view -bS ${sample_id}.bam > /dev/null
		fi
		rm file.list
		"""
}
