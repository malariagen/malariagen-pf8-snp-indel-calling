process process_reference {
	// This process creates the reference index and dict files, generates a bwa index file
	// (for alignment) and packages the index into a tarball, creates one .bed file 
	// per chromosome split by intervals, then split further using split_intervals script. 
	
    input:
		path(reference_fasta)
		val window_size
		val upstream_padding
	
	output:
		path fai_file                 , emit: fai_file
		path dict_file                , emit: dict_file
		path bwa_index_tarball        , emit: bwa_index_tarball
		path chromosome_interval_files, emit: chromosome_interval_files
		path db_import_interval_files , emit: db_import_interval_files
		path join_interval_files      , emit: join_interval_files
	
	script:
		base_name                 = reference_fasta.getBaseName()
		fai_file                  = reference_fasta + ".fai"

		dict_file                 = base_name + ".dict"
		bwa_index_tarball         = base_name + ".bwa.tar"
		chromosome_intervals      = base_name + ".chromosome_intervals.bed"
		chromosome_interval_files = base_name + ".chromosome_interval_files.list"
		db_import_intervals       = base_name + ".db_import_intervals.bed"
		db_import_interval_files  = base_name + ".db_import_interval_files.list"
		join_intervals            = base_name + ".join_intervals.bed"
		join_interval_files       = base_name + ".join_interval_files.list"

		"""
		set -x
		samtools dict ${reference_fasta} > ${dict_file} || exit 101
		samtools faidx ${reference_fasta} || exit 102
		mkdir index
		bwa index -a bwtsw -p index/${base_name} ${reference_fasta} || exit 103
		echo "${base_name}" > index/prefix.txt
		tar cf ${bwa_index_tarball} index/* --transform 's/^index\\///' || exit 104
		rm -Rf index

		cat > ./split_intervals <<-"EOF"
		#!/usr/bin/env perl
		use strict;
		use warnings;
		use File::Spec;
		my \$list = shift;
		my \$outlist = shift;
		my \$dir = shift;
	
		mkdir \$dir or die "could not create output dir \$dir";
		open (my \$intervals, "<", "\$list") or die "could not open intervals file \$list";
		open (my \$olist, ">" , "\$outlist") or die "could not open output list file \$outlist";
		my \$seq = 0;
		my \$last_chr = "";
		while (my \$line = <\$intervals>) {
			my (\$chr, \$start, \$end) = split(/\\t/, \$line);
			\$seq = \$last_chr eq \$chr ? \$seq + 1 : 1;
			open (my \$out, ">", "\$dir/\${chr}_\${seq}.bed") or die "could not open interval output file \${chr}_\$seq";
			print \$out \$line or die "could not write content in \${chr}_\$seq";
			print \$olist File::Spec->rel2abs( "\$dir/\${chr}_\${seq}.bed" ),"\\n" or die "could not write in output file list \$outlist";
			\$last_chr = \$chr;
			close \$out;
		}
		close \$olist or die "problem closing output list file \$outlist";
		close \$intervals;
		EOF

		chmod a+x ./split_intervals

		cat ${fai_file} | awk '{ print \$1"\\t0\\t"\$2 }' > ${chromosome_intervals}
		
		# TODO need to add bedtools to the singularity image so that we don't need to download or rely on modules

		wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary -O ./bedtools-static
		chmod a+x ./bedtools-static
		./bedtools-static makewindows -g ${fai_file} -w ${window_size} > ${join_intervals}
		cat ${join_intervals} | \
			awk '{
			if (\$2 > 500)
				print \$1"\\t"(\$2 - ${upstream_padding})"\\t"\$3;
			else 
				print \$1"\\t0\\t"\$3;
			}' \
			> ${db_import_intervals} 
		rm bedtools-static
		
		./split_intervals ${chromosome_intervals} ${chromosome_interval_files} chromosome_interval_files || exit 111
		./split_intervals ${join_intervals} ${join_interval_files} join_interval_files || exit 112
		./split_intervals ${db_import_intervals} ${db_import_interval_files} db_import_interval_files || exit 113
		"""
}