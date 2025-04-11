#!/usr/bin/env nextflow
params.generate_bams = false
params.generate_gvcfs = false
params.joint_genotyping = false
params.help = false
params.bam_files = ''

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { bcftools_targets      } from './modules/bcftools_targets.nf'
include { bqsr                  } from './modules/bqsr.nf'
include { gatk_haplotype_caller } from './modules/gatk_haplotype_caller.nf'
include { genotype_gvcfs        } from './modules/genotype_gvcfs.nf'
include { get_file_from_irods   } from './modules/irods.nf'
include { import_gvcfs          } from './modules/import_gvcfs.nf'
include { mark_duplicates       } from './modules/mark_duplicates.nf'
include { merge_bams            } from './modules/merge_bams.nf'
include { process_reference     } from './modules/process_reference.nf'
include { samtools_stats        } from './modules/samtools_stats.nf'
include { remap_bam_file        } from './modules/remap_bam_file.nf'
include { variant_annotation    } from './modules/variant_annotation.nf'
include { variant_filtration    } from './modules/variant_filtration.nf'
include { variant_recalibration_apply_model } from './modules/variant_recalibration_apply_model.nf'
include { variant_recalibration_build_model as   variant_recalibration_build_snp_model } from './modules/variant_recalibration_build_model.nf'
include { variant_recalibration_build_model as variant_recalibration_build_indel_model } from './modules/variant_recalibration_build_model.nf'

// helper functions
def printHelp() {
    log.info """
    Usage:
        nextflow run main.nf (--generate_bams | --generate_gvcfs | --joint_genotyping) (--manifest [manifest path] | --bam_files [bam files directory path] | --gvcf_map_directory [gvcf maps directory path] ) --reference [ref genome fasta path] --results_dir [output dir]

    Description:
        Genotyping and variant calling pipeline for Plasmodium falciparum genomic data (part of the Pf8 malariaGEN project)

    Workflows:
        --generate_bams                   Run sample BAM generation steps of pipeline
        --generate_gvcfs		          Run gvcf generation steps of pipeline
        --joint_genotyping                Run joint genotyping

    Options:
        Inputs:
        --manifest <path>                 Path to sample manifest file (Must be used if --generate_bams is specified)
        --bam_files <path>                Path to sample bam files (Must be used if --generate_bams is not specified but --generate_gvcfs is specified)
        --gvcf_map_directory <path>       Path to gvcf map directory (Must be used if --joint_genotyping is specified)

        Outputs:
        --results_dir <path>              Path to a directory for output files

        Additional options:
        -N <email address>                For email notifications when the pipeline completes/fails        
        -resume                           If specified, pipeline can be restarted from last saved point
        -with-trace                       Enables the generation of a trace.txt file, containing detailed execution information for all processes in the pipeline. Useful for monitoring and debugging workflow performance. 
        -with-report                      Generates an HTML report that provides a high-level summary of the pipeline execution
        --help                            If true, print help information then quit (default: false)
    """.stripIndent()
}

def get_num_samples(file) {
    def samples = [] as Set
    def line_num = 0
    def index = -1
    File my_file = new File(file)

    my_file.splitEachLine("\t") { fields ->
        if (++line_num == 1) {
            index = fields.findIndexOf { it -> it == "sample" }
            if (index == -1)  {
                exit 1, 'The header of the given file does not contain the column "sample"'
            }
        } else {
            samples << fields[index]
        }
    }
    return samples.size()
}

// Subworkflows
workflow improve_lane_bams {

    take:
        reference_bundle 
        bqsr_known_sites
        sample_and_input_file
    
    main:
        reference_fa   = reference_bundle.fasta
        reference_fai  = reference_bundle.index
        reference_dict = reference_bundle.dictionary

        // Mapping
        remap_bam_file(reference_bundle.fasta, reference_bundle.bwa_index, sample_and_input_file)
        remapped_bam = remap_bam_file.out
        
        // Mark duplicates
        mark_duplicates(remapped_bam)
        mark_duplicated = mark_duplicates.out


        bqsr(reference_fa, reference_fai, reference_dict, 
            bqsr_known_sites, mark_duplicated)

    emit:
        bqsr.out.recalibrated
}

workflow generate_sample_bams {
    take:
        reference_bundle
        bqsr_known_sites
        sample_and_stagged_path

    main:
        // Run pipeline for every sample
        improve_lane_bams(reference_bundle, bqsr_known_sites, sample_and_stagged_path)

        grouped_by_sample_id = improve_lane_bams.out
            .groupTuple()

        // Merge BAM files
        merge_bams(grouped_by_sample_id)

        // Stats
        merge_bams_alignments =  merge_bams.out.map { o -> o[1] }
        samtools_stats(merge_bams_alignments)

    emit:
        merge_bams_out = merge_bams.out
        samtools_stats_out = samtools_stats.out
}

workflow variant_recalibration {
    take: 
        input_intervals_vcf_and_index

    main:
        vqsr_input_vcfs_list = input_intervals_vcf_and_index
            .map { o -> "" + o[1] }
            .toList()
        
        compose_file_list(vqsr_input_vcfs_list, "all-vqsr-input-vcfs.list")
        input_vcf_list_file = compose_file_list.out
        
        // Variant quality score recalibration (based on core regions of genome)
        // SNP model
        variant_recalibration_build_snp_model(  input_vcf_list_file, 'VQSR',   'SNP', params.vqsr_default_parameters_file)

        // Indel model
        variant_recalibration_build_indel_model(input_vcf_list_file, 'VQSR', "INDEL", params.vqsr_default_parameters_file)

        // Apply VQSR model
        variant_recalibration_apply_model(input_intervals_vcf_and_index,
            variant_recalibration_build_snp_model.out,
            variant_recalibration_build_indel_model.out,
            99, 99)

    emit:
        variant_recalibration_apply_model.out
}

workflow generate_sample_gvcfs {
    take:
        sample_data
        reference_fa
        reference_fai
        reference_dict
        chromosome_interval_file_list

    main:        
        chromosomes_and_samples = chromosome_interval_file_list
            .splitText()
            .map { it -> file(it.trim()) }
            .combine(sample_data)
            .map { it -> it.flatten() }
        
        chromosomes_and_samples
            | subscribe { it -> println "This is chromosomes and samples: ${it}"}
        
        gatk_haplotype_caller(chromosomes_and_samples, reference_fa, reference_fai, reference_dict)
        
        gatk_haplotype_caller.out
            .map { it -> [it[0].getBaseName().replaceAll(/_[0-9]+$/,'') , [it[1], "" + it[2]]] } 
            .groupTuple() // order is not important.
            .set { gvcf_maps_input }
        
        gatk_haplotype_caller.out
            .map { file -> [ file[2], file[3] ] }
            .set { output_gvcfs }
        
        compose_gvcf_map(gvcf_maps_input)
        gvcf_maps_by_chromosome = compose_gvcf_map.out

    emit:
        gvcf_maps_by_chromosome
        output_gvcfs
}

workflow joint_genotyping {
	take:
        gvcf_maps_by_chromosome
        reference_fa
        reference_fai
        reference_dict
        db_import_interval_file_list
        join_interval_file_list
        chromosome_interval_file_list

	main:
        db_import_interval_file_list
            .splitText(by: 1)
            .map { it -> file(it.trim()) }
            .set { db_import_interval_files }

        join_interval_file_list
            .splitText(by: 1)
            .map { it -> file(it.trim()) }
            .set { join_interval_files }
        
        gvcf_maps_by_chromosome.map { it -> [ it.getBaseName().replaceAll(/_[0-9]+$/,""), it ] }
            .cross( db_import_interval_files.map { it -> [ it.getBaseName().replaceAll(/_[0-9]+$/,"").replaceAll("db_import_","") , it ] })
            .map { it -> tuple(it[1][1], it[0][1]) }
            .set { db_import_intervals_and_maps }
        
        // Create intervals

        // Import GVCFs into GenomicDB
        import_gvcfs(db_import_intervals_and_maps)

        // Joint Genotyping
        genotype_gvcfs(import_gvcfs.out, reference_fa, reference_fai, reference_dict)

        genotype_gvcfs.out.map { it -> [ it[0].getBaseName(), it ] }
            .cross( join_interval_files.map { it -> [ it.getBaseName() , it ] })
            .map { it -> tuple(it[1][1], it[0][1][1], it[0][1][2]) }
            .set { genotype_gvcfs_with_unpadded_intervals }

        // BCFTools Target
        bcftools_targets(genotype_gvcfs_with_unpadded_intervals)

        bcftools_targets.out
            .map { it -> {
                    def chr = it[0].getBaseName().replaceAll(/_[0-9]+$/, "")
                    return [chr, [ it[0].getBaseName().replaceAll(/^.*_([0-9]+)$/, '$1').toInteger() , it[1] ]] 
                }
            }
            .groupTuple(sort : { it -> it[0] })
            .map { it -> [ it[0], it[1].collect{ it1 -> it1[1] }.flatten()  ] }
            .set { vcfs_sorted_and_grouped_by_chromosome }
    
        chromosome_interval_file_list
            .splitText()
            .map { it -> file(it.trim()) }
            .map { it -> [ it.getBaseName().replaceAll(/_([0-9]+)$/, "") , it ] }
            .cross( vcfs_sorted_and_grouped_by_chromosome )
            .map { it -> tuple(it[0][1], it[0][0], it[1][1]) } // <chr_int.bed, chr, list-of-vcfs >
            .set { chromosome_and_vcfs_sorted }
        
        gather_vcfs(chromosome_and_vcfs_sorted)

        variant_recalibration(gather_vcfs.out)
        recalibrated_vcfs = variant_recalibration.out.map { it -> [it[1], it[2]] } // discard the intervals file.

        variant_annotation(recalibrated_vcfs)
        annotated_vcfs = variant_annotation.out

        variant_filtration(annotated_vcfs, reference_fa, reference_fai, reference_dict)
        filtered_vcfs = variant_filtration.out 
    emit:
        out = filtered_vcfs

}

process compose_file_list {
    // This process creates a text file listing paths from file_list, names it
    // with out_name, and outputs it — useful for tools that require an input list
    // of files.

    input:
        val(file_list)
        val(out_name)

    output:
        path(out_full_name)

    script:
        out_full_name = "${out_name}" // implicitly is prefixed with the out-dir path.
        
        // Creating a list of file paths to write to the output
        file_list_lines = file_list.join('\n')

        """
        echo -e "${file_list_lines}" > ${out_full_name}
        """
}


process gather_vcfs {
    // This process gathers multiple VCFs into a merged, compressed VCF, using GATK GatherVcfs 
    // when needed, and ensures the output is indexed with tabix.

    memory '4G'

    input:
        tuple(path(intervals), val(out_name), val(sorted_list_of_vcfs))
    
    output:
        tuple path(intervals), path(out_vcf), path(out_vcf_idx) 

    script:
        out_vcf = "${out_name}.vcf.gz"
        out_vcf_idx = "${out_vcf}.tbi"

        if (sorted_list_of_vcfs.size() == 1) {
            """
            ln -s ${sorted_list_of_vcfs[0]} ${out_vcf}
            if [ -f ${sorted_list_of_vcfs[0]}.tbi ]; then
                ln -s ${sorted_list_of_vcfs[0]}.tbi ${out_vcf_idx}
            else 
                tabix ${out_vcf}
            fi
            """
        } else {
            args_file = "${out_vcf}.args"
            args_lines = sorted_list_of_vcfs.collect { "-I ${it}" }.join('\n')

            """
            echo -e "${args_lines}" > ${args_file}

            gatk GatherVcfs --java-options "-Xmx3g" \\
                --arguments_file ${args_file} \\
                -O ${out_vcf}
            tabix ${out_vcf}
            """
        }
}

process compose_gvcf_map {
    tag "$chromosome" 

    input:
        tuple(val(chromosome), val(list_of_samples_and_vcfs))

    output:
        path(map_file)

    script:
        map_file = "${chromosome}.gvcf_map"
        tmp_file = "${map_file}.tmp"

        lines = list_of_samples_and_vcfs.collect { "${it[0]}\t${it[1]}" }.join('\n')

        """
        echo -e "${lines}" > ${tmp_file}
        mv ${tmp_file} ${map_file}
        """
}


// Parameter checking
// Strategy is to list all the errors/problems with the command invocation at once, rather than one error at a time
// Variable to check if errors were produced
workflow initialising_pipeline {
    def errors = 0

    if (params.help) {
        printHelp()
        System.exit(0)
    }

    // Create results directory
    if (params.results_dir) {
        def results_path = file(params.results_dir)
        if (!results_path.exists() && !results_path.mkdir()) {
            System.exit(1), "Failed to make results dir: $results_path. Check you have permissions to create the dir."
        } else {
            println(results_path.exists() ? "Using existing results dir: $results_path" : "Successfully made results dir: $results_path")
        }
    } else {
        log.error("Please use --results_dir to supply a path to a directory that will store pipeline results (if the directory does not exist, it will be created)")
        errors = 1
    }

    // Validate reference genome
    if (!params.reference || !file(params.reference).exists()) {
        log.error(params.reference ? "The reference genome supplied using --reference cannot be found at the given path" : "Please use --reference to supply a path to a Plasmodium falciparum reference genome (FASTA format)")
        errors = 1
    }

    // Validate manifest file if required
    if (params.generate_bams) {
        if (!params.manifest || !file(params.manifest).exists()) {
            log.error(params.manifest ? "The manifest supplied using --manifest cannot be found at the given path" : "Please use --manifest to supply a path to a manifest file")
            errors = 1
        }
    }

    if (errors == 1) {
        System.exit(1)
    }
}


// Main entry-point workflow
workflow {
    initialising_pipeline()

    reference_fa = params.reference

    reference_outputs = process_reference(reference_fa, params.genotyping_interval_size, params.db_import_upstream_padding)

    fasta_index                   = reference_outputs.fai_file
    dict_file                     = reference_outputs.dict_file
    bwa_index                     = reference_outputs.bwa_index_tarball
    chromosome_interval_file_list = reference_outputs.chromosome_interval_files
    db_import_interval_file_list  = reference_outputs.db_import_interval_files
    join_interval_file_list       = reference_outputs.join_interval_files

    reference_bundle = [
        fasta      : reference_fa,
        index      : fasta_index,
        dictionary : dict_file,
        bwa_index  : bwa_index
    ]

    // Create channel of input sample bam files or parse a manifest file
    if (params.generate_bams) {
        files_from_irods = Channel.fromPath(params.manifest)
            .splitCsv(header:true, sep:"\t")
            .map { row -> tuple(row.sample, row.irods_path) }

        groupings_ch = Channel.fromPath(params.manifest)
            .splitCsv(header: true, sep:"\t")
            .map { row -> [ row.sample, row.irods_path ] }
            .groupTuple()
            .map { sample_key, file_list -> [ sample_key, groupKey( sample_key, file_list.size() ) ] }

        samples_and_irods_paths = groupings_ch.combine(files_from_irods, by:0)
            .map { it -> [ it.get(1), it.get(2) ] } 

        get_file_from_irods(samples_and_irods_paths)
        samples_and_stagged_paths = get_file_from_irods.out 

        bqsr_known_sites = tuple(params.bqsr_known_sites, params.bqsr_known_sites_index)           
        
        generate_sample_bams(reference_bundle, bqsr_known_sites, samples_and_stagged_paths)

        if (params.generate_gvcfs) {
            generate_sample_bams.out.merge_bams_out
                .set { bam_files_ch }
        }

    } else if (params.generate_gvcfs) {
        Channel.fromPath(params.bam_files, checkIfExists:true)
            .splitText()
            .map { it.trim().split("\t") }
            .map { row -> tuple(row[0], row[1], row[2]) }
            .set { bam_files_ch }
    }

    if (params.generate_gvcfs) {
        // Call variants
        generate_sample_gvcfs(bam_files_ch, reference_fa, fasta_index,
            dict_file, chromosome_interval_file_list)
    }

    if (params.joint_genotyping) {
        Channel.fromPath("${params.gvcf_map_directory}/*.gvcf_map")
            .set { gvcf_maps }

        joint_genotyping(gvcf_maps, reference_fa, fasta_index, dict_file,
            db_import_interval_file_list, join_interval_file_list, chromosome_interval_file_list)
    }

}
