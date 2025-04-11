process variant_annotation {
    // This process annotates a VCF file using SnpEff for variant annotation, 
    // adds CDS and region information using bcftools annotate, and outputs 
    // the annotated VCF and its index. It also handles special cases like 
    // setting VQSLOD to "Infinity" if it contains inf.

    input:
        tuple path(in_vcf), path(in_vcf_index)

    output:
        tuple path(out_vcf), path(out_vcf_index)

    script:
        snpeff_jar_path       = params.snpeff_jar_path
        snpeff_db_name        = params.variant_annotation_snpeff_db_name
        config_file           = params.variant_annotation_config_file
        snpeff_extra_options  = params.variant_annotation_snpeff_extra_options
        cds_gff_fn            = params.variant_annotation_cds_gff_fn
        annotations_header_fn = params.variant_annotation_header_fn
        regions_bed_fn        = params.variant_annotation_regions_bed_fn

        jvm_memory = task.memory.mega-200
        jvm_args   = "-Xmx${jvm_memory}m"

        base_name = in_vcf.getName().replaceAll(/(\.recal)?\.vcf(\.gz)?$/,"")
        if (base_name == in_vcf.getName()) {
            throw new Exception("the input vcf does not finish with the standard extensions .vcf or .vcf.gz")
        }
        out_vcf = base_name + ".ann.vcf.gz"
        out_vcf_index = out_vcf + ".tbi"
        """
        java -jar ${jvm_args} \\
            ${snpeff_jar_path} \\
            ${snpeff_db_name} \\
            -config ${config_file} \\
            ${snpeff_extra_options} \\
            ${in_vcf} | \\
        bcftools annotate \\
            -a ${cds_gff_fn} \\
            -c CHROM,FROM,TO,CDS \\
            -h ${annotations_header_fn} | \\
        bcftools annotate -a ${regions_bed_fn} \\
            -c CHROM,FROM,TO,RegionType | \\
        sed -r 's/VQSLOD=([-+]?)inf[^;]*(;?.*)\$/VQSLOD=\\1Infinity\\2/' | \\
        bgzip > ${out_vcf}

        bcftools index --tbi ${out_vcf}

        """
}
