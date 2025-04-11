process bcftools_targets {
    // This process uses bcftools to extract variants from a VCF within specific
    // target intervals, writing a filtered VCF and its index for that region.

    tag { intervals.getBaseName() }

    input:
        tuple path(intervals), path(vcf), path(vcf_index)

    output:
        tuple path(intervals), path(out_vcf), path(out_index)

    script:

    out_vcf = intervals.getBaseName() + ".bcftools.vcf.gz"
    out_index = out_vcf + ".tbi"

    """
    set -exo pipefail

    cat ${intervals} | awk '{ print \$1"\\t"(\$2 + 1)"\\t"\$3"\\t"(\$2 == 1 ? 1 : \$2 - ${params.db_import_upstream_padding}) }' > intervals.txt

    bcftools view \\
        --targets-file intervals.txt \\
        --output-type z \\
        --output-file ${out_vcf} \\
        ${vcf}

    bcftools index -t ${out_vcf}

    """
}
