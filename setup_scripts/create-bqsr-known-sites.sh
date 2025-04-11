#!/usr/bin/env bash 
PF7_DIR="<INSERT PATH HERE>pf_70_build/pf_70_internal_release/vcf/"

zless -S $PF7_DIR/*.vcf.gz | awk '/^##/ { print $0; next} {exit}' > bqsr-known-sites-pf7-pass.vcf  

zless -S $PF7_DIR/*.vcf.gz | grep -m1 ^#CHROM | cut -f 1-8 >> bqsr-known-sites-pf7-pass.vcf

zless -S "<INSERT PATH HERE>/pf_70_build/pf_70_internal_release/vcf/Pf3D7_*_v3.pf7.vcf.gz" | grep -v '^#' | awk '($7 == "PASS")' | cut -f 1-8 >> bqsr-known-sites-pf7-pass.vcf
