BUNDLE=$1

PF_PLASMDB_BASEURL=https://plasmodb.org/common/downloads
PF_PLASMDB_RELEASE=54
PF_CROSSES_BASEURL=ftp://ngs.sanger.ac.uk/production/malaria/pf-crosses/1.0/
PF_CROSSES=(7g8_gb4 hb3_dd2 3d7_hb3)
SNP_EFF_URL=https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

mkdir -p $BUNDLE

if false; then
wget $PF_PLASMDB_BASEURL/release-$PF_PLASMODB_RELEASE/Pfalciparum3D7/fasta/data/PlasmoDB-${PF_PLASMODB_RELEASE}_Pfalciparum3D7_Genome.fasta \
	-P $BUNDLE

for CROSS in "${PF_CROSSES[@]}"; do 	
  wget $PF_CROSSES_BASEURL/${CROSS}.combined.final.vcf.gz $PF_CROSSES_BASEURL/${CROSS}.combined.final.vcf.gz.tbi -P $BUNDLE
done

wget $PF_CROSSES_BASEURL/regions-20130225.onebased.txt -O - | \
	awk '{ print $1"\t"($2 - 1)"\t"$3"\t"$4 } ' \
	> $BUNDLE/regions-20130225-plus-api-mt.bed
echo -e	"Pf3D7_API_v3\t0\t34250\tApicoplast" >> $BUNDLE/regions-20130225-plus-api-mt.bed
echo -e "Pf_M76611\t0\t5967\tMitochondrion" >> $BUNDLE/regions-20130225-plus-api-mt.bed
grep -F "Core" < $BUNDLE/regions-20130225-plus-api-mt.bed > $BUNDLE/regions-20130225.core.bed 
gzip $BUNDLE/regions-20130225-plus-api-mt.bed


wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip -P $BUNDLE

unzip $BUNDLE/snpEff_latest_core.zip snpEff/snpEff.config -d $BUNDLE
mv $BUNDLE/snpEff/snpEff.config $BUNDLE
rm -Rf $BUNDLE/snpEff snpEff_latest_core.zip

fi

### Files that I don't know how to generate from scratch, so we simply copy over the existing ones:

cp '<INSERT PATH HERE>/resources/Pfalciparum_replace_Pf3D7_MIT_v3_with_Pf_M76611.gff.cds.gz' $BUNDLE
cp '<INSERT PATH HERE>/resources/annotations.hdr' $BUNDLE

