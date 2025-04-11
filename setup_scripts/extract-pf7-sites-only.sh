OUTDIR=$1
REF_IDX=$2
FILES=`find "<INSERT PATH HERE>/pf_70_build/pf_70_internal_release/vcf/" | grep ".vcf.gz$"`

mkdir -p $OUTDIR

for f in $FILES
do
	bn=$(basename $f .vcf.gz)
	cat > $OUTDIR/${bn}.sh <<EOF
        if [ -f "$OUTDIR/${bn}.done" ]; then
           exit 0
        fi
gunzip < $f | grep -v '^#' | cut -f 1-8 | awk '(\$7 == "PASS"){ print \$0 }' | gzip > $OUTDIR/${bn}.vcf.gz
        touch $OUTDIR/${bn}.done
EOF
	chmod a+x $OUTDIR/${bn}.sh
	bsub -R"select[mem>4000] rusage[mem=4000]"  -M4000 -J ${bn} -o $OUTDIR/%J.out -e $OUTDIR/%J.err $OUTDIR/${bn}.sh
done

