ANNOVARDIR=../../opt/annovar
YENGODIR=../data/raw/sum_stats/yengo2018
TARGETDIR=../data/processed/annovar

awk '{print $2,$3,$3,$5,$4,$1,$7,$8,$9}' $YENGODIR/Meta-analysis_Locke_et_al+UKBiobank_2018_top_941_from_COJO_analysis_UPDATED.txt \
	> $TARGETDIR/yengo_941_variants.avinput

perl $ANNOVARDIR/annotate_variants.pl \
	-out $TARGETDIR/yengo_941_variants \
	-build hg19 \
	$TARGETDIR/yengo_941_variants.avinput \
	$ANNOVARDIR/humandb/
