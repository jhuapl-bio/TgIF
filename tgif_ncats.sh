#!/bin/bash
# author:	Robert Player

#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF

This tgif algorithm is tailored for read data from a modified NCATS enriched library.

NOTES:
	- all reference files (-r and -i) should have only 2 lines per sequence (i.e. no newlines [\n or \r] within the sequence)
	- fastq input should be from an NCATS enriched library, however may be from a WGS sheared library
	- all 'process comments' output to command line are coming from stderr output (i.e. will not redirect with stdout '>')
	- searching for valleys in alignment pileup
	- there is still softclipping happening in alignments for those reads from the nCATS enrichment
	- output tsv is in directory where -f file exists 'tgif_ncats-\$(basename \$f)'

DEPENDENCIES:
	GNU Parallel
	minimap2
		please symlink the minimap2 to 'tgif/bin/' (see Installation section of README)
		e.g. ln -s $fullpathto/minimap2/minimap2 $fullpathto/tgif/bin/

HELP/OUTFMT:
	-h      help	show this message
	-o		format	format of tsv output
REQUIRED:
	-t	INT	number of threads to GNU parallel over
	-f	READS	sequencing reads fasta/q file run NCATS enriched library
	-r	FASTA	fasta reference of target organism
OPTIONAL (but HIGHLY recommended):
	-i	FASTA	fasta of plasmid/inserted gene(s)
	-s	y/n		output sorted bam of downselected reads for IGV use [n]
	-p	y/n		generate read pileup png with R/ggplot2 per insertion site [n]

USAGE:
f="/data/project/reads.fastq"
r="/data/project/org_reference.fna"
i="/data/project/vector.fasta"
bash tgif_ncats.sh -t 10 -f "\$f" -r "\$r" -i "\$i"
# the output directory will be '/data/project/tgif_ncats-reads.fastq/'

EOF
}

# format of main output from this script
outfmt()
{
cat << EOF
TGIF OUTPUT FORMAT (insertions_filtered.tgif), line 1 contains column headers
	col1	HEADER			header(s), from reference [-r] input file
	col2	GAP_START		position in reference where gap starts (5' end of gap)
	col3	PREGAP_DEPTH	depth of ("position in col2" - 1)
	col4	PREGAP_PSTRANDEDNESS	comma separated values: positive strand count, total pregap read count, pos proportion
	col5	PREGAP_NSTRANDEDNESS	comma separated values: negative strand count, total pregap read count, neg proportion
	col6	GAP_END			position in reference where gap ends (3' end of gap)
	col7	POSTGAP_DEPTH	depth of ("position in col6" + 1)
	col8	POSTGAP_PSTRANDEDNESS	comma separated values: positive strand count, total postgap read count, pos proportion
	col9	POSTGAP_NSTRANDEDNESS	comma separated values: negative strand count, total postgap read count, neg proportion
	col10	GAP_LENGTH		length (bp) of deletion (gap, depth zero to zero; col3-col2)
	col11	CONFIDENCE		[low], [medium], or [high] probability of being an insertion site
								low - sites where at least one flanking depth is >1
								medium - flanking depths are both >1
								high - flanking depths are both >1, and both flanks have 100% opposing strandedness, e.g. if one flank has 100% of reads aligning to the +strand, then 100% of reads from the other align to the -strand (depending on depth, these are almost definitively insertion sites)
						*if using shear (WGS) data, 'high' confidence is not as meaningful

EOF
}


#	DEFAULTS & INPUTS & CHECKS
#===============================================================================
#	notes:
#		echo $? (0 = successful execution)
# absolute path to script dir
absolute_path_x="$(readlink -fn -- "$0"; echo x)"
absolute_path_of_script="${absolute_path_x%x}"
scriptdir=$(dirname "$absolute_path_of_script")
if [[ $? != 0 ]]; then
	>&2 echo "Please locate the function 'dirname' and symlink it to this script's bin."
	>&2 echo "example: ln -s /usr/bin/dirname /full/path/to/cinder/bin/dirname"
	exit
fi
bin="$scriptdir/bin"

# parse args
while getopts "hot:f:r:i:s:p:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		o) outfmt; exit 1 ;;
		t) THREADS=$OPTARG ;;
		f) FASTFILE=$OPTARG ;;
		r) REF=$OPTARG ;;
		i) INSERT=$OPTARG ;;
		s) SAMTOOLS=$OPTARG ;;
		p) PLOT=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# check args
if [[ -z "$THREADS" ]]; then printf "%s\n" "Please specify number of threads (-t)."; exit; fi
if [[ -z "$FASTFILE" ]]; then printf "%s\n" "Please specify input fastq (-f)."; exit; fi
if [[ ! -f "$FASTFILE" ]]; then printf "%s\n" "The input (-f) $FASTFILE file does not exist."; exit; fi
if [[ -z "$REF" ]]; then printf "%s\n" "Please specify reference fasta (-r)."; exit; fi
if [[ ! -f "$REF" ]]; then printf "%s\n" "The input (-r) $REF file does not exist."; exit; fi
# make sure user knows that running without the insert sequence as a filter is NOT advised
if [[ -z "$INSERT" ]]; then
	>&2 echo "You did not specify an insert sequence for filtering (-i)."
	>&2 echo "This is not advised, and will likely lead to incorrect insertion point identification."
	read -p "Do you want to continue? [y/n] " cont
	if [[ "$cont" == "y" ]]; then >&2 echo "too bad, terminating"; exit; fi
	if [[ "$cont" == "n" ]]; then >&2 echo "terminating"; exit; fi
	if [[ "$cont" != "y" && "$cont" != "n" ]]; then >&2 echo "that's not an expected answer, terminating"; exit; fi
else
	if [[ ! -f "$INSERT" ]]; then printf "%s\n" "The input (-i) '$INSERT' file does not exist."; exit; fi
fi
if [[ -z "$SAMTOOLS" ]]; then SAMTOOLS="n"; fi
case "$SAMTOOLS" in
	y) echo "testing" > /dev/null;;
	n) echo "testing" > /dev/null;;
	*) SAMTOOLS="n";;
esac
if [[ -z "$PLOT" ]]; then PLOT="n"; fi
case "$PLOT" in
	y) echo "testing" > /dev/null;;
	n) echo "testing" > /dev/null;;
	*) PLOT="n";;
esac




# setup other variables
runtime=$(date +"%Y%m%d%H%M%S%N")
#outdir="$scriptdir/tmp-$runtime"
bni=$(basename "$INSERT")
bnr=$(basename "$REF")
bn=$(basename "$FASTFILE")
fastdir=$(dirname "$FASTFILE")



# move outputs to FASTFILE input dir
outdir="$fastdir/tgif_ncats-$bn"




# only run if the 'mv' directory does not already exist
#	this assumes the '.dep' files are already present for filtering
#	please run 'rm -r "$fastdir/tgif_ncats-$bn"' if you want to re-run alignments etc 
if [[ ! -d "$outdir" ]]; then
	mkdir -p "$outdir/alignments"
	mkdir -p "$outdir/sortmp"
fi

>&2 echo "Tg finding $bni (-i) in $bnr (-r) with $bn (-f)"
echo "Tg finding $bni (-i) in $bnr (-r) with $bn (-f)" > "$outdir/log"


# get total reads in FASTFILE or exit
if [[ ! -f "$outdir/reads.fasta" ]]; then
	check=$(head -c1 "$FASTFILE")
	if [[ "$check" == ">" ]]; then		# is fasta
		# copy fasta input into out (fix headers to reduce file size)
		sed 's/ .*//' "$FASTFILE" > "$outdir/reads.fasta"
		total_reads=$(awk 'END{print(NR/2)}' "$outdir/reads.fasta")
	elif [[ "$check" == "@" ]]; then	# is fastq
		>&2 echo "converting fastq2fasta"
		# and convert to fasta (fix headers to reduce file size)
		$bin/fastq2fasta.sh -i "$FASTFILE" > "$outdir/reads.tmp"
		sed 's/ .*//' "$outdir/reads.tmp" > "$outdir/reads.fasta"; rm "$outdir/reads.tmp"
		total_reads=$(awk 'END{print(NR/2)}' "$outdir/reads.fasta")
	else
		printf "Fatal Error: $FASTFILE does not appear to be fasta for fastq formatted! Exiting."
		exit
	fi
else
	>&2 echo "counting reads..."
	total_reads=$(awk 'END{print(NR/2)}' "$outdir/reads.fasta")
fi




#	concatenating INSERT and REF file for mapping
target="$outdir/combined_insert_ref.fa"
cp "$REF" "$target"
cat "$INSERT" >> "$target"
insert_name=$(head -1 "$INSERT" | sed 's/^>//')
#	map reads to insertion sequence
query="$outdir/reads.fasta"
bnt=$(basename "$target")
if [[ ! -f "$target.idx" ]]; then	# make index of target
	echo "indexing $bnt" >> "$outdir/log"
	>&2 echo "indexing $bnt"
	$bin/minimap2 -t "$THREADS" "$target" -d "$target.idx" 2> /dev/null
fi
echo "mapping $total_reads reads to $bnt" >> "$outdir/log"
>&2 echo "mapping $total_reads reads to $bnt"
if [[ ! -f "$outdir/alignments/reads_to_both.paf" ]]; then
	$bin/minimap2 -x map-ont -t "$THREADS" "$target.idx" "$query" > "$outdir/alignments/reads_to_both.paf" 2> /dev/null
fi
uniq_count=$(awk -F'\t' '{i[$1]=1}END{print(length(i))}' "$outdir/alignments/reads_to_both.paf")
percent=$(printf "$uniq_count" | awk -v total_reads="$total_reads" '{printf("%.2f",100*($0/total_reads))}')
echo "	mapped $uniq_count reads uniquely" >> "$outdir/log"
echo "	alignment rate = ${percent}%" >> "$outdir/log"
>&2 echo "	mapped $uniq_count reads"
>&2 echo "	alignment rate = ${percent}%"

# get reads mapping to both a target genome sequence AND the insert sequence
echo "	finding reads aligning to both (-i) and (-r)" >> "$outdir/log"
>&2 echo "	finding reads aligning to both (-i) and (-r)"
grep -P "\t$insert_name\t" "$outdir/alignments/reads_to_both.paf" > "$outdir/alignments/reads_to_i.paf"
i_count=$(awk -F'\t' '{i[$1]=1}END{print(length(i))}' "$outdir/alignments/reads_to_i.paf")
echo "	mapped ${i_count} to (-i)" >> "$outdir/log"
>&2 echo "	mapped ${i_count} to (-i)"
grep -vP "\t$insert_name\t" "$outdir/alignments/reads_to_both.paf" > "$outdir/alignments/reads_to_r.paf"
r_count=$(awk -F'\t' '{i[$1]=1}END{print(length(i))}' "$outdir/alignments/reads_to_i.paf")
echo "	mapped ${r_count} to (-r)" >> "$outdir/log"
>&2 echo "	mapped ${r_count} to (-r)"
# index insert aligned read headers (col1), then only print alignments of these reads found in reference
awk -F'\t' '{if(NR==FNR){i[$1]=$0}else{if($1 in i){print($0)}}}' "$outdir/alignments/reads_to_i.paf" "$outdir/alignments/reads_to_r.paf" > "$outdir/alignments/overlap.paf"


#	DOWNSELECT ALIGNMENTS
# get single best alignment per read based on highest MAPQ
echo "downselecting overlap reads (MAPQ>=30)..." >> "$outdir/log"
>&2 echo "downselecting overlap reads (MAPQ>=30)..."
# on MAPQ>=30
awk -F'\t' '{if($12>=30){print($0)}}' "$outdir/alignments/overlap.paf" > "$outdir/alignments/downselect_q30.paf"
q30_count=$(wc -l <"$outdir/alignments/downselect_q30.paf")
echo "	$q30_count reads with MAPQ >=30" >> "$outdir/log"
>&2 echo "	$q30_count reads with MAPQ >=30"
# single best alignment per read on MAPQ
awk -F'\t' '{if($12>=30){if($1 in b){if($12>best[$1]){best[$1]=$12; b[$1]=$0;}}else{best[$1]=$12; b[$1]=$0;}}}END{for(r in b){print(b[r]);}}' "$outdir/alignments/downselect_q30.paf" > "$outdir/alignments/downselect.paf"
best_count=$(wc -l <"$outdir/alignments/downselect.paf")
if [[ "$best_count" -gt 0 ]]; then
	echo "	found single best alignments for $best_count reads" >> "$outdir/log"
	>&2 echo "	found single best alignments for $best_count reads"
else
	echo "	found $best_count reads with MAPQ>=30, exiting" >> "$outdir/log"
	>&2 echo "	found $best_count reads with MAPQ>=30, exiting"
	exit
fi




if [[ "$SAMTOOLS" == "y" ]]; then
	# added 20201124; generate sam file for downselected reads, 1 per ref (-i and -r)
	query="$outdir/downselected_reads.fasta"
	cut -f1 "$outdir/alignments/downselect.paf" | while read readheader; do
		grep -A1 -m1 "^>$readheader" "$outdir/reads.fasta"
	done > "$query"

	$bin/minimap2 -x map-ont -a -t "$THREADS" "$REF" "$query" > "$outdir/alignments/downselected_reads_to_r.sam" 2> /dev/null
	$bin/minimap2 -x map-ont -a -t "$THREADS" "$INSERT" "$query" > "$outdir/alignments/downselected_reads_to_i.sam" 2> /dev/null
	# convert to bam and sort
	#samtools="/data/apps/bin/samtools-1.2"
	if [[ ! -f "$outdir/alignments/downselected_reads_to_r.sorted.bam.bai" ]]; then
		$bin/samtools view -@ "$mT" -hb "$outdir/alignments/downselected_reads_to_r.sam" > "$outdir/alignments/downselected_reads_to_r.bam"
		# index the sam for tablet/IGV
		$bin/samtools sort -@ "$mT" "$outdir/alignments/downselected_reads_to_r.bam" "$outdir/alignments/downselected_reads_to_r.sorted"
		$bin/samtools index "$outdir/alignments/downselected_reads_to_r.sorted.bam"
	fi
	if [[ ! -f "$outdir/alignments/downselected_reads_to_i.sorted.bam.bai" ]]; then
		$bin/samtools view -@ "$mT" -hb "$outdir/alignments/downselected_reads_to_i.sam" > "$outdir/alignments/downselected_reads_to_i.bam"
		# index the sam for tablet/IGV
		$bin/samtools sort -@ "$mT" "$outdir/alignments/downselected_reads_to_i.bam" "$outdir/alignments/downselected_reads_to_i.sorted"
		$bin/samtools index "$outdir/alignments/downselected_reads_to_i.sorted.bam"
	fi
	# clean up
	rm "$outdir/alignments/downselected_reads_to_r.sam" 2> /dev/null
	rm "$outdir/alignments/downselected_reads_to_r.bam" 2> /dev/null

	rm "$outdir/alignments/downselected_reads_to_i.sam" 2> /dev/null
	rm "$outdir/alignments/downselected_reads_to_i.bam" 2> /dev/null
fi










# get depth of each position from PAF, in both directions (forward [5'-3'], then reverse [3'-5'])
# PAF
#	1	seqname
#	6	refname
#	7	reflength
#	8	ref align start position (1-based)
#	9	ref align end position (1-based)
# forward/reverse.dep
#	1	chr
#	2	position of gap (1-based) (3' side for forward, 5' side for reverse)
#	3	depth
#	4	number of immediately preceeding positions with depth=0 (i.e. gap length)

# forward depth file per unique reference sequence
fdep()
{
	bn=$(basename "$1")
	echo "	$bn" >> "$outdir/log"
	>&2 echo "	$bn"
	#		ANALYSIS DIRECTION: 5'-3'
	awk -F'\t' 'BEGIN{
		z=0;
	};{
		ilen[$6]=$7;
		ireads[$6]+=1;
		for(i=$8+1;i<=$9+1;i++){
			dep[$6][i]+=1;
		};
	}END{
		for(chr in ilen){
			for(p=1;p<=ilen[chr]+1;p++){
				if(dep[chr][p]!=""){
					printf("%s\t%s\t%s\t%s\n",chr,p,dep[chr][p],z);
					z=0;
				}else{
					z++;
				}
			}
		}
	}' "$1"
}
export outdir
export -f fdep

# reverse depth file per unique reference sequence
rdep()
{
	bn=$(basename "$1")
	echo "	$bn" >> "$outdir/log"
	>&2 echo "	$bn"
	#		ANALYSIS DIRECTION: 3'-5'
	awk -F'\t' 'BEGIN{
		z=0;
	};{
		ilen[$6]=$7;
		ireads[$6]+=1;
		for(i=$8+1;i<=$9+1;i++){
			dep[$6][i]+=1
		};
	}END{
		for(chr in ilen){
			for(p=ilen[chr]+1;p>=1;p--){
				if(dep[chr][p]!=""){
					printf("%s\t%s\t%s\t%s\n",chr,p,dep[chr][p],z);
					z=0;
				}else{
					z++;
				}
			}
		}
	}' "$1"
}
export outdir
export -f rdep






# split paf per unique reference sequence name
cut -f6 "$outdir/alignments/downselect.paf" | sort | uniq | while read rname; do
	awk -F'\t' -v rname="$rname" '{if($6==rname){print($0)}}' "$outdir/alignments/downselect.paf" > "$outdir/alignments/downselect-$rname.paf"
done
# parallel forward/reverse depth functions (fdep/rdep) per unique sequence in -r
find "$outdir/alignments" -name "downselect-*.paf" -exec du -b {} + | sort -r | cut -f2 > "$outdir/alignments/parallel.dep"
echo "generating 5'-3' depth file for..." >> "$outdir/log"
>&2 echo "generating 5'-3' depth file for..."
parallel --arg-file "$outdir/alignments/parallel.dep" --jobs="$THREADS" fdep > "$outdir/alignments/forward.dep"
echo "generating 3'-5' depth file for..." >> "$outdir/log"
>&2 echo "generating 3'-5' depth file for..."
parallel --arg-file "$outdir/alignments/parallel.dep" --jobs="$THREADS" rdep > "$outdir/alignments/reverse.dep"







echo "finding all insertion sites from '.dep' files..." >> "$outdir/log"
>&2 echo "finding all insertion sites from '.dep' files..."
# detect inconsistencies and print most probable gapped Tg insertion locations
# example data/info with 20200212 (basically a positive control dataset)
#	awk -F'\t' '{if($4!=0){if($3>1){print($0)}}}' forward.dep
#	Chr01	7080341	121	7080340
#	Chr07	2128032	6	27
#	Chr08	33451905	18	365
#	awk -F'\t' '{if($4!=0){if($3>1){print($0)}}}' reverse.dep
#	Chr01	7081470	114	6236
#	Chr07	2128004	11	27
#	Chr08	33451539	8	365

# LOGIC
#	check if each forward.dep gap location has a match in reverse.dep
# 		find rows where gap length != 0 and flanking depths are both >0
#		these are ALL potential insertions, likely most will have depth=1 flanking positions (i.e. col 3 and 5 both == 1)
#		and the gap length is >100kb (possibly? not sure about this)
awk -F'\t' '{
	if(NR==FNR){
		if($4!=0){
			if($3>0){
				fpos[$1]=$2; fdep[$1][$2]=$3; fgap[$1][$2]=$4;
			}
		}
	}else{
		if($4!=0){
			if($3>0){
				rdep[$1][$2]=$3;
				rgap[$1][$2]=$4;
			}
		}
	}
}END{
	for(chr in fpos){
		for(pos in fdep[chr]){
			expected_rpos=pos-(fgap[chr][pos]+1);
			if(fgap[chr][pos]==rgap[chr][expected_rpos]){
				printf("%s\t%s\t%s\t%s\t%s\t%s\n",chr,expected_rpos,rdep[chr][expected_rpos],pos,fdep[chr][pos],fgap[chr][pos]);
			}
		}
	}
}'  "$outdir/alignments/forward.dep" "$outdir/alignments/reverse.dep" > "$outdir/insertions_all.tsv"
# exit if file is empty
if [[ ! -s "$outdir/insertions_all.tsv" ]]; then
	echo "	no probable insertion sites found, check your .dep files if you think this is a mistake" >> "$outdir/log"
	>&2 echo "	no probable insertion sites found, check your .dep files if you think this is a mistake"
	exit
fi

# filter based on depth evidence
echo "filter1, checking depths of gap flanks..." >> "$outdir/log"
>&2 echo "filter1, checking depths of gap flanks..."
# INTERMEDIATE OUTPUT FORMAT (Tg insertions are assumed to be homozygous)
#	col1	SEQ				chr (or ref seq header)
#	col2	GAP_START		position in reference where gap starts (5' end of gap, where depth drops precipitously to zero) {from reverse.dep}
#	col3	PREGAP_DEPTH	depth of 'position in col2 minus 1'
#	col4	GAP_END			position in reference where gap ends (3' end of gap, where depth jumps dramatically) {from forward.dep}
#	col5	POSTGAP_DEPTH	depth of 'position in col4 plus 1'
#	col6	GAP_LENGTH		length (bp) of deletion (gap, depth zero to zero; col3-col2)

########################	OLD START	#################################
#	col7	CONFIDENCE		[wee] or [low] probability of being an insertion site
#								wee - probably not worth following up (if only one flanking position has >2 depth)
#								low - sites where both flanking depths are >2
########################	OLD END		#################################

# 20210420 change
#	col7	CONFIDENCE		[wee] or [low] probability of being an insertion site
#								wee - both flanking positions has depth==1
#								low - sites where at least one flanking depth is >1
#								experimental evidence shows that even a site where one of the flanking depths == 1 is a validated insertion site (20201029_15-16_merge.fastq)

# assign [wee] and [low] prob to sites
awk -F'\t' '{
	# wee prob sites
	if($3==1){
		if($5==1){
			printf("%s\t%s\n",$0,"wee"); next;
		}
	};
	# low prob sites (will include medium and high prob sites at this point)
	#	5-prime evidence
	if($3>1){printf("%s\t%s\n",$0,"low"); next};
	# 	3-prime evidence
	if($5>1){printf("%s\t%s\n",$0,"low"); next};
}' "$outdir/insertions_all.tsv" > "$outdir/insertions_filter1.tsv"
# exit if file is empty
if [[ ! -s "$outdir/insertions_filter1.tsv" ]]; then
	echo "	no probable insertion sites after filter1, exiting" >> "$outdir/log"
	>&2 echo "	no probable insertion sites after filter1, exiting"
	exit
fi


echo "filter2, checking strandedness of gap flanks..." >> "$outdir/log"
>&2 echo "filter2, checking strandedness of gap flanks..."
# grab reads from flanks of low prob sites only (filtering out 'wee' prob sites), and check strandedness
# 5' flank (target end ($9) + 1)
#	testing
#		pos reads
#while read f; do gapstart=$(printf "$f" | cut -f2); awk -F'\t' -v gs="$gapstart" '{if(gs==$9+1){print($0)}}' alignments/downselect.paf; done < insertions_filter1.tsv
#		neg reads
#while read f; do gapend=$(printf "$f" | cut -f4); awk -F'\t' -v ge="$gapend" '{if(ge==$8+1){print($0)}}' alignments/downselect.paf; done < insertions_filter1.tsv
awk -F'\t' '{
	if(FNR==NR){
		a[NR]=$1; b[NR]=$2; c[NR]=$3;
		d[NR]=$4; e[NR]=$5; f[NR]=$6;
	}else{
		for(x in b){if(b[x]==$9+1){if($5=="+"){tot1[x]++; pos1[x]++}; if($5=="-"){tot1[x]++; neg1[x]++}}};
		for(x in d){if(d[x]==$8+1){if($5=="+"){tot2[x]++; pos2[x]++}; if($5=="-"){tot2[x]++; neg2[x]++}}};
	}
}END{
	for(x in b){
		if(pos1[x]==""){pos1[x]=0}; if(tot1[x]==""){tot1[x]=0; ppr1[x]=0}else{ppr1[x]=pos1[x]/tot1[x]};
		if(pos2[x]==""){pos2[x]=0}; if(tot2[x]==""){tot2[x]=0; ppr2[x]=0}else{ppr2[x]=pos2[x]/tot2[x]};
		if(neg1[x]==""){neg1[x]=0}; if(tot1[x]==""){tot1[x]=0; npr1[x]=0}else{npr1[x]=neg1[x]/tot1[x]};
		if(neg2[x]==""){neg2[x]=0}; if(tot2[x]==""){tot2[x]=0; npr2[x]=0}else{npr2[x]=neg2[x]/tot2[x]};
		printf("%s\t%s\t%s\t%s,%s,%.4f\t%s,%s,%.4f\t%s\t%s\t%s,%s,%.4f\t%s,%s,%.4f\t%s\n", a[x], b[x], c[x], pos1[x], tot1[x], ppr1[x], neg1[x], tot1[x], npr1[x], d[x], e[x], pos2[x], tot2[x], ppr2[x], neg2[x], tot2[x], npr2[x], f[x]);
	}
}' <(grep -v "wee$" "$outdir/insertions_filter1.tsv") "$outdir/alignments/downselect.paf" > "$outdir/insertions_filter2.tsv"
# check output
#cat "$outdir/insertions_filter2.tsv"


echo "generating filtered list with confidence assignments..." >> "$outdir/log"
>&2 echo "generating filtered list with confidence assignments..."
# add to final output, with [medium] and [high] confidence
#	get lowest prob sites
printf "HEADER\tGAP_START\tPREGAP_DEPTH\tPREGAP_PSTRANDEDNESS\tPREGAP_NSTRANDEDNESS\tGAP_END\tPOSTGAP_DEPTH\tPOSTGAP_PSTRANDEDNESS\tPOSTGAP_NSTRANDEDNESS\tGAP_LENGTH\tCONFIDENCE\n" > "$outdir/insertions_filtered.tgif"
#	check strandedness info, assign [low], [medium], or [high] confidence
awk -F'\t' '{

	# get proportion from each strand for each gap flank
	split($4,p1,",");
	split($5,n1,",");
	split($8,p2,",");
	split($9,n2,",");

	# high probability insertion site
	if(p1[3]>=1){if(n2[3]>=1){printf("%s\t%s\n",$0,"high"); next}};
	if(n1[3]>=1){if(p2[3]>=1){printf("%s\t%s\n",$0,"high"); next}};

	# medium probability insertion site (both flank depth >2
	if($3>2){if($7>2){printf("%s\t%s\n",$0,"medium"); next}};

	# if no conditions above are met, print line (low confidence)
	printf("%s\t%s\n",$0,"low");

}' "$outdir/insertions_filter2.tsv" >> "$outdir/insertions_filtered.tgif"
# OUTPUT FORMAT (insertions_filtered.tgif)
#	col1	HEADER			header (from reference [-r] input file)
#	col2	GAP_START		position in reference where gap starts (5' end of gap, where depth drops precipitously to zero) {from reverse.dep}
#	col3	PREGAP_DEPTH	depth of 'position in col2' minus 1
#	col4	PREGAP_PSTRANDEDNESS	comma separated values: positive strand count, total pregap read count, pos proportion
#	col5	PREGAP_NSTRANDEDNESS	comma separated values: negative strand count, total pregap read count, neg proportion
#	col6	GAP_END			position in reference where gap ends (3' end of gap, where depth jumps dramatically) {from forward.dep}
#	col7	POSTGAP_DEPTH	depth of 'position in col4' plus 1
#	col8	POSTGAP_PSTRANDEDNESS	comma separated values: positive strand count, total postgap read count, pos proportion
#	col9	POSTGAP_NSTRANDEDNESS	comma separated values: negative strand count, total postgap read count, neg proportion
#	col10	GAP_LENGTH		length (bp) of deletion (gap, depth zero to zero; col3-col2)
#	col11	CONFIDENCE		[low], [medium], or [high] probability of being an insertion site
#								low - sites where at least one flanking depth is >2
#								medium - flanking depths are both >2
#								high - flanking depths are both >2, and both flanks have 100% opposing strandedness, e.g. if one flank has 100% of reads aligning to the +strand, then 100% of reads from the other align to the-strand (depending on depth, these are almost definitively insertion sites)
#						*if using shear (WGS) data, 'high' confidence is not meaningful



if [[ "$PLOT" == "y" ]]; then
	# get reads per group for plotting dataframe
	bash $scriptdir/tgif_plot.sh "$outdir"
fi


# print whitespace to log, print primary output to screen
echo "" >> "$outdir/log"
>&2 echo ""
cat "$outdir/insertions_filtered.tgif"













