#!/bin/bash
# author:	Robert Player

#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF

This tgif module is used to design primers for validation of potential insertion sites identified by the primary tgif algorithm (tgif_ncats.sh).

DEPENDENCIES:
	GNU Parallel
		O. Tange (2011): GNU Parallel - The Command-Line Power Tool,
		;login: The USENIX Magazine, February 2011:42-47.
	Primer3
		please symlink 'primer3_core' to 'tgif/bin/' (see Installation section of README)
		http://primer3.org/manual.html

HELP/OUTFMT:
	-h      help	show this message
	-f		format	format of input
REQUIRED:
	-t	INT	number of threads to GNU parallel over
	-i	TGIF	full path to output from primary tgif algorithm (insertions_filtered.tgif)
	-r	FASTA	fasta reference of insert target organism (used to generate -i)
OPTIONAL:
	-g	INT		if GAP_LENGTH of insertion site is greater than INT bp (-g), no primers will be generated for the site [default 5000]

USAGE:
i="/data/project/tgif_ncats-reads.fastq/insertions_filtered.tgif"
r="/data/project/org_reference.fna"
bash tgif_primer3.sh -t 10 -i "\$i" -r "\$r"
EOF
}

# required format of input file (main output of primary script 'tgif_ncats.sh')
infmt()
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
								low - sites where at least one flanking depth is >2
								medium - flanking depths are both >2
								high - flanking depths are both >2, and both flanks have 100% opposing strandedness, e.g. if one flank has 100% of reads aligning to the +strand, then 100% of reads from the other align to the -strand (depending on depth, these are almost definitively insertion sites)
						*if using shear (WGS) data, 'high' confidence is not meaningful

EOF
}



p3()
{
	insert="$1"

	bn=$(basename "$insert")

	HEADER=$(cut -f1 "$insert")
	GAP_START=$(cut -f2 "$insert")
	GAP_END=$(cut -f6 "$insert")
	GAP_LENGTH=$(cut -f10 "$insert")
	SIZE_RANGE_MAX=$(printf "$GAP_LENGTH" | awk '{print($0+1000)}')

	# grab template sequence (500 bp upstream/downstream of gap start/end)
	temp_start=$(printf "$GAP_START" | awk '{print($0-500)}')
	temp_end=$(printf "$GAP_END" | awk '{print($0+500)}')

	#	if GAP_LENGTH is >$GAPCUTOFF bp, skip
	if [[ "$GAP_LENGTH" -gt "$GAPCUTOFF" ]]; then 
		echo "	skipping site at $HEADER between $temp_start and $temp_end: gap length of $GAP_LENGTH > 5000" >> "$workdir/log"
		>&2 echo "	skipping site at $HEADER between $temp_start and $temp_end: gap length of $GAP_LENGTH > 5000"
		return
	fi
	echo "	finding primers for site at $HEADER between $temp_start and $temp_end" >> "$workdir/log"
	>&2 echo "	finding primers for site at $HEADER between $temp_start and $temp_end"
	TEMPLATE=$(grep -A1 "^>$HEADER" "$REF" | tail -1 | cut -c${temp_start}-${temp_end})


# do NOT tab over EOF or anything between them...
cat > "$workdir/p3record.$bn" << EOF
SEQUENCE_ID=$HEADER-$GAP_START-$GAP_END
SEQUENCE_TEMPLATE=$TEMPLATE
SEQUENCE_TARGET=500,$GAP_LENGTH
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=16
PRIMER_MAX_SIZE=22
PRIMER_MAX_NS_ACCEPTED=1
PRIMER_PRODUCT_SIZE_RANGE=22-$SIZE_RANGE_MAX
P3_FILE_FLAG=1
SEQUENCE_INTERNAL_EXCLUDED_REGION=500,$GAP_LENGTH
PRIMER_EXPLAIN_FLAG=1
=
EOF

	cd $workdir
	$bin/primer3_core < "$workdir/p3record.$bn" > /dev/null


}
export bin workdir REF GAPCUTOFF
export -f p3






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
while getopts "hft:i:r:g:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		f) infmt; exit 1 ;;
		t) THREADS=$OPTARG ;;
		i) TGIF=$OPTARG ;;
		r) REF=$OPTARG ;;
		g) GAPCUTOFF=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# check args
if [[ -z "$THREADS" ]]; then printf "%s\n" "Please specify number of threads (-t)."; exit; fi
# tgif output as input
if [[ -z "$TGIF" ]]; then printf "%s\n" "Please specify a TGIF output for designing primers against (-i)."; exit; fi
if [[ ! -f "$TGIF" ]]; then printf "%s\n" "The input (-i) '$TGIF' file does not exist."; exit; fi
if [[ -d "$TGIF" ]]; then printf "%s\n" "The input (-i) '$TGIF' is a directory, exiting."; exit; fi
if [[ -z "$REF" ]]; then printf "%s\n" "Please specify reference fasta (-r)."; exit; fi
if [[ ! -f "$REF" ]]; then printf "%s\n" "The input (-r) '$REF' file does not exist."; exit; fi
if [[ -z "$GAPCUTOFF" ]]; then
	>&2 echo "a gap cutoff length (-g) was not specified, setting to 5000 bp"
	GAPCUTOFF=5000;
fi



# setup other variables
runtime=$(date +"%Y%m%d%H%M%S%N")
#outdir="$scriptdir/tmp-$runtime"
bni=$(basename "$TGIF")
bnr=$(basename "$REF")
indir=$(dirname "$TGIF")
workdir="$indir/primer3_files"
mkdir -p "$workdir"






# remove header of tgif output file
tail -n+2 "$TGIF" > "$workdir/tgif"
TOTAL_SITES=$(wc -l <"$workdir/tgif")
# split into 1 insert (row) per file for parallel jobs
split -l1 "$workdir/tgif" "$workdir/tgif_"
# make parallel input file and submit in parallel to p3 function
find "$workdir" -name "tgif_*" -exec du -b {} + | sort -r | cut -f2 > "$workdir/parallel.p3"
echo "finding primers with primer3_core for $TOTAL_SITES total sites:" > "$workdir/log"
>&2 echo "finding primers with primer3_core for $TOTAL_SITES total sites:"
parallel --arg-file "$workdir/parallel.p3" --jobs="$THREADS" p3












