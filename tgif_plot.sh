#!/bin/bash
# This script will call R and generate a pileup visualization (png) of reads per identified insertion site (from the output dir of the primary script 'tgif_ncats.sh') using ggplot2.
#	arg1 is the "$outdir" from primary TGIF script
#		eg. if (-f) was "/data/project/reads.fastq" for 'tgif_ncats.sh', to plot results:
#		bash tgif_plot.sh /data/project/tgif_ncats-reads.fastq/
if [[ "$1" == "" ]]; then echo "ARG1 is required"; exit; fi
if [[ ! -d "$1" ]]; then echo "The input dir '$1' does not exist."; exit; fi
outdir="$1"

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




# get reads per group for plotting dataframe
echo "making alignment plot per detected insertion site:" >> "$outdir/log"
>&2 echo "making alignment plot per detected insertion site:"
mkdir -p "$outdir/plots"
while read f; do
	header=$(printf "$f" | awk -F'\t' '{print($1)}')
	# get start position reads
	gapstart=$(printf "$f" | awk -F'\t' '{print($2)}')
	# get end position reads
	gapend=$(printf "$f" | awk -F'\t' '{print($6)}')
	echo "	${header}:	${gapstart} - ${gapend}" >> "$outdir/log"
	>&2 echo "	${header}:	${gapstart} - ${gapend}"
	# PAF target end ($9+1) of read alignment is where gap would start
	awk -F'\t' -v g="$gapstart" '{if(g==$9+1){if($5=="+"){p++; printf("%s\t%s\t%s\t%s\t%s\n","(+)",$8,p,$9,p)}; if($5=="-"){n--; printf("%s\t%s\t%s\t%s\t%s\n","(-)",$8,n,$9,n)}}}' "$outdir/alignments/downselect.paf" > "$outdir/plots/${header}_${gapstart}_${gapend}.data"
	# PAF target start ($8+1) of read alignment is where gap would end
	awk -F'\t' -v g="$gapend" '{if(g==$8+1){if($5=="+"){p++; printf("%s\t%s\t%s\t%s\t%s\n","(+)",$8,p,$9,p)}; if($5=="-"){n--; printf("%s\t%s\t%s\t%s\t%s\n","(-)",$8,n,$9,n)}}}' "$outdir/alignments/downselect.paf" >> "$outdir/plots/${header}_${gapstart}_${gapend}.data"
	# run plot script
	$bin/Rscript $scriptdir/tgif_plot.R "$outdir/plots/${header}_${gapstart}_${gapend}.data"
done < <(tail -n+2 "$outdir/insertions_filtered.tgif")
