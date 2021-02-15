#!/bin/bash

# display usage info
usage()
{
cat << EOF

Prints converted fastq (in fasta format) to STDOUT.

usage:
fastq2fasta.sh -i <file.fastq>

OPTIONS:
	-h      show this message
	-i	fastq to convert

EOF
}
#===============================================================================
#======				 		 ARGUMENTS & DEFAULTS						  ======
#===============================================================================
# parse input arguments
while getopts "hi:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		i) FASTQ=$OPTARG ;;
		?) usage; exit ;;
	esac
done

if [[ -z "$FASTQ" ]]; then
	usage; exit
fi

#sed '$!N;s/\n/\t/' "$FASTQ" | sed '$!N;s/\n/\t/' | awk -F $'\t' '{sub("^@",">",$1); printf "%s\n%s\n", $1, $2}'

# playera1 added 20190201
sed -n '1~4s/^@/>/p;2~4p' "$FASTQ"
