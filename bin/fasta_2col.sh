#!/bin/bash

sed $'$!N;s/\\\n/\t/' "$@"


# alternative
#awk '{printf /^>/ ? "\n"$0"\t" : $0}' "$@" | tail -n+2
#printf "\n" ""
