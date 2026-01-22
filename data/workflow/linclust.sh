#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[   -f "$2.dbtype" ] && echo "$2.dbtype exists already!" && exit 1;
[ ! -d "$3" ] && echo "tmp directory $3 not found!" && mkdir -p "$3";

INPUT="$1"
TMP_PATH="$3"
# SOURCE="$INPUT"

# 1. Finding exact $k$-mer matches.
if notExists "${TMP_PATH}/pref.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" kmermatcher "$INPUT" "${TMP_PATH}/pref" ${KMERMATCHER_PAR} \
        || fail "kmermatcher died"
fi
RESULTDB="${TMP_PATH}/pref"

# 2. Ungapped and gapped alignment(w/ banded-block aligner)
if notExists "${TMP_PATH}/aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" align2clust "$INPUT" "$INPUT" "$RESULTDB" "$2" ${ALIGN2CLUST_PAR} \
        || fail "align2clust step died"
fi

if [ -n "$REMOVE_TMP" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY}
    rm -f "${TMP_PATH}/linclust.sh"
fi
