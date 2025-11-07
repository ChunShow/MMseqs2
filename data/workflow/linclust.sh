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

# 5. Clustering using greedy set cover.
if notExists "${TMP_PATH}/clust.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" alignblock "$INPUT" "$INPUT" "${TMP_PATH}/pref" "${TMP_PATH}/aln" ${ALIGNBLOCK_PAR} \
        || fail "AlignBlock step died"
fi

RESULTDB="${TMP_PATH}/aln"
# 5. Clustering using greedy set cover.
if notExists "${TMP_PATH}/clust.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" clust "$INPUT" "$RESULTDB" "$2" ${CLUSTER_PAR} \
        || fail "Clustering step died"
fi

if [ -n "$REMOVE_TMP" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref_filter1" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref_rescore1" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pre_clust" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input_step_redundancy" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/input_step_redundancy_h" ${VERBOSITY}
    rm -f "${TMP_PATH}/order_redundancy"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref_filter2" ${VERBOSITY}
    if [ -n "$FILTER" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_rescore2" ${VERBOSITY}
    fi
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aln" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/clust" ${VERBOSITY}
    rm -f "${TMP_PATH}/linclust.sh"
fi
