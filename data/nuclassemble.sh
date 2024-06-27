#!/bin/sh -e

# Assembler workflow script
fail() {
    echo "Error: $1"
    exit 1
}

deleteIncremental() {
    if [ -n "$REMOVE_INCREMENTAL_TMP" ] &&  [ -n "$1" ]; then
         "$MMSEQS" rmdb "$1"
    fi
}

notExists() {
	[ ! -f "$1" ]
}

cyclecheck() {
	if [ -n "$CALL_CYCLE_CHECK" ]; then
        if notExists "${1}_cycle.done"; then
            # shellcheck disable=SC2086
            "$MMSEQS" cyclecheck "$1" "${1}_cycle" ${CYCLE_CHECK_PAR} \
                || fail "Cycle check step died"


            if [ -s "${1}_cycle" ]; then

                if notExists "${1}_noneCycle"; then
                    awk 'NR==FNR { a[$1]=$0; next } !($1 in a) {print $0}' "${1}_cycle.index" \
                    "${1}.index" > "${1}_noneCycle.index"
                    ln -s "$1" "${1}_noneCycle"
                    ln -s "${1}.dbtype" "${1}_noneCycle.dbtype"
                fi

                if [ -z "$PREV_CYCLE_ALL" ]; then
                    # shellcheck disable=SC2086
                    "$MMSEQS" mvdb "${1}_cycle" "${1}_cycle_all"
                else
                    # shellcheck disable=SC2086
                    "$MMSEQS" concatdbs "${PREV_CYCLE_ALL}" "${1}_cycle" "${1}_cycle_all" --preserve-keys
                fi

            else
                ln -s "$1" "${1}_noneCycle"
                ln -s "${1}.index" "${1}_noneCycle.index"
                ln -s "${1}.dbtype" "${1}_noneCycle.dbtype"
            fi
            touch "${1}_cycle.done"
            deleteIncremental "$PREV_CYCLE"
            PREV_CYCLE="${1}_cycle"
        fi

        if [ -s "${1}_cycle_all" ]; then
            deleteIncremental "${PREV_CYCLE_ALL}"
            PREV_CYCLE_ALL="${1}_cycle_all"
        fi

        PREV_CONTIG_ASSEMBLY="${1}_noneCycle"
    fi
}


# check input variables
[ -z "${OUT_FILE}" ] && echo "Please provide OUT_FILE" && exit 1
[ -z "${TMP_PATH}" ] && echo "Please provide TMP_PATH" && exit 1

# check if files exist
[   -f "${OUT_FILE}" ] &&  echo "${OUT_FILE} exists already!" && exit 1
[ ! -d "${TMP_PATH}" ] &&  echo "tmp directory ${TMP_PATH} not found!" && mkdir -p "${TMP_PATH}"

if [ -n "${DB_MODE}" ]; then
    INPUT="$1"
    [ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
else
    INPUT="${TMP_PATH}/nucl_reads"
    if notExists "${TMP_PATH}/nucl_reads"; then
         if [ -n "${PAIRED_END}" ]; then
            echo "PAIRED END MODE"
            # shellcheck disable=SC2086
            "$MMSEQS" mergereads "$@" "${TMP_PATH}/nucl_reads" ${VERBOSITY_PAR} \
              || fail "mergereads failed"
         else
            # shellcheck disable=SC2086
            "$MMSEQS" createdb "$@" "${TMP_PATH}/nucl_reads" ${CREATEDB_PAR} \
                || fail "createdb failed"
         fi
    fi
fi

SOURCE=${INPUT}
STEP=0
if [ -z "$NUM_IT" ]; then
    NUM_IT=1;
fi

while [ $STEP -lt $NUM_IT ]; do
    echo "STEP: $STEP"

    if [ $STEP -lt $NUM_IT_READS ]; then

        # 1. Finding exact $k$-mer matches.
        if notExists "${TMP_PATH}/pref_${STEP}.done"; then
            # shellcheck disable=SC2086
            "$MMSEQS" kmermatcher "$INPUT" "${TMP_PATH}/pref_${STEP}" ${KMERMATCHER_READS_PAR} \
                || fail "Kmer matching step died"
            deleteIncremental "$PREV_KMER_PREF"
            touch "${TMP_PATH}/pref_${STEP}.done"
            PREV_KMER_PREF="${TMP_PATH}/pref_${STEP}"
        fi

        # 2. Ungapped alignment
        if notExists "${TMP_PATH}/aln_${STEP}.done"; then
            # shellcheck disable=SC2086
            "$MMSEQS" rescorediagonal "$INPUT" "$INPUT" "${TMP_PATH}/pref_${STEP}" "${TMP_PATH}/aln_${STEP}" ${UNGAPPED_ALN_PAR_ANCIENT_READS} \
                || fail "Ungapped alignment step died"
            touch "${TMP_PATH}/aln_${STEP}.done"
            deleteIncremental "$PREV_ALN"
            PREV_ALN="${TMP_PATH}/aln_${STEP}"
        fi

        # Louis was here
        # 2.5 Deamination correction
        if notExists "${TMP_PATH}/correction.done"; then
            # shellcheck disable=SC2086
            "$MMSEQS" ancient_correction "$INPUT" "${TMP_PATH}/aln_${STEP}" "${TMP_PATH}/correction_${STEP}" ${ASSEMBLE_RESULT_PAR_ANCIENT_READS} \
                || fail "ancient_correction died"
            touch "${TMP_PATH}/correction_${STEP}.done"
            deleteIncremental "$PREV_CORR"
            PREV_CORR="${TMP_PATH}/correction_${STEP}"
        fi

        # 3. Assemble only with reads
        if notExists "${TMP_PATH}/assembly_reads_${STEP}.done"; then
            # shellcheck disable=SC2086
            "$MMSEQS" ancient_read_assemble "${TMP_PATH}/correction_${STEP}" "${TMP_PATH}/aln_${STEP}" "${TMP_PATH}/assembly_reads_${STEP}" ${ASSEMBLE_RESULT_PAR_ANCIENT_READS} \
                || fail "ancient_read_assemble step died"
            touch "${TMP_PATH}/assembly_reads_${STEP}.done"
            deleteIncremental "$PREV_ASSEMBLY_READS"
            deleteIncremental "$PREV_ASSEMBLY_READS_STEP"
            PREV_ASSEMBLY_READS="${TMP_PATH}/assembly_reads_${STEP}"
            PREV_ASSEMBLY_READS_STEP="${TMP_PATH}/assembly_reads_${STEP}"
        fi

        INPUT="${PREV_ASSEMBLY_READS}"
        STEP="$((STEP+1))"

    else

        # 1. Finding exact $k$-mer matches.
        if notExists "${TMP_PATH}/pref_asm_${STEP}.done"; then
            # shellcheck disable=SC2086
            "$MMSEQS" kmermatcher "$INPUT" "${TMP_PATH}/pref_asm_${STEP}" ${KMERMATCHER_CONTIGS_PAR} \
            || fail "Kmer matching step died"
            deleteIncremental "$PREV_KMER_ASM_PREF"
            touch "${TMP_PATH}/pref_asm_${STEP}.done"
            PREV_KMER_ASM_PREF="${TMP_PATH}/pref_asm_${STEP}"
        fi

        # 5. Ungapped alignment
        if notExists "${TMP_PATH}/aln_asm_${STEP}.done"; then
            # shellcheck disable=SC2086
            "$MMSEQS" rescorediagonal "$INPUT" "$INPUT" "${TMP_PATH}/pref_asm_${STEP}" "${TMP_PATH}/aln_asm_${STEP}" ${UNGAPPED_ALN_PAR_ANCIENT_CONTIGS} \
                || fail "Ungapped alignment step died"
            touch "${TMP_PATH}/aln_asm_${STEP}.done"
            deleteIncremental "$PREV_ALN_ASM"
            PREV_ALN_ASM="${TMP_PATH}/aln_asm_${STEP}"
        fi

        # Louis was here
        # 5.5 Deamination correction
        if notExists "${TMP_PATH}/correction.done"; then
            # shellcheck disable=SC2086
            "$MMSEQS" ancient_correction "$INPUT" "${TMP_PATH}/aln_asm_${STEP}" "${TMP_PATH}/correction_${STEP}" ${ASSEMBLE_RESULT_PAR_ANCIENT_CONTIGS} \
                || fail "ancient_correction died"
            touch "${TMP_PATH}/correction_${STEP}.done"
            deleteIncremental "$PREV_CORR"
            PREV_CORR="${TMP_PATH}/correction_${STEP}"
        fi

        # 6. Assemble only with contigs
        if notExists "${TMP_PATH}/assembly_contigs_${STEP}.done"; then
            # shellcheck disable=SC2086
            "$MMSEQS" ancient_contig_merge "${TMP_PATH}/correction_${STEP}" "${TMP_PATH}/aln_asm_${STEP}" "${TMP_PATH}/assembly_contigs_${STEP}" ${ASSEMBLE_RESULT_PAR_ANCIENT_CONTIGS} \
                || fail "ancient_contig_merge step died"
            touch "${TMP_PATH}/assembly_contigs_${STEP}.done"
            deleteIncremental "$PREV_CONTIG_ASSEMBLY"
            deleteIncremental "$PREV_CONTIG_ASSEMBLY_STEP"
        fi

        PREV_CONTIG_ASSEMBLY="${TMP_PATH}/assembly_contigs_${STEP}"
        PREV_CONTIG_ASSEMBLY_STEP="${TMP_PATH}/assembly_contigs_${STEP}"
        cyclecheck "${PREV_CONTIG_ASSEMBLY}"

        INPUT="${PREV_CONTIG_ASSEMBLY}"
        STEP="$((STEP+1))"

    fi

    done
    STEP="$((STEP-1))"
    RESULT="${TMP_PATH}/assembly_reads_${STEP}"

if [ -n "$PREV_CYCLE_ALL" ]; then

    RESULT="${TMP_PATH}/assembly_merged"
    if notExists "${TMP_PATH}/assembly_merged"; then
        # shellcheck disable=SC2086
        "$MMSEQS" concatdbs "${PREV_CONTIG_ASSEMBLY}" "${PREV_CYCLE_ALL}" "${TMP_PATH}/assembly_merged" --preserve-keys \
             || fail "Concatenation of non cyclic and cyclic contigs died"
    fi
fi

# select only assembled sequences
if notExists "${RESULT}_only_assembled.index"; then
    awk 'NR == FNR { f[$1] = $0; next } $1 in f { print f[$1], $0 }' "${RESULT}.index" "${SOURCE}.index" > "${RESULT}_tmp.index"
    awk '$3 > $7 { print }' "${RESULT}_tmp.index" > "${RESULT}_only_assembled.index"
fi

# select only sequences fullfilling a minimum length threshold
if notExists "${RESULT}_only_assembled_filtered.index"; then
    # shellcheck disable=SC208
    awk -v thr="${MIN_CONTIG_LEN}" '$3 > (thr+1) { print }' "${RESULT}_only_assembled.index" > "${RESULT}_only_assembled_filtered.index"
fi

# create db outfile
if notExists "${OUT_FILE}.dbtype"; then
    "$MMSEQS" createsubdb "${RESULT}_only_assembled_filtered.index" "${RESULT}" "${TMP_PATH}/assembly" --subdb-mode 0 \
        || fail "Create filtered contig db died"
    if [ -n "$PREV_CYCLE_ALL" ]; then
        awk 'NR == FNR { f[$1] = $0; next } $1 in f { print $0 }' "${PREV_CYCLE_ALL}.index" "${TMP_PATH}/assembly.index" > "${TMP_PATH}/assembly_cycle.index"
    fi
fi

if [ -z "${DB_MODE}" ]; then

    if notExists "${TMP_PATH}/assembly_h.dbtype"; then
        # shellcheck disable=SC2086
        if [ -f "${TMP_PATH}/assembly_cycle.index" ]; then
            "$MMSEQS" createhdb "${TMP_PATH}/assembly" "${TMP_PATH}/assembly_cycle" "${TMP_PATH}/assembly" ${VERBOSITY_PAR} \
                || fail "createhdb failed"
        else
            "$MMSEQS" createhdb "${TMP_PATH}/assembly" "${TMP_PATH}/assembly" ${VERBOSITY_PAR} \
                || fail "createhdb failed"
        fi
    fi

    if notExists "${TMP_PATH}/assembly.fasta"; then
        # shellcheck disable=SC2086
        "$MMSEQS" convert2fasta "${TMP_PATH}/assembly" "${TMP_PATH}/assembly.fasta" ${VERBOSITY_PAR} \
            || fail "convert2fasta died"
    fi

    mv -f "${TMP_PATH}/assembly.fasta" "$OUT_FILE" \
        || fail "Could not move result to $OUT_FILE"

else
    "$MMSEQS" mvdb "${TMP_PATH}/assembly" "$OUT_FILE"\
       || fail "Could not move result to $OUT_FILE"

    if [ -f "${TMP_PATH}/assembly_cycle.index" ]; then
        mv "${TMP_PATH}/assembly_cycle.index" "${OUT_FILE}_cycle.index"
    fi
fi

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    rm -f "${TMP_PATH}/pref_"*
    rm -f "${TMP_PATH}/aln_"*
    rm -f "${TMP_PATH}/assembly"*
    rm -f "${TMP_PATH}/nuclassemble.sh"
fi








