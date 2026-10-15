#!/bin/bash
set -euo pipefail
trap 'echo -e "\n[ERROR] Script failed at line $LINENO"; exit 1' ERR

WORK_DIR="$(pwd)/03_simulation"
GMX="gmx"

cd "${WORK_DIR}"
echo "Initiating automated MM-GBSA calculations (400-500ns window)..."

for REP in 1 2 3; do
    MMPBSA_DIR="${WORK_DIR}/mmpbsa_rep${REP}"
    if [ -f "${MMPBSA_DIR}/FINAL_RESULTS_MMPBSA.dat" ]; then
        echo "  [INFO] Replicate ${REP} already completed. Skipping."
        continue
    fi

    echo "--------------------------------------------------"
    echo "Processing Replicate ${REP}"
    
    TPR="${WORK_DIR}/prod_500ns_rep${REP}.tpr"                 
    FIT_XTC="${WORK_DIR}/prod_500ns_rep${REP}_fit.xtc"         
    NDX="${WORK_DIR}/ana_rep${REP}/ana.ndx"                    
    mkdir -p "${MMPBSA_DIR}"
    
    echo "  -> 1. Isolating protein topology..."
    PROT_TPR="${MMPBSA_DIR}/tpr_protein.tpr"
    echo "1" | ${GMX} convert-tpr -s "${TPR}" -o "${PROT_TPR}" -quiet

    echo "  -> 2. Extracting equilibrium window (400-500ns)..."
    THIN_XTC="${MMPBSA_DIR}/thinned.xtc"
    echo "0" | ${GMX} trjconv -s "${PROT_TPR}" -f "${FIT_XTC}" -o "${THIN_XTC}" -b 400000 -skip 10 -quiet

    cd "${MMPBSA_DIR}"
    echo "  -> 3. Generating gmx_MMPBSA input..."
    cat << EOF > mmpbsa.in
&general
sys_name="Rep${REP}_Complex", 
startframe=1, 
endframe=100, 
forcefields="leaprc.protein.ff14SB",
/
&gb
igb=5, saltcon=0.150,
/
EOF

    echo "  -> 4. Executing gmx_MMPBSA engine..."
    
    gmx_MMPBSA -O \
        -i mmpbsa.in \
        -cs "${PROT_TPR}" \
        -ct "${THIN_XTC}" \
        -ci "${NDX}" \
        -cg 17 18 \
        -nogui > mmpbsa.log 2>&1 || {
            echo "[ERROR] gmx_MMPBSA execution failed. Last 30 lines of log:"
            tail -n 30 mmpbsa.log
            exit 1
        }
        
    echo "  [OK] MM-GBSA calculation completed for Replicate ${REP}."
    cd ..
done