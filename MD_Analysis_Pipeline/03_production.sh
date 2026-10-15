#!/bin/bash
set -euo pipefail
trap 'echo -e "\n[ERROR] Production simulation interrupted at line $LINENO"; exit 1' ERR

WORK_DIR="$(pwd)"
TOP_FILE="${WORK_DIR}/02_topology/topol.top"
NDX_FILE="${WORK_DIR}/02_topology/final_ndx.ndx"
OUTPUT_DIR="${WORK_DIR}/03_simulation"
MDP_DIR="${OUTPUT_DIR}/md_parameters"
LOG_DIR="${WORK_DIR}/logs"

GMX="gmx"
NT_MPI=1
NT_OMP=14
GPU_ID=0
export OMP_PROC_BIND=close
export OMP_PLACES=cores

mkdir -p "${LOG_DIR}"

echo "Generating 500ns Production MDP parameters..."
cat > "${MDP_DIR}/production_classic.mdp" <<EOF
integrator = md
dt = 0.002
nsteps = 250000000 ; 500 ns
comm-mode = linear
tcoupl = V-rescale
tc-grps = Protein Water_and_ions
tau_t = 1.0 1.0
ref_t = 310 310
pcoupl = C-rescale
pcoupltype = isotropic
tau_p = 2.0
ref_p = 1.0
compressibility = 4.5e-5
constraints = h-bonds
constraint-algorithm = LINCS
cutoff-scheme = Verlet
nstlist = 100
pbc = xyz
rlist = 1.2
rcoulomb = 1.2
rvdw = 1.2
coulombtype = PME
pme_order = 4
fourierspacing = 0.12
nstxout-compressed = 50000     ; Save every 100 ps
compressed-x-grps  = Protein   ; Save only protein coordinates to reduce storage
EOF

for REP in 1 2 3; do
    echo "--------------------------------------------------"
    echo "Starting 500ns Production: Replicate ${REP}"
    
    GRO_FILE="${OUTPUT_DIR}/production_input_rep${REP}.gro"
    CPT_FILE="${OUTPUT_DIR}/production_input_rep${REP}.cpt"
    SIM_BASE="${OUTPUT_DIR}/prod_500ns_rep${REP}"
    SIM_LOG="${LOG_DIR}/prod_500ns_rep${REP}.log"
    
    ${GMX} grompp -f "${MDP_DIR}/production_classic.mdp" -c "${GRO_FILE}" -t "${CPT_FILE}" \
               -p "${TOP_FILE}" -n "${NDX_FILE}" -o "${SIM_BASE}.tpr" -maxwarn 1
               
    RUN_CMD="${GMX} mdrun -v -deffnm ${SIM_BASE} -ntmpi ${NT_MPI} -ntomp ${NT_OMP} -gpu_id ${GPU_ID} -pin auto -g ${SIM_LOG}"
    
    # Auto-resume logic if interrupted
    if [ -f "${SIM_BASE}.cpt" ]; then
        echo "[INFO] Checkpoint detected. Resuming simulation in append mode..."
        RUN_CMD="${RUN_CMD} -cpi ${SIM_BASE}.cpt -append"
    fi

    $RUN_CMD
    echo "[OK] Replicate ${REP} production completed."
done

echo -e "\n=================================================="
echo "Initiating Trajectory Post-Processing (PBC Correction)"
echo "=================================================="

for REP in 1 2 3; do
    SIM_BASE="${OUTPUT_DIR}/prod_500ns_rep${REP}"
    if [ -f "${SIM_BASE}.xtc" ]; then
        # Center trajectory (Group 1: Protein, Group 1: Protein)
        echo "1 1" | ${GMX} trjconv -s "${SIM_BASE}.tpr" -f "${SIM_BASE}.xtc" -n "${NDX_FILE}" -o "${SIM_BASE}_center.xtc" -pbc mol -center -quiet
        
        # Fit rotational and translational motions (Group 4: Backbone, Group 1: Protein)
        echo "4 1" | ${GMX} trjconv -s "${SIM_BASE}.tpr" -f "${SIM_BASE}_center.xtc" -n "${NDX_FILE}" -o "${SIM_BASE}_fit.xtc" -fit rot+trans -quiet
        
        echo "  [OK] Rep ${REP} trajectory centered and fitted."
    fi
done