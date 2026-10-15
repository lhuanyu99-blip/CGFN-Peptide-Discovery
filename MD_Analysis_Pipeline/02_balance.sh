#!/bin/bash
set -euo pipefail
trap 'echo -e "\n[ERROR] Equilibration failed at line $LINENO"; tail -40 "${LOG_FILE}"; exit 1' ERR

WORK_DIR="$(pwd)"
TOP_FILE="${WORK_DIR}/02_topology/topol.top"
NDX_FILE="${WORK_DIR}/02_topology/final_ndx.ndx"
EM_GRO="${WORK_DIR}/03_simulation/em_steep.gro"
OUTPUT_DIR="${WORK_DIR}/03_simulation"
MDP_DIR="${OUTPUT_DIR}/md_parameters"
LOG_DIR="${WORK_DIR}/logs"
LOG_FILE="${LOG_DIR}/equilibration_3reps.log"

GMX="gmx"
NT_MPI=1
NT_OMP=14
GPU_ID=0

export OMP_PROC_BIND=close
export OMP_PLACES=cores

mkdir -p "${MDP_DIR}" "${LOG_DIR}" >/dev/null 2>&1
> "${LOG_FILE}"

echo "Generating NVT/NPT equilibration parameters..."

# NVT Parameters with Simulated Annealing
cat > "${MDP_DIR}/nvt_classic.mdp" <<EOF
integrator = md
dt = 0.002
nsteps = 500000  ; 1ns
comm-mode = linear
tcoupl = V-rescale
tc-grps = Protein Water_and_ions
tau_t = 1.0 1.0
ref_t = 310 310
gen_vel = yes
gen_temp = 0
gen_seed = -1    ; Randomized velocities for independent replicates
annealing = single single
annealing-npoints = 2 2
annealing-time = 0 1000 0 1000
annealing-temp = 0 310 0 310
pcoupl = no
constraints = h-bonds
constraint-algorithm = LINCS
define = -DPOSRES ; Apply position restraints to heavy atoms
cutoff-scheme = Verlet
nstlist = 20
pbc = xyz
rlist = 1.2
rcoulomb = 1.2
rvdw = 1.2
coulombtype = PME
pme_order = 4
fourierspacing = 0.12
EOF

# NPT Parameters
cat > "${MDP_DIR}/npt_classic.mdp" <<EOF
integrator = md
dt = 0.002
nsteps = 500000  ; 1ns
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
define = -DPOSRES 
cutoff-scheme = Verlet
nstlist = 20
pbc = xyz
rlist = 1.2
rcoulomb = 1.2
rvdw = 1.2
coulombtype = PME
pme_order = 4
fourierspacing = 0.12
EOF

echo "Executing 3 independent replicates for equilibration..."

for REP in 1 2 3; do
    echo "--------------------------------------------------"
    echo "Processing Replicate: ${REP}"
    
    echo "  -> Phase 1: NVT Equilibration (1ns)..."
    ${GMX} grompp -f "${MDP_DIR}/nvt_classic.mdp" -c "${EM_GRO}" -r "${EM_GRO}" -p "${TOP_FILE}" -n "${NDX_FILE}" -o "${OUTPUT_DIR}/nvt_rep${REP}.tpr" -maxwarn 1 >>"${LOG_FILE}" 2>&1
    ${GMX} mdrun -v -deffnm "${OUTPUT_DIR}/nvt_rep${REP}" -ntmpi ${NT_MPI} -ntomp ${NT_OMP} -gpu_id ${GPU_ID} -pin auto >>"${LOG_FILE}" 2>&1
    
    echo "  -> Phase 2: NPT Equilibration (1ns)..."
    ${GMX} grompp -f "${MDP_DIR}/npt_classic.mdp" -c "${OUTPUT_DIR}/nvt_rep${REP}.gro" -r "${OUTPUT_DIR}/nvt_rep${REP}.gro" -t "${OUTPUT_DIR}/nvt_rep${REP}.cpt" -p "${TOP_FILE}" -n "${NDX_FILE}" -o "${OUTPUT_DIR}/npt_rep${REP}.tpr" -maxwarn 1 >>"${LOG_FILE}" 2>&1
    ${GMX} mdrun -v -deffnm "${OUTPUT_DIR}/npt_rep${REP}" -ntmpi ${NT_MPI} -ntomp ${NT_OMP} -gpu_id ${GPU_ID} -pin auto >>"${LOG_FILE}" 2>&1

    cp "${OUTPUT_DIR}/npt_rep${REP}.gro" "${OUTPUT_DIR}/production_input_rep${REP}.gro"
    cp "${OUTPUT_DIR}/npt_rep${REP}.cpt" "${OUTPUT_DIR}/production_input_rep${REP}.cpt"
    echo "  [OK] Replicate ${REP} equilibration completed."
done

echo -e "\n=================================================="
echo "Automated Quality Control (QC) Extraction"
echo "=================================================="
for REP in 1 2 3; do
    echo "[REP ${REP}] Thermodynamic properties:"
    if [ -f "${OUTPUT_DIR}/nvt_rep${REP}.edr" ]; then
        TEMP_RAW=$(echo "Temperature" | ${GMX} energy -f "${OUTPUT_DIR}/nvt_rep${REP}.edr" -o "${OUTPUT_DIR}/temp_rep${REP}.xvg" 2>/dev/null | grep -E "^Temperature")
        echo "  NVT Average Temp: $(echo "${TEMP_RAW}" | awk '{print $2 " ± " $3 " K"}')"
    fi
    if [ -f "${OUTPUT_DIR}/npt_rep${REP}.edr" ]; then
        DENS_RAW=$(echo "Density" | ${GMX} energy -f "${OUTPUT_DIR}/npt_rep${REP}.edr" -o "${OUTPUT_DIR}/density_rep${REP}.xvg" 2>/dev/null | grep -E "^Density")
        echo "  NPT Average Density: $(echo "${DENS_RAW}" | awk '{print $2 " ± " $3 " kg/m³"}')"
    fi
done