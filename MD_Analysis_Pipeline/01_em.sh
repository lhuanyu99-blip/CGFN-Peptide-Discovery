#!/bin/bash
set -uo pipefail
trap 'echo -e "\n[ERROR] Energy minimization failed at line $LINENO"; tail -40 "${LOG_FILE}"; exit 1' ERR

# ==============================================================================
# CONFIGURATION
# ==============================================================================
WORK_DIR="$(pwd)"
TOP_FILE="${WORK_DIR}/02_topology/topol.top"
NDX_FILE="${WORK_DIR}/02_topology/final_ndx.ndx"
INPUT_GRO="${WORK_DIR}/03_simulation/solvated_ions.gro"
OUTPUT_DIR="${WORK_DIR}/03_simulation"
MDP_DIR="${OUTPUT_DIR}/md_parameters"
LOG_DIR="${WORK_DIR}/logs"
LOG_FILE="${LOG_DIR}/em_final.log"

GMX="gmx"           # Change to gmx_mpi if using MPI version
NT_MPI=1            # Number of MPI threads
NT_OMP=14           # Number of OpenMP threads per MPI rank
GPU_ID=0            # Target GPU ID

export OMP_PROC_BIND=close
export OMP_PLACES=cores

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================
log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo -e "${msg}"
    echo "${msg}" >> "${LOG_FILE}"
}

verify_file() { [ -f "$1" ] || { log "[ERROR] Missing required file: $1"; exit 1; }; }

extract_fmax() {
    local log_file=$1
    local fmax_raw=$(grep -E "^[[:space:]]*Maximum force[[:space:]]+=[[:space:]]+[0-9.]+(e[+-]?[0-9]+)?[[:space:]]+on" "${log_file}" | tail -1 | awk '{print $4}')
    
    if [ -z "${fmax_raw}" ] || ! [[ "${fmax_raw}" =~ ^[0-9.]+(e[+-]?[0-9]+)?$ ]]; then
        echo "UNKNOWN"
        return 0
    fi
    
    local fmax=$(echo "${fmax_raw}" | awk '{
        if ($0 ~ /e/) {
            split($0, parts, "e");
            printf "%.0f", parts[1] * (10 ^ parts[2]);
        } else {
            printf "%.0f", $0;
        }
    }')
    echo "${fmax}"
}

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================
mkdir -p "${MDP_DIR}" "${LOG_DIR}" >/dev/null 2>&1
> "${LOG_FILE}"

log "Starting Steepest Descent Energy Minimization"
log "Hardware Allocation: CPU=${NT_MPI}x${NT_OMP}, GPU_ID=${GPU_ID}"

verify_file "${TOP_FILE}"
verify_file "${NDX_FILE}"
verify_file "${INPUT_GRO}"
log "Input files verified successfully."

# Generate MDP File
cat > "${MDP_DIR}/em_steep_standard.mdp" <<EOF
integrator  = steep
emtol       = 1000.0      ; Maximum force convergence threshold
emstep      = 0.01        ; Initial step size (nm)
nsteps      = 50000       ; Maximum number of steps
constraints = none        ; No constraints during EM
cutoff-scheme = Verlet
nstlist     = 20
rlist       = 1.2
rcoulomb    = 1.2
rvdw        = 1.2
coulombtype = PME
pme_order   = 4
fourierspacing = 0.16
EOF

log "Executing Energy Minimization..."
${GMX} grompp -f "${MDP_DIR}/em_steep_standard.mdp" -c "${INPUT_GRO}" -p "${TOP_FILE}" -n "${NDX_FILE}" -o "${OUTPUT_DIR}/em_steep.tpr" -maxwarn 1 >>"${LOG_FILE}" 2>&1
${GMX} mdrun -v -deffnm "${OUTPUT_DIR}/em_steep" -ntmpi ${NT_MPI} -ntomp ${NT_OMP} -gpu_id ${GPU_ID} -pin auto >>"${LOG_FILE}" 2>&1

verify_file "${OUTPUT_DIR}/em_steep.gro"
FINAL_FMAX=$(extract_fmax "${OUTPUT_DIR}/em_steep.log")

log "Energy Minimization Completed."
log "Final Fmax = ${FINAL_FMAX} kJ/mol/nm"
if [ "${FINAL_FMAX}" != "UNKNOWN" ] && [ "${FINAL_FMAX}" -lt 1000 ]; then
    log "Convergence criteria met (Fmax < 1000)."
else
    log "[WARNING] System may not have strictly converged below 1000 kJ/mol/nm, verify structurally before proceeding."
fi