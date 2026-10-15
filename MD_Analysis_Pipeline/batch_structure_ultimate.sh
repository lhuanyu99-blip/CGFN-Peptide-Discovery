#!/bin/bash
set -euo pipefail
trap 'echo -e "\n[ERROR] Analysis failed at line $LINENO"; exit 1' ERR

WORK_DIR="$(pwd)/03_simulation"
GMX="gmx"

# ==============================================================================
# USER CONFIGURATION (Adapt these atom indices to your specific complex)
# ==============================================================================
RECEPTOR_ATOMS=2977          # Total atoms in the receptor
TOTAL_PROTEIN_ATOMS=3430     # Total atoms (receptor + peptide ligand)
CORE_RESIDUES="54-127"       # Residue index for specific domain analysis (e.g., Core)
LATCH_RESIDUES="130-165"     # Residue index for specific domain analysis (e.g., Latch)
# ==============================================================================

cd "${WORK_DIR}"
echo "Initiating Automated Structural & Interaction Analysis Pipeline"

for REP in 1 2 3; do
    echo "--------------------------------------------------"
    echo "Analyzing Replicate: ${REP}"
    TPR="${WORK_DIR}/prod_500ns_rep${REP}.tpr"
    FIT_XTC="${WORK_DIR}/prod_500ns_rep${REP}_fit.xtc"
    ANA_DIR="${WORK_DIR}/ana_rep${REP}"
    mkdir -p "${ANA_DIR}"
    cd "${ANA_DIR}"
    
    echo "  -> 1. Constructing dynamic index groups..."
    PEP_START_ATOM=$((RECEPTOR_ATOMS + 1))
    
    ${GMX} make_ndx -f "${TPR}" -o ana.ndx -quiet <<EOF
a 1-${RECEPTOR_ATOMS}
a ${PEP_START_ATOM}-${TOTAL_PROTEIN_ATOMS}
r ${CORE_RESIDUES}
r ${LATCH_RESIDUES}
q
EOF

    echo "  -> 2. Extracting basic dynamics (RMSD, RMSF, DSSP, Gyration)..."
    # Complex RMSD (Group 1: Protein)
    echo "1 1" | ${GMX} rms -s "${TPR}" -f "${FIT_XTC}" -n ana.ndx -o rmsd_complex.xvg -quiet
    
    # Receptor RMSF
    echo "17" | ${GMX} rmsf -s "${TPR}" -f "${FIT_XTC}" -n ana.ndx -o rmsf_receptor.xvg -res -quiet
    echo "17" | ${GMX} gyrate -s "${TPR}" -f "${FIT_XTC}" -n ana.ndx -o gyrate.xvg -quiet

    echo "  -> 3. Extracting Conformational specific metrics (SASA)..."
    ${GMX} sasa -s "${TPR}" -f "${FIT_XTC}" -surface "atomnr 1 to ${RECEPTOR_ATOMS}" -o sasa_receptor.xvg -quiet

    echo "  -> 4. Generating stripped TPR for Python interaction analysis..."
    DRY_TPR="dry_protein_rep${REP}.tpr"
    echo "1" | ${GMX} convert-tpr -s "${TPR}" -o "${DRY_TPR}" -quiet

    echo "  -> 5. Executing MDAnalysis interaction profiler (Hydrogen Bonds & Hydrophobics)..."
    cat << EOF_PY > global_interaction_profiler.py
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.lib.distances import distance_array
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

u = mda.Universe("${DRY_TPR}", "${FIT_XTC}")
total_frames = u.trajectory.n_frames

rec_sel = "bynum 1:${RECEPTOR_ATOMS}"
lig_sel = "protein and not bynum 1:${RECEPTOR_ATOMS}"

print("       Extracting Hydrogen Bond Networks...")
hbonds = HydrogenBondAnalysis(u, between=[rec_sel, lig_sel], d_a_cutoff=3.5, d_h_a_angle_cutoff=150)
hbonds.run(verbose=False)
df_hb = pd.DataFrame(hbonds.results.hbonds, columns=['frame', 'd_idx', 'h_idx', 'a_idx', 'dist', 'angle'])

hb_over_time = df_hb.groupby('frame').size().reset_index(name='count')
all_frames = pd.DataFrame({'frame': range(total_frames)})
hb_over_time = pd.merge(all_frames, hb_over_time, on='frame', how='left').fillna(0)
hb_over_time.to_csv("hbnum_total.csv", index=False)

d_atoms = u.atoms[df_hb['d_idx'].values.astype(int)]
a_atoms = u.atoms[df_hb['a_idx'].values.astype(int)]
df_hb['pair'] = d_atoms.resnames + d_atoms.resids.astype(str) + " - " + a_atoms.resnames + a_atoms.resids.astype(str)

hb_occ = (df_hb['pair'].value_counts() / total_frames * 100).reset_index()
hb_occ.columns = ['Interaction Pair', 'Occupancy(%)']
hb_occ[hb_occ['Occupancy(%)'] > 10.0].to_csv("top_hbonds.csv", index=False)

print("       Extracting Hydrophobic Contact Networks...")
hydro_sel = "resname ALA VAL ILE LEU MET PHE TYR TRP and not name C CA N O"
rec_hydro = u.select_atoms(f"{rec_sel} and {hydro_sel}")
lig_hydro = u.select_atoms(f"{lig_sel} and {hydro_sel}")
contact_counts = {}
for ts in u.trajectory:
    dist = distance_array(rec_hydro.positions, lig_hydro.positions)
    idx_r, idx_l = np.where(dist < 4.0)
    frame_pairs = set()
    for r, l in zip(idx_r, idx_l):
        res_r = rec_hydro[r].resname + str(rec_hydro[r].resid)
        res_l = lig_hydro[l].resname + str(lig_hydro[l].resid)
        frame_pairs.add(f"{res_r} - {res_l}")
    for pair in frame_pairs:
        contact_counts[pair] = contact_counts.get(pair, 0) + 1

hydro_data = [{"Interaction Pair": pair, "Occupancy(%)": (count/total_frames)*100} for pair, count in contact_counts.items() if (count/total_frames)*100 > 10.0]
df_hydro = pd.DataFrame(hydro_data)
if not df_hydro.empty:
    df_hydro = df_hydro.sort_values(by="Occupancy(%)", ascending=False)
df_hydro.to_csv("top_hydrophobic.csv", index=False)
print("       [OK] Interaction profiling completed.")
EOF_PY
    python global_interaction_profiler.py
    cd ..
done
echo "Global structural extraction completed."