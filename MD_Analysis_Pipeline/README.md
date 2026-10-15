# Automated Molecular Dynamics & Structural Analysis Toolkit

This repository contains the automated bash and Python pipelines utilized to simulate and extract the structural and thermodynamic metrics for the candidate pro-apoptotic peptides (as presented in Figures 6–9 of the manuscript).

## Prerequisites
- **GROMACS** (>= 2021.x)
- **CPPTRAJ** / **gmx_MMPBSA**
- **Python 3.8+** with `MDAnalysis`, `pandas`, `numpy`, `matplotlib`, `seaborn`

## Pipeline Execution Logic
The pipeline is designed to be executed sequentially on a prepared system (containing standard topology and solvated coordinates):

1. **`01_em.sh`**: Performs steepest-descent energy minimization to resolve steric clashes.
2. **`02_balance.sh`**: Executes standard canonical NVT and NPT equilibrations with automated thermodynamic convergence extraction.
3. **`03_production.sh`**: Runs explicit-solvent production MD and automates periodic boundary condition (PBC) corrections and trajectory fitting.
4. **`batch_structure_ultimate.sh`**: The core analytical engine. Extracts RMSD, RMSF, SASA, and utilizes `MDAnalysis` to systematically mine interface interaction fingerprints (Hydrogen bonds & Hydrophobic occupancies). *(Note: User must configure atom index parameters at the top of the script based on their specific complex).*
5. **`batch_mmpbsa_fixed.sh`**: Automates the calculation of binding free energies over the equilibrium trajectory window.
6. **`04_plot_reproduction.py`**: A demonstrative visualization script to automatically generate graphical representations of the extracted `.xvg` and `.csv` matrices.

## Note on Figures
While the final manuscript figures were manually assembled and vector-rendered in professional graphing software (e.g., OriginPro) for publication aesthetics, the raw data structures and identical quantitative trends can be programmatically reproduced using this toolkit.