## Meshes

This folder contains four versions of the computational mesh with different resolutions. The **baseline** mesh is the default one used in the main case (already placed in `constant/polyMesh`). The other meshes are provided for convergence studies.

| Mesh | Number of Cells | Description |
|------|-----------------|------------|
| Baseline | 250000 | Default mesh (used in the main case) |
| Coarse | 108900 | Coarse mesh for quick tests |
| Medium | 168100 |  Medium resolution |
| Fine | 372100 | Fine mesh for convergence analysis 

Each mesh is stored in its own subdirectory with the full OpenFOAM structure:
*   `baseline/polyMesh`
*   `coarse/polyMesh`
*   `medium/polyMesh`
*   `fine/polyMesh`

To run the case with a different mesh, copy the desired mesh to `constant/polyMesh/` before running `./Allrun`.
