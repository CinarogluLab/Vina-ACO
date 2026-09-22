# Vina-ACO

Vina-ACO is AutoDock Vina 1.2.7 with an ant colony optimization (ACO) global search. Vina's scoring
functions, grid maps, BFGS local optimizer, pose clustering and PDBQT output are unchanged; the global
search can be switched between Vina's Monte Carlo and two ACO variants:

- `--search aco` (continuous ACO_R, the default ACO variant): a ranked solution archive over
  translations, torsions and orientations, with Gaussian kernels for the Euclidean coordinates and an
  SO(3)-correct kernel for the orientation.
- `--search aco --aco_variant mmas`: a PLANTS-like MAX-MIN Ant System with discretized pheromone trails,
  kept for comparison.

Without `--search aco` the program runs Vina's Monte Carlo search and gives the same poses as upstream
Vina 1.2.7 built with the same compiler and seed.

The repository also contains `tools/ring_ensemble.py` and `tools/dock_ensemble.py`, which dock up to four
ring conformers of a ligand and pool the poses. In standard Meeko preparation, 5- and 6-membered
non-aromatic rings are held at their input pucker during docking; the ensemble gives the search other
puckers to choose from.

## What to expect

On PoseBusters and Astex Diverse (501 complexes, three seeds; details in the paper):

| search | top-1 success (RMSD ≤ 2 Å) | time per ligand (one core) |
|---|---:|---:|
| Vina-ACO `fast` | 36.7% | 0.63 s |
| Vina-ACO `balanced` | 41.7% | 0.83 s |
| Vina-ACO `accurate` | 47.0% | 2.19 s |
| Vina-ACO `accurate` + ring ensemble | 52.7% | 8.55 s |
| Vina MC, exhaustiveness 1 | 38.7% | 4.81 s |
| Vina MC, exhaustiveness 8 | 50.4% | 37.9 s |
| Vina MC, exhaustiveness 8 + ring ensemble | 56.2% | 111.1 s |

Times are for one ligand in a fresh process, including grid construction (and, for the ensembles,
conformer generation and docking every member).
The `balanced` preset uses about nine times fewer scoring evaluations than Vina at exhaustiveness 1 for
similar pose accuracy. In batch mode (grids built once per receptor), `fast` and `balanced` docked
6.7 and 3.9 ligands/s on one core in a 24-ligand test. These numbers come from one machine
(AMD Ryzen 7 7745HX) and will differ on yours.

## Build

Requirements: a C++14 compiler, CMake ≥ 3.16, and Boost (program_options, filesystem, system, thread,
serialization).

```bash
./build.sh                      # -> bin/vina_aco
./bin/vina_aco --version
```

`build.sh` uses `$BOOST_ROOT` if set, otherwise `$CONDA_PREFIX`, otherwise the Boost that CMake finds.
A conda environment with everything needed is in `environment.yml`:

```bash
conda env create -f environment.yml
conda activate vina-aco
./build.sh
```

The build targets the current CPU (`-march=native`). For a portable binary use
`./build.sh -DVINA_ACO_NATIVE=OFF`.

## Usage

`vina_aco` accepts every option of the `vina` command line; see the
[AutoDock Vina documentation](https://autodock-vina.readthedocs.io). The ACO options are:

| option | meaning |
|---|---|
| `--search aco` | use the ant colony instead of Monte Carlo |
| `--aco_preset fast\|balanced\|accurate` | parameter set (default `balanced`) |
| `--aco_variant acor\|mmas` | continuous ACO_R (default) or MAX-MIN Ant System |
| `--exhaustiveness N` | with ACO: number of independent colonies (default 1) |
| `--max_evals N` | evaluation budget per colony (overrides the preset) |

Expert options (`--aco_budget`, archive size, ants, kernel width and others) are listed by
`./bin/vina_aco --help`. An explicit expert option overrides the preset.

| preset | intended use |
|---|---|
| `fast` | high-throughput screening |
| `balanced` | general use |
| `accurate` | pose prediction, flexible ligands |

Single ligand:

```bash
./bin/vina_aco --receptor receptor.pdbqt --ligand ligand.pdbqt --config box.txt \
               --search aco --aco_preset balanced --cpu 1 --seed 1 --out ligand_out.pdbqt
```

Screening a library (the receptor is read and the grids are built once):

```bash
./bin/vina_aco --receptor receptor.pdbqt --batch ligands/*.pdbqt --config box.txt \
               --search aco --aco_preset fast --cpu 1 --dir poses/
```

### Ring-conformer ensembles

The tools need Python 3 with RDKit and Meeko (both in `environment.yml`). The input ligand is an SDF
with explicit hydrogens and 3D coordinates.

```bash
# one step: build the ensemble and dock it
python tools/dock_ensemble.py --receptor receptor.pdbqt --ligand ligand.sdf --config box.txt \
       --out ligand_out.pdbqt -- --search aco --aco_preset accurate --cpu 1

# for a library: build each ensemble once, as part of ligand preparation ...
python tools/ring_ensemble.py ligand.sdf ens/ligand      # conf<k>.pdbqt + manifest.json
# ... then dock the stored members (one vina_aco process, grids built once)
python tools/dock_ensemble.py --receptor receptor.pdbqt --members ens/ligand --config box.txt \
       --out ligand_out.pdbqt -- --search aco --aco_preset accurate --cpu 1
```

Member 0 is always the input conformer, and a ligand without a non-aromatic ring is docked once, exactly
as plain `vina_aco` would dock it. Poses are ranked by Vina score and tagged `REMARK RING_CONFORMER k`.
Meeko opens rings of 7 or more atoms (macrocycles) with glue atoms, so Vina already samples those; the
ensemble matters most for 5- and 6-membered rings.

## Example

`example/aco_docking/` holds one PoseBusters complex (5S8I) with its box, the prepared receptor and
ligand, a crystal reference pose, and a stored two-member ring ensemble:

```bash
cd example/aco_docking
../../bin/vina_aco --receptor receptor.pdbqt --ligand ligand.pdbqt --config box.txt \
                   --search aco --aco_preset balanced --cpu 1 --seed 1 --out out.pdbqt
python ../../tools/dock_ensemble.py --receptor receptor.pdbqt --members ensemble --config box.txt \
       --out ens_out.pdbqt --vina ../../bin/vina_aco -- --search aco --aco_preset accurate --cpu 1
```

`example/basic_docking/` and `example/flexible_docking/` are the upstream Vina examples (without the
precomputed AutoDock4 maps); add `--search aco` to any of their `vina` commands.

## Tests

```bash
tests/run_math_test.sh                 # ACO internals: archive weights, SO(3) sampling, kernel scaling
python tests/test_ring_ensemble.py     # ensemble generator: determinism, stereo, member 0 = direct prep
python tests/test_ensemble_batch.py    # --batch --ensemble reproduces separate runs byte for byte
```

## Documentation

[docs/ALGORITHM.md](docs/ALGORITHM.md) describes the search: the solution archive, the sampling kernels,
the budget rule, the presets, the MMAS variant and the speed changes.

## Citation

If you use Vina-ACO, please cite the paper (in review; the reference will be added here) and AutoDock Vina:

- J. Eberhardt, D. Santos-Martins, A. F. Tillack, S. Forli. AutoDock Vina 1.2.0: New docking methods,
  expanded force field, and Python bindings. *J. Chem. Inf. Model.* 2021, 61, 3891–3898.
- O. Trott, A. J. Olson. AutoDock Vina: improving the speed and accuracy of docking with a new scoring
  function, efficient optimization, and multithreading. *J. Comput. Chem.* 2010, 31, 455–461.

`CITATION.cff` has the machine-readable form.

## License

Apache License 2.0, as for AutoDock Vina. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
