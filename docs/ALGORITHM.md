# Vina-ACO search algorithm

Vina-ACO adds an Ant Colony Optimization (ACO) global search to AutoDock Vina 1.2.x
(`--search aco`). Everything else is unchanged Vina: the scoring functions (vina, vinardo, ad4),
the affinity grids, the BFGS local optimizer, the final refinement with explicit receptor atoms,
output clustering and PDBQT writing. `--search mc` (the default) is the original Monte Carlo search
and produces byte-identical results to stock Vina built with the same compiler flags.

Source: `vina-src/src/lib/aco.{h,cpp}`; integration in `vina.{h,cpp}` (`Vina::set_search`,
`Vina::global_search`) and `src/main/main.cpp` (CLI).

## Search space

A pose is a Vina `conf`: ligand position (3), orientation (unit quaternion), ligand torsions, and
flexible-residue torsions. Every candidate built by an ant is refined with Vina's BFGS (Lamarckian
search: the refined pose replaces the constructed one).

## Default variant: ACO_R (continuous ACO, Socha & Dorigo 2008)

Pheromone is a ranked archive of the `k` best poses found so far (`--aco_archive`, 10).

1. **Initialization** — `k` random poses (uniform position in the box, Haar-uniform orientation,
   uniform torsions), each locally optimized.
2. **Construction** — each ant picks a guide pose from the archive with rank weight
   `w_l ∝ exp(-(l-1)^2 / (2 q^2 k^2))` (`--aco_q`, 0.1), then samples every DOF from a Gaussian kernel
   centred on the guide, with width `σ = ξ · mean distance of that DOF to the other archive members`
   (`--aco_xi`, 0.85). Torsions use periodic distances; orientation uses the rotation angle between
   quaternions and a random rotation-vector kernel; positions are clamped to the box.
   Sampling all DOFs around one guide preserves correlations between position, orientation and
   torsions, which independent per-DOF sampling destroys.
3. **Exploration floor** — with probability `p = explore / (#torsions + 2)` (`--aco_explore`, 1.0,
   i.e. about one DOF per ant) a DOF ignores the archive: a uniform torsion, a uniform random
   orientation, or a position jump of up to 2 Å. This is the continuous counterpart of the MMAS
   lower pheromone bound τ_min: every value always keeps a non-zero construction probability.
4. **Local search** — each ant: `--aco_local_steps` (20) BFGS steps on the capped ("hunt") potential,
   then rescored with the authentic potential (so every archive energy is on the same objective).
   The `--aco_elite` (3) best ants of the iteration get `--aco_elite_steps` (30) more BFGS steps on the
   authentic potential and are offered to the output pose container (Vina's RMSD-diverse list).
5. **Pheromone update** — the archive keeps the best `k` of archive ∪ ants.
6. **Restarts** — after `--aco_patience` (15) iterations without improvement, the archive is
   re-initialized (the output container keeps everything found). This repeats while budget remains.
7. **Stopping** — evaluation budget, or *consensus*: when the best pose of a finished restart phase
   matches the best pose of all earlier phases (within 1 kcal/mol and 2 Å) `--aco_consensus` times.
   Easy ligands therefore stop early; hard ones use the full budget.

## Budget and presets

Evaluations per colony = `aco_budget × (50 + movable atoms + 10 × DOF)`, the same size heuristic
Vina uses for its Monte Carlo step count. `--max_evals N` overrides it.

| preset | aco_budget | consensus | intent |
|---|---|---|---|
| `fast` | 300 | 1 | virtual screening throughput |
| `balanced` (default) | 600 | 2 | general use; held-out Astex mean 0.66 s cold incl. grid setup |
| `accurate` | 3000 | 2 | pose prediction for difficult/flexible ligands |

With `--search aco`, `--exhaustiveness` defaults to 1 (the number of independent colonies); larger
values run several colonies and merge their poses (parallel over `--cpu`).

## Alternative variant: MMAS (`--aco_variant mmas`)

A PLANTS-like MAX-MIN Ant System (Korb, Stützle & Exner 2006/2007): every DOF has a discretized
pheromone trail (0.375 Å position bins; 72 bins for torsions and for the Shoemake orientation
parameters `u1,u2,u3`, which map uniform trails to Haar-uniform rotations). Ants sample each DOF
independently (roulette wheel + uniform jitter inside the bin), are locally optimized, and the global
best plus ranked iteration-best ants deposit pheromone (with ±1-bin smoothing). Trails evaporate
(`--aco_evap`) and are clamped to `[τ_max/(tau_ratio · bins), τ_max]`. It is kept for comparison.
At an equal 150k-evaluation budget on the leniently prepared PoseBusters validation split
(118-complex development split, 1 seed), its top score came within 0.5 kcal/mol of stock Vina
exhaustiveness 8 in 66% of runs vs 75% for ACO_R and 56% for Vina MC, consistent with independent
per-DOF sampling breaking the correlations between DOFs. Top-1 pose success at that budget was
similar for the two ACO variants (37.3% MMAS, 39.0% ACO_R; overlapping intervals) and 21.2% for MC.

## Speed work (all bit-identical to stock Vina arithmetic)

- Intramolecular pair tables are shared per XS type pair (vina/vinardo) instead of per atom pair
  (AD4 keeps per-atom tables because of charges): the ~26 MB of per-pair spline tables for a
  40-atom ligand shrink to a few hundred KB, which removes setup time and cache misses.
- Grid construction (`cache::populate`) uses per-cell contiguous receptor atom arrays and
  pre-resolved lookup tables, with the same summation order.
- The receptor–ligand type-pair tables are only built for the ligand's atom types in single-ligand
  mode.
Net effect for 1iep (single-case measurements): `--score_only` 0.73 s → 0.27 s, and an MC exh 1
docking ~18% faster with byte-identical output.

## References

- K. Socha, M. Dorigo. Ant colony optimization for continuous domains. EJOR 185 (2008) 1155–1173.
  doi:10.1016/j.ejor.2006.06.046 — our ACO_R is an adaptation (SO(3) kernel, periodic torsions,
  Lamarckian BFGS, exploration floor, restarts, consensus stop), not the exact published algorithm.
- O. Korb, T. Stützle, T. E. Exner. An ant colony optimization approach to flexible protein–ligand
  docking. Swarm Intelligence 1 (2007) 115–134. doi:10.1007/s11721-007-0006-9
- O. Korb, T. Stützle, T. E. Exner. PLANTS: Application of ant colony optimization to structure-based
  drug design. ANTS 2006, LNCS 4150, 247–258.
- O. Trott, A. J. Olson. AutoDock Vina. J Comput Chem 31 (2010) 455–461; J. Eberhardt et al.,
  AutoDock Vina 1.2.0. J Chem Inf Model 61 (2021) 3891–3898.
- K. Shoemake. Uniform random rotations. Graphics Gems III (1992) 124–132.
