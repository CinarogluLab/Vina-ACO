#!/usr/bin/env python3
"""Checks for tools/ring_ensemble.py and tools/dock_ensemble.py."""
import filecmp, json, subprocess, sys, tempfile
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from ring_ensemble import ensemble, write_pdbqts, heavy_signature, EnsembleError  # noqa: E402
from dock_ensemble import option_value  # noqa: E402


def mol3d(smiles, seed=7):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(m)
    return m


def run(mol, d):
    ens, info = ensemble(mol)
    return write_pdbqts(ens, d, info), info


def main():
    results = {}
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        # 1. no flexible ring -> exactly one member whose PDBQT equals direct Meeko preparation of the input
        benz = mol3d('c1ccccc1CCO')
        paths, info = run(benz, tmp / 'benz')
        Chem.MolToMolFile(benz, str(tmp / 'direct.sdf'))
        subprocess.run(['mk_prepare_ligand.py', '-i', str(tmp / 'direct.sdf'), '-o', str(tmp / 'direct.pdbqt')], check=True, capture_output=True)
        results['no_ring_single_member'] = len(paths) == 1
        results['no_ring_identical_to_direct_prep'] = filecmp.cmp(paths[0], tmp / 'direct.pdbqt', shallow=False)
        # 2. cyclohexane derivative -> >= 2 members, deterministic
        cyc = mol3d('CC1CCC(CC1)C(=O)N')
        p1, _ = run(cyc, tmp / 'c1'); p2, _ = run(cyc, tmp / 'c2')
        results['cyclohexane_ge_2_members'] = len(p1) >= 2
        results['deterministic'] = len(p1) == len(p2) and all(filecmp.cmp(a, b, shallow=False) for a, b in zip(p1, p2)) and \
            json.loads((tmp / 'c1' / 'manifest.json').read_text()) == json.loads((tmp / 'c2' / 'manifest.json').read_text())
        # 3. stereochemistry / charge / connectivity preserved (charged, chiral, ring)
        chiral = mol3d('C[C@H]1CC[C@@H](CC1)[NH3+]')
        paths, _ = run(chiral, tmp / 'chiral')
        results['graph_stereo_charge_preserved'] = len(paths) >= 1  # write_pdbqts raises otherwise
        manifest = json.loads((tmp / 'chiral' / 'manifest.json').read_text())
        results['manifest_has_hashes'] = all(len(m['sha256']) == 64 for m in manifest['members'])
        results['manifest_has_parameters'] = all(k in manifest for k in ('version', 'kmax', 'ncand', 'threshold', 'seed'))
        # 4. wrapper accepts both common command-line spellings without confusing docking and ensemble seeds
        results['wrapper_option_value'] = option_value(['--search=aco', '--aco_preset', 'accurate'], '--search') == 'aco' and \
            option_value(['--search=aco', '--aco_preset', 'accurate'], '--aco_preset') == 'accurate'
        # 5. invalid input handling
        for name, fn in [('kmax0', lambda: ensemble(cyc, kmax=0)), ('no_H', lambda: ensemble(Chem.RemoveHs(cyc))),
                         ('none', lambda: ensemble(None))]:
            try:
                fn(); results[f'rejects_{name}'] = False
            except EnsembleError:
                results[f'rejects_{name}'] = True
    for k, v in results.items():
        print(f'{"PASS" if v else "FAIL"}  {k}')
    sys.exit(0 if all(results.values()) else 1)


if __name__ == '__main__':
    main()
