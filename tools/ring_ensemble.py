#!/usr/bin/env python3
"""Ring-conformer ensemble for docking ordinary non-aromatic rings whose prepared puckers remain fixed.

From an input ligand (SDF with explicit hydrogens and 3D coordinates) generate up to K conformers
that differ in the geometry of their non-aromatic rings and write one PDBQT per conformer (Meeko).
The input conformer is always member 0. Ligands without a non-aromatic ring of size >= 5 give only
the input conformer. Only the ligand's own topology and geometry are used (no receptor, no
reference pose).

Selection (deterministic for a fixed seed):
  1. ETKDGv3 embeds `ncand` candidates of the input's stereoisomer (stereo assigned from the input 3D
     geometry; pruneRmsThresh 0.1 A), MMFF94 optimisation; candidates whose MMFF optimisation did not
     converge, or whose 3D stereochemistry differs from the input, are discarded. Without MMFF
     parameters, only the input is used.
  2. Greedy in generation order: keep a candidate if, against every kept member, (a) the RMSD of the
     non-aromatic ring atoms after superposition on those atoms (fixed atom indices, not symmetry-aware)
     exceeds `threshold`, and (b) the symmetry-aware whole heavy-atom RMSD (RDKit GetBestRMS) exceeds
     `threshold` (duplicate guard for symmetric / fused ring systems).
Every emitted PDBQT is read back with Meeko and must reproduce the input's heavy-atom elements, formal
charges and molecular graph (isomeric canonical SMILES; Meeko may renumber atoms, so the comparison is
order-independent), otherwise an error is raised.

usage: ring_ensemble.py input.sdf outdir [--max K] [--ncand N] [--threshold A] [--seed S]
       writes outdir/conf<k>.sdf, outdir/conf<k>.pdbqt and outdir/manifest.json
"""
import argparse, hashlib, json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

VERSION = 1


class EnsembleError(RuntimeError):
    pass


def flexible_ring_atoms(mol):
    atoms = set()
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) >= 5 and not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            atoms.update(ring)
    return sorted(atoms)


def heavy_signature(mol):
    """Order-independent identity: sorted (element, formal charge) of heavy atoms + isomeric canonical
    SMILES (graph incl. stereochemistry). Meeko may renumber atoms, so atom order is not compared."""
    heavy = Chem.RemoveHs(mol)
    if heavy.GetNumConformers():
        Chem.AssignStereochemistryFrom3D(heavy)  # both sides judged by their 3D geometry
    return (sorted((a.GetSymbol(), a.GetFormalCharge()) for a in heavy.GetAtoms()),
            Chem.MolToSmiles(heavy, isomericSmiles=True))


def ensemble(mol, kmax=4, ncand=30, threshold=0.25, seed=42):
    """Returns (molecule whose conformers are the selected ensemble, input first; info dict)."""
    if mol is None:
        raise EnsembleError('input molecule could not be read')
    if kmax < 1:
        raise EnsembleError('kmax must be >= 1')
    if mol.GetNumConformers() != 1 or not mol.GetConformer().Is3D():
        raise EnsembleError('input must contain exactly one 3D conformer')
    if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        raise EnsembleError('input must contain explicit hydrogens')
    ring_atoms = flexible_ring_atoms(mol)
    info = dict(version=VERSION, kmax=kmax, ncand=ncand, threshold=threshold, seed=seed,
                flexible_ring_atoms=len(ring_atoms), candidates=0, converged=0, note='')
    out = Chem.Mol(mol)
    if not ring_atoms or kmax == 1:
        info['note'] = 'no non-aromatic ring' if not ring_atoms else 'kmax 1'
        return out, info
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        info['note'] = 'no MMFF parameters: input conformer only'
        return out, info
    work = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(work)  # embed the input's own stereoisomer
    work.RemoveAllConformers()
    params = AllChem.ETKDGv3(); params.randomSeed = seed; params.pruneRmsThresh = 0.1
    cids = list(AllChem.EmbedMultipleConfs(work, numConfs=ncand, params=params))
    info['candidates'] = len(cids)
    if not cids:
        info['note'] = 'embedding failed: input conformer only'
        return out, info
    status = AllChem.MMFFOptimizeMoleculeConfs(work, maxIters=2000)
    ref_sig = heavy_signature(mol)
    keep = [cid for cid, (flag, _) in zip(cids, status) if flag == 0]
    info['converged'] = len(keep)
    # a candidate must be the same stereoisomer as the input (e.g. cis/trans ring substitution)
    keep = [cid for cid in keep if heavy_signature(Chem.Mol(work, confId=cid)) == ref_sig]
    info['same_stereo'] = len(keep)
    pool = Chem.Mol(mol)
    for cid in keep:
        pool.AddConformer(Chem.Conformer(work.GetConformer(cid)), assignId=True)
    single = lambda cid: Chem.Mol(pool, confId=cid)
    heavy = lambda cid: Chem.RemoveHs(single(cid))
    ids = [c.GetId() for c in pool.GetConformers()]
    chosen = [ids[0]]
    ring_map = [(a, a) for a in ring_atoms]
    for cid in ids[1:]:
        if len(chosen) >= kmax:
            break
        distinct = True
        for j in chosen:
            if rdMolAlign.AlignMol(single(cid), single(j), atomMap=ring_map) <= threshold or \
               rdMolAlign.GetBestRMS(heavy(cid), heavy(j)) <= threshold:
                distinct = False
                break
        if distinct:
            chosen.append(cid)
    out.RemoveAllConformers()
    for cid in chosen:
        out.AddConformer(Chem.Conformer(pool.GetConformer(cid)), assignId=True)
    return out, info


def write_pdbqts(mol, outdir, info=None):
    """Writes one SDF + PDBQT per conformer, validates each PDBQT against the input graph, and writes
    manifest.json. Raises EnsembleError on any preparation or validation failure."""
    from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    ref_sig = heavy_signature(mol)
    members = []
    for k, conf in enumerate(mol.GetConformers()):
        sdf, pdbqt = outdir / f'conf{k}.sdf', outdir / f'conf{k}.pdbqt'
        Chem.MolToMolFile(Chem.Mol(mol, confId=conf.GetId()), str(sdf))
        # Meeko Python API on the SDF round-trip: byte-identical to `mk_prepare_ligand.py -i sdf -o pdbqt`
        # (checked on 60 ensemble members) without a Python start-up per conformer
        try:
            setups = MoleculePreparation().prepare(Chem.MolFromMolFile(str(sdf), removeHs=False))
            text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        except Exception as e:
            raise EnsembleError(f'Meeko failed for conformer {k}: {e}')
        if not ok:
            raise EnsembleError(f'Meeko failed for conformer {k}: {err}')
        pdbqt.write_text(text)
        back = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule.from_file(str(pdbqt), skip_typing=True))
        if len(back) != 1 or back[0] is None or heavy_signature(back[0]) != ref_sig:
            raise EnsembleError(f'conformer {k}: PDBQT does not reproduce the input molecular graph')
        members.append(dict(file=pdbqt.name, sha256=hashlib.sha256(pdbqt.read_bytes()).hexdigest()))
    manifest = dict(info or {}, smiles=ref_sig[1], members=members)
    (outdir / 'manifest.json').write_text(json.dumps(manifest, indent=1))
    return [outdir / m['file'] for m in members]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input'); ap.add_argument('outdir')
    ap.add_argument('--max', type=int, default=4); ap.add_argument('--ncand', type=int, default=30)
    ap.add_argument('--threshold', type=float, default=0.25); ap.add_argument('--seed', type=int, default=42)
    a = ap.parse_args()
    mol = Chem.MolFromMolFile(a.input, removeHs=False)
    ens, info = ensemble(mol, a.max, a.ncand, a.threshold, a.seed)
    for p in write_pdbqts(ens, a.outdir, info):
        print(p)


if __name__ == '__main__':
    main()
