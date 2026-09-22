#!/usr/bin/env python3
"""Ring-conformer ensemble docking with Vina-ACO.

For ordinary non-aromatic rings whose prepared puckers remain fixed during docking, this
workflow docks up to K ring conformers (tools/ring_ensemble.py) and merges their poses, ranked by
Vina energy, into one output PDBQT. Ligands without non-aromatic rings are docked once, exactly as
the plain vina_aco command would.

usage: dock_ensemble.py --receptor R.pdbqt (--ligand L.sdf | --members DIR) --config box.txt --out out.pdbqt
                        [--kmax 4] [--num_modes 9] [--split] [--ensemble-seed 42]
                        [--vina bin/vina_aco] -- <extra vina_aco args>
--members uses an ensemble written beforehand by `ring_ensemble.py L.sdf DIR` (for a library, generate the
ensembles once, like any other ligand preparation); the member files are checked against DIR/manifest.json.
All members are docked by one vina_aco process (--batch --ensemble), which reads the receptor and builds the
grid maps once; the poses are byte-identical to docking every member in its own process.
--split divides the preset's nominal evaluation budget among the K members.
"""
import argparse, hashlib, json, re, subprocess, sys, tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

PRESET_BUDGET = {'fast': 300, 'balanced': 600, 'accurate': 3000}


def option_value(args, name, default=None):
    """Return a Vina option written as either `--name value` or `--name=value`."""
    for i, arg in enumerate(args):
        if arg == name:
            if i + 1 >= len(args):
                raise ValueError(f'{name} requires a value')
            return args[i + 1]
        if arg.startswith(name + '='):
            return arg.split('=', 1)[1]
    return default


def models(text):
    """Split a PDBQT with MODEL records into (energy, block) pairs."""
    out, cur = [], []
    for line in text.splitlines(keepends=True):
        if line.startswith('MODEL'):
            cur = []
        elif line.startswith('ENDMDL'):
            block = ''.join(cur)
            e = float(re.search(r'REMARK VINA RESULT:\s+([-+\d.eE]+)', block).group(1))
            out.append((e, block))
        else:
            cur.append(line)
    return out


def pregenerated(d):
    """Member PDBQTs of an ensemble directory, verified against its manifest."""
    d = Path(d)
    man = json.loads((d / 'manifest.json').read_text())
    members = []
    for m in man['members']:
        f = d / m['file']
        if hashlib.sha256(f.read_bytes()).hexdigest() != m['sha256']:
            sys.exit(f'{f} does not match {d / "manifest.json"}')
        members.append(f)
    return members


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--receptor', required=True)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument('--ligand', help='SDF with hydrogens (the ensemble is generated on the fly)')
    src.add_argument('--members', help='directory written by ring_ensemble.py (pre-generated ensemble)')
    ap.add_argument('--config', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--kmax', type=int, default=4); ap.add_argument('--num_modes', type=int, default=9)
    ap.add_argument('--split', action='store_true'); ap.add_argument('--seed', type=int, default=42,
                                                                    help='Vina docking seed (default: 42)')
    ap.add_argument('--ensemble-seed', type=int, default=42, help='ETKDG ensemble seed (default: 42)')
    ap.add_argument('--vina', default=str(HERE.parent / 'bin' / 'vina_aco'))
    ap.add_argument('extra', nargs=argparse.REMAINDER)
    a = ap.parse_args()
    if a.num_modes < 1:
        ap.error('--num_modes must be at least 1')
    extra = a.extra[1:] if a.extra[:1] == ['--'] else a.extra
    try:
        search = option_value(extra, '--search', 'aco')
    except ValueError as e:
        ap.error(str(e))
    if search == 'aco' and option_value(extra, '--search') is None:
        extra = ['--search', 'aco'] + extra
    if a.split and search != 'aco':
        ap.error('--split applies only to --search aco')
    with tempfile.TemporaryDirectory() as tmp:
        if a.members:
            members = pregenerated(a.members)
        else:
            from rdkit import Chem
            from ring_ensemble import ensemble, write_pdbqts
            mol = Chem.MolFromMolFile(a.ligand, removeHs=False)
            if mol is None:
                sys.exit(f'cannot read {a.ligand}')
            ens, info = ensemble(mol, kmax=a.kmax, seed=a.ensemble_seed)
            members = write_pdbqts(ens, tmp, info)
        if not members:
            sys.exit('ligand preparation failed')
        run_extra = list(extra)
        if a.split and len(members) > 1:
            try:
                preset = option_value(run_extra, '--aco_preset', 'balanced')
            except ValueError as e:
                ap.error(str(e))
            if preset not in PRESET_BUDGET:
                ap.error('--aco_preset must be fast, balanced or accurate')
            run_extra += ['--aco_budget', str(PRESET_BUDGET[preset] / len(members))]
        common = [a.vina, '--receptor', a.receptor, '--config', a.config, '--seed', str(a.seed),
                  '--num_modes', str(a.num_modes), '--verbosity', '0']
        if len(members) == 1:
            outs = [Path(tmp) / 'out0.pdbqt']
            cmd = common + ['--ligand', str(members[0]), '--out', str(outs[0])] + run_extra
        else:
            # one process for all members; batch mode names its outputs <input stem>_out.pdbqt
            copies = []
            for k, lig in enumerate(members):
                copies.append(Path(tmp) / f"m{k}.pdbqt")
                copies[-1].write_bytes(Path(lig).read_bytes())
            (Path(tmp) / 'out').mkdir()
            outs = [Path(tmp) / 'out' / f'm{k}_out.pdbqt' for k in range(len(members))]
            cmd = common + ['--batch'] + [str(x) for x in copies] + ['--dir', str(Path(tmp) / 'out'), '--ensemble'] + run_extra
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0 or not all(o.exists() for o in outs):
            sys.exit(p.stderr or p.stdout[-500:])
        pooled = []
        for k, out in enumerate(outs):
            pooled += [(e, k, block) for e, block in models(out.read_text())]
        if not pooled:
            sys.exit('vina_aco returned no parseable poses')
        pooled.sort(key=lambda t: t[0])
        with open(a.out, 'w') as f:
            for i, (e, k, block) in enumerate(pooled[:a.num_modes]):
                f.write(f'MODEL {i + 1}\nREMARK RING_CONFORMER {k} of {len(members)}\n{block}ENDMDL\n')
        print(f'{len(members)} ring conformer(s); best {pooled[0][0]:.3f} kcal/mol from conformer {pooled[0][1]}; wrote {a.out}')


if __name__ == '__main__':
    main()
