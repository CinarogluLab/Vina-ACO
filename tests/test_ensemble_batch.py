#!/usr/bin/env python3
"""Verify that one-process ensemble docking reproduces separate member processes byte for byte."""
import argparse
import json
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--binary', type=Path, default=ROOT / 'bin' / 'vina_aco')
    ap.add_argument('--core', default='0')
    a = ap.parse_args()
    binary = a.binary.resolve()
    data = ROOT / 'example' / 'aco_docking'
    ens = data / 'ensemble'
    members = [ens / m['file'] for m in json.loads((ens / 'manifest.json').read_text())['members']]
    if len(members) < 2:
        raise SystemExit('test fixture needs at least two ensemble members')
    common = ['--receptor', str(data / 'receptor.pdbqt'), '--config', str(data / 'box.txt'),
              '--cpu', '1', '--seed', '19', '--verbosity', '0', '--search', 'aco',
              '--aco_preset', 'fast', '--max_evals', '10000']
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        separate = []
        for i, member in enumerate(members):
            out = tmp / f'separate{i}.pdbqt'
            subprocess.run(['taskset', '-c', a.core, str(binary), '--ligand', str(member),
                            '--out', str(out)] + common, check=True, capture_output=True)
            separate.append(out.read_bytes())
        batch_dir = tmp / 'batch'
        batch_dir.mkdir()
        subprocess.run(['taskset', '-c', a.core, str(binary), '--batch', *map(str, members),
                        '--dir', str(batch_dir), '--ensemble'] + common, check=True, capture_output=True)
        batch = [(batch_dir / f'{member.stem}_out.pdbqt').read_bytes() for member in members]
        if batch != separate:
            raise SystemExit('batch ensemble poses differ from separate member processes')

        pooled = tmp / 'pooled.pdbqt'
        subprocess.run(['taskset', '-c', a.core, 'python3', str(ROOT / 'tools' / 'dock_ensemble.py'),
                        '--receptor', str(data / 'receptor.pdbqt'), '--members', str(ens),
                        '--config', str(data / 'box.txt'), '--out', str(pooled), '--seed', '19',
                        '--vina', str(binary), '--', '--search', 'aco', '--aco_preset', 'fast',
                        '--max_evals', '10000'], check=True, capture_output=True)
        text = pooled.read_text()
        if 'REMARK RING_CONFORMER' not in text or 'REMARK VINA RESULT' not in text:
            raise SystemExit('wrapper did not produce a parseable pooled pose file')
    print(f'PASS  ensemble batch byte equivalence ({len(members)} members) and pooled wrapper output')


if __name__ == '__main__':
    main()
