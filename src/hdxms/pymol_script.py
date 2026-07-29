from __future__ import annotations
import argparse, os
from collections import defaultdict
from pymol import cmd, util
from hdxms.colors import PYMOL_RGB, bin_name
from hdxms.io import load_hdx_csv, apply_offsets


def _ranges(values):
    values = sorted(set(values))
    if not values: return []
    out=[]; start=end=values[0]
    for value in values[1:]:
        if value == end + 1: end=value
        else: out.append(f"{start}-{end}" if start != end else str(start)); start=end=value
    out.append(f"{start}-{end}" if start != end else str(start))
    return out


def run(args):
    obj = os.path.splitext(os.path.basename(args.structure))[0]
    cmd.load(args.structure, obj)
    for name, rgb in PYMOL_RGB.items(): cmd.set_color(name, rgb)
    df = apply_offsets(load_hdx_csv(args.csv, args.value_col), args.offset, args.offset_per_chain)
    if args.chain_filter:
        allowed=set(args.chain_filter.split(',')); df=df[(df.chain=='') | df.chain.isin(allowed)]
    sums=defaultdict(float); counts=defaultdict(int)
    for row in df.itertuples(index=False):
        for residue in range(row.start, row.end+1): sums[(row.chain,residue)] += row.value; counts[(row.chain,residue)] += 1
    model=cmd.get_model(f"{obj} and polymer.protein")
    present=defaultdict(set)
    for atom in model.atom:
        try: present[atom.chain].add(int(atom.resi))
        except ValueError: pass
    mapped={}
    for chain,residues in present.items():
        if args.chain_filter and chain not in set(args.chain_filter.split(',')): continue
        for residue in residues:
            key=(chain,residue) if (chain,residue) in sums else ('',residue)
            if counts.get(key,0) >= args.mincov: mapped[(chain,residue)] = sums[key]/counts[key]
    cmd.hide('everything','all'); cmd.show('cartoon',f'{obj} and polymer.protein')
    cmd.color('hdx_unmapped',f'{obj} and polymer.protein')
    grouped=defaultdict(lambda: defaultdict(list))
    for (chain,residue),value in mapped.items(): grouped[bin_name(value)][chain].append(residue)
    for color,chain_map in grouped.items():
        if color is None: continue
        for chain,residues in chain_map.items():
            cmd.color(color, f"{obj} and chain {chain} and polymer.protein and resi {'+'.join(_ranges(residues))}")
    if not args.protein_only:
        cmd.show('sticks',f'{obj} and polymer.nucleic'); util.cbaw(f'{obj} and polymer.nucleic')
    cmd.set('cartoon_transparency',args.cartoon_trans); cmd.bg_color('white'); cmd.orient(f'{obj} and polymer.protein')
    if args.out_bfactor_pdb:
        cmd.alter(f'{obj} and polymer.protein','b=0.0')
        for (chain,residue),value in mapped.items(): cmd.alter(f'{obj} and chain {chain} and resi {residue}',f'b={float(value)}')
        cmd.save(args.out_bfactor_pdb,obj)
    if args.session: cmd.save(args.session)
    if args.out_png:
        cmd.viewport(args.width,args.height); cmd.ray(); cmd.png(args.out_png,dpi=args.dpi)


def main(argv=None):
    parser=argparse.ArgumentParser()
    parser.add_argument('structure'); parser.add_argument('csv')
    parser.add_argument('--value-col',default='pctD_diff'); parser.add_argument('--chain-filter')
    parser.add_argument('--offset',type=int,default=0); parser.add_argument('--offset-per-chain')
    parser.add_argument('--mincov',type=int,default=1); parser.add_argument('--cartoon-trans',type=float,default=0.0)
    parser.add_argument('--protein-only',action='store_true'); parser.add_argument('--out-png'); parser.add_argument('--session')
    parser.add_argument('--out-bfactor-pdb'); parser.add_argument('--width',type=int,default=1800); parser.add_argument('--height',type=int,default=1400)
    parser.add_argument('--dpi',type=int,default=300)
    run(parser.parse_args(argv))


if __name__ == "__main__":
    main()
