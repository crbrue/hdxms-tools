from __future__ import annotations
import argparse, subprocess, sys
from importlib.util import find_spec
from pathlib import Path
from .back_exchange import calculate_back_exchange_stats
from .consensus import consensus_from_csv
from .deduplicate import remove_duplicates
from .strip import (
    export_svg_strip,
    export_consensus_full_manual_fixed_size,
    export_consensus_zoom_manual_fixed_size,
)
from .diffplots import plot_diff_bars, plot_diff_map, plot_aligned_diff
from .workflow import run_manifest


def main(argv=None):
    parser=argparse.ArgumentParser(prog='hdxms',description='HDX-MS processing, consensus, and structure coloring')
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('consensus'); p.add_argument('csv'); p.add_argument('-o','--output',required=True); p.add_argument('--value-col',default='pctD_diff')
    p.add_argument('--chain',default=''); p.add_argument('--offset',type=int,default=0); p.add_argument('--offset-per-chain'); p.add_argument('--stat',choices=['mean','median','len_weighted_mean'],default='mean')
    p.add_argument('--mincov',type=int,default=1); p.add_argument('--smooth',type=int,default=0); p.add_argument('--resi-min',type=int); p.add_argument('--resi-max',type=int)
    p=sub.add_parser('strip'); p.add_argument('csv'); p.add_argument('-o','--output',required=True); p.add_argument('--value-col',default='pctD_diff'); p.add_argument('--resi-col',default='resi')
    p.add_argument('--resi-min',type=int); p.add_argument('--resi-max',type=int); p.add_argument('--width-in',type=float); p.add_argument('--height-in',type=float)
    p=sub.add_parser('strip-full'); p.add_argument('csv'); p.add_argument('-o','--output',required=True); p.add_argument('--value-col',default='pctD_diff'); p.add_argument('--resi-col',default='resi'); p.add_argument('--width-in',type=float,default=7.4558); p.add_argument('--height-in',type=float,default=0.0751)
    p=sub.add_parser('strip-zoom'); p.add_argument('csv'); p.add_argument('-o','--output',required=True); p.add_argument('--value-col',default='pctD_diff'); p.add_argument('--resi-col',default='resi'); p.add_argument('--resi-min',type=int,required=True); p.add_argument('--resi-max',type=int,required=True); p.add_argument('--width-in',type=float,default=1.6089); p.add_argument('--height-in',type=float,default=0.0392)
    p=sub.add_parser('deduplicate'); p.add_argument('csv',nargs='+'); p.add_argument('--no-backup',action='store_true'); p.add_argument('--backup-suffix',default='.bak')
    p=sub.add_parser('back-exchange'); p.add_argument('csv'); p.add_argument('--per-output',required=True); p.add_argument('--global-output',required=True)
    p=sub.add_parser('diff-bars'); p.add_argument('csv'); p.add_argument('-o','--output',required=True); p.add_argument('--value-col',default='pctD_diff'); p.add_argument('--label',default='A - B'); p.add_argument('--x-label-every',type=int,default=20); p.add_argument('--no-charge',action='store_true')
    p=sub.add_parser('diff-map'); p.add_argument('csv'); p.add_argument('-o','--output',required=True); p.add_argument('--value-col',default='pctD_diff'); p.add_argument('--label',default='A - B'); p.add_argument('--offset',type=int,default=0); p.add_argument('--resi-min',type=int); p.add_argument('--resi-max',type=int); p.add_argument('--compact',action='store_true')
    p=sub.add_parser('diff-figure'); p.add_argument('csv'); p.add_argument('-o','--output',required=True); p.add_argument('--value-col',default='pctD_diff'); p.add_argument('--offset',type=int,default=0); p.add_argument('--stat',choices=['mean','median','len_weighted_mean'],default='mean'); p.add_argument('--mincov',type=int,default=1); p.add_argument('--smooth',type=int,default=0); p.add_argument('--resi-min',type=int); p.add_argument('--resi-max',type=int); p.add_argument('--consensus-output'); p.add_argument('--title',default='Δ%HDX — consensus and peptide difference map')
    p=sub.add_parser('run'); p.add_argument('manifest')
    p=sub.add_parser('color'); p.add_argument('structure'); p.add_argument('csv'); p.add_argument('--pymol',default='pymol'); p.add_argument('--value-col',default='pctD_diff')
    p.add_argument('--chain-filter'); p.add_argument('--offset',type=int,default=0); p.add_argument('--offset-per-chain'); p.add_argument('--mincov',type=int,default=1)
    p.add_argument('--cartoon-trans',type=float,default=0.0); p.add_argument('--protein-only',action='store_true'); p.add_argument('--out-png'); p.add_argument('--session'); p.add_argument('--out-bfactor-pdb')
    args=parser.parse_args(argv)
    if args.command=='run':
        out = run_manifest(args.manifest)
        print(f'Workflow complete: {out}')
    elif args.command=='consensus':
        result=consensus_from_csv(args.csv,value_col=args.value_col,chain=args.chain,offset=args.offset,offset_per_chain=args.offset_per_chain,stat=args.stat,mincov=args.mincov,smooth=args.smooth,resi_min=args.resi_min,resi_max=args.resi_max,output=args.output)
        print(f'Wrote {args.output} ({len(result)} residues)')
    elif args.command=='strip':
        export_svg_strip(args.csv,args.output,value_col=args.value_col,resi_col=args.resi_col,resi_min=args.resi_min,resi_max=args.resi_max,width_in=args.width_in,height_in=args.height_in); print(f'Wrote {args.output}')
    elif args.command=='strip-full':
        export_consensus_full_manual_fixed_size(args.csv,args.output,width_in=args.width_in,height_in=args.height_in,resi_col=args.resi_col,val_col=args.value_col); print(f'Wrote {args.output}')
    elif args.command=='strip-zoom':
        export_consensus_zoom_manual_fixed_size(args.csv,args.resi_min,args.resi_max,args.output,zoom_width_in=args.width_in,zoom_height_in=args.height_in,resi_col=args.resi_col,val_col=args.value_col); print(f'Wrote {args.output}')
    elif args.command=='deduplicate':
        for path in args.csv:
            before=sum(1 for _ in open(path,encoding='utf-8'))-1
            result=remove_duplicates(path,backup=not args.no_backup,backup_suffix=args.backup_suffix)
            print(f'{path}: {before} -> {len(result)} rows')
    elif args.command=='back-exchange':
        calculate_back_exchange_stats(args.csv,args.per_output,args.global_output); print(f'Wrote {args.per_output} and {args.global_output}')
    elif args.command=='diff-bars':
        plot_diff_bars(args.csv,args.output,value_col=args.value_col,label=args.label,x_label_every=args.x_label_every,show_charge=not args.no_charge); print(f'Wrote {args.output}')
    elif args.command=='diff-map':
        plot_diff_map(args.csv,args.output,value_col=args.value_col,label=args.label,offset=args.offset,resi_min=args.resi_min,resi_max=args.resi_max,compact=args.compact); print(f'Wrote {args.output}')
    elif args.command=='diff-figure':
        plot_aligned_diff(args.csv,args.output,value_col=args.value_col,offset=args.offset,stat=args.stat,mincov=args.mincov,smooth=args.smooth,resi_min=args.resi_min,resi_max=args.resi_max,consensus_output=args.consensus_output,title=args.title); print(f'Wrote {args.output}')
    elif args.command=='color':
        module=Path(__file__).with_name('pymol_script.py')
        cmd=[args.pymol,'-cq',str(module),'--',args.structure,args.csv,'--value-col',args.value_col,'--offset',str(args.offset),'--mincov',str(args.mincov),'--cartoon-trans',str(args.cartoon_trans)]
        for flag,value in [('--chain-filter',args.chain_filter),('--offset-per-chain',args.offset_per_chain),('--out-png',args.out_png),('--session',args.session),('--out-bfactor-pdb',args.out_bfactor_pdb)]:
            if value: cmd += [flag,str(value)]
        if args.protein_only: cmd.append('--protein-only')
        raise SystemExit(subprocess.call(cmd))
