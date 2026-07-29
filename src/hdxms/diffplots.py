from __future__ import annotations

from pathlib import Path
from typing import Iterable
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

from .colors import HDX_BINS, hex_color
from .consensus import consensus_from_csv

PARP1_DOMAINS = [
    ("Zn1", 7, 97), ("Zn2", 102, 201), ("Zn3", 214, 372),
    ("BRCT", 373, 486), ("WGR", 523, 656),
    ("HD", 678, 787), ("ART", 788, 1014),
]


def _load(path: str | Path, value_col: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"Start", "End", value_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    for c in ["Start", "End", "Charge", value_col]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    if "Charge" not in df.columns:
        df["Charge"] = 1
    return df.dropna(subset=["Start", "End", value_col]).copy()


def plot_diff_bars(csv, output, value_col="pctD_diff", label="A - B", x_label_every=20,
                   show_charge=True, dpi=300, figsize=(14, 4.5)):
    df = _load(csv, value_col).sort_values(["Start", "End", "Charge"]).reset_index(drop=True)
    x = np.arange(len(df))
    if show_charge:
        labels = [f"{int(s)}-{int(e)}(z{int(z)})" for s,e,z in zip(df.Start,df.End,df.Charge)]
    else:
        labels = [f"{int(s)}-{int(e)}" for s,e in zip(df.Start,df.End)]
    shown = [lab if i % max(1, x_label_every) == 0 else "" for i,lab in enumerate(labels)] if len(labels)>x_label_every else labels
    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(x, df[value_col].to_numpy(), width=0.35, color="black")
    ax.axhline(0, linewidth=1, color="gray")
    ax.set_xticks(x, shown, rotation=90, fontsize=8)
    ax.set_ylabel(f"%D difference ({label})")
    ax.set_xlabel("Peptide (Start-End" + (", z)" if show_charge else ")"))
    ax.set_title("Peptide-charge %D difference", weight="bold")
    fig.tight_layout()
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return df


def _legend_handles():
    return [Patch(facecolor=h, edgecolor="none", label=label) for _,h,label in [
        (0,"#282a73","< -30"),(0,"#4478bb","-30 to -15"),(0,"#bde6f5","-15 to -5"),
        (0,"#afb0b0","-5 to +5"),(0,"#fac9c9","+5 to +15"),(0,"#ed2024","+15 to +30"),(0,"#7f1416","> +30")]]


def _assign_lanes(rows: list[dict], gap=0):
    ends=[]; lanes=[]
    for r in rows:
        for i,end in enumerate(ends):
            if r["start"] > end + gap:
                ends[i]=r["end"]; lanes.append(i); break
        else:
            ends.append(r["end"]); lanes.append(len(ends)-1)
    return lanes, len(ends)


def _draw_diffmap(ax, df, value_col, offset=0, xmin=None, xmax=None, compact=True):
    rows=[]
    for r in df.itertuples(index=False):
        s=int(getattr(r,"Start"))+offset; e=int(getattr(r,"End"))+offset
        if e<s: s,e=e,s
        if xmin is not None and e<xmin: continue
        if xmax is not None and s>xmax: continue
        if xmin is not None: s=max(s,xmin)
        if xmax is not None: e=min(e,xmax)
        rows.append({"start":s,"end":e,"value":float(getattr(r,value_col))})
    rows.sort(key=lambda x:(x["start"],-(x["end"]-x["start"])))
    if compact:
        lanes,nlanes=_assign_lanes(rows)
    else:
        lanes=list(range(len(rows))); nlanes=len(rows)
    for row,lane in zip(rows,lanes):
        ax.add_patch(Rectangle((row["start"]-.5,lane),row["end"]-row["start"]+1,.8,
                               linewidth=0,facecolor=hex_color(row["value"])))
    ax.set_ylim(-.2,max(1,nlanes))
    ax.set_yticks([])
    return rows,nlanes


def plot_diff_map(csv, output, value_col="pctD_diff", offset=0, resi_min=None, resi_max=None,
                  compact=False, label="A - B", dpi=300, figsize=(12,6)):
    df=_load(csv,value_col)
    fig,ax=plt.subplots(figsize=figsize)
    rows,_=_draw_diffmap(ax,df,value_col,offset,resi_min,resi_max,compact)
    if not rows: raise ValueError("No peptides overlap the requested residue range")
    ax.set_xlim((resi_min if resi_min is not None else min(r['start'] for r in rows))-.5,
                (resi_max if resi_max is not None else max(r['end'] for r in rows))+.5)
    ax.set_xlabel("Residue index")
    ax.set_title(f"Peptide uptake difference map: {label}")
    ax.legend(handles=_legend_handles(),title="Δ%D",loc="upper right")
    fig.tight_layout(); Path(output).parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=dpi,bbox_inches="tight",facecolor="white"); plt.close(fig)
    return df


def plot_aligned_diff(csv, output, value_col="pctD_diff", offset=0, stat="mean", mincov=1,
                      smooth=0, domains: Iterable[tuple[str,int,int]]|None=None,
                      resi_min=None, resi_max=None, consensus_output=None,
                      title="Δ%HDX — consensus and peptide difference map", dpi=300,
                      figsize=(16,8.5), tick_step=50):
    df=_load(csv,value_col)
    cons=consensus_from_csv(csv,value_col=value_col,offset=offset,stat=stat,mincov=mincov,
                            smooth=smooth,resi_min=resi_min,resi_max=resi_max,
                            output=consensus_output)
    valid=cons.dropna(subset=[value_col])
    if valid.empty: raise ValueError("No consensus residues available")
    xmin=int(resi_min if resi_min is not None else valid.resi.min())
    xmax=int(resi_max if resi_max is not None else valid.resi.max())
    fig=plt.figure(figsize=figsize)
    gs=fig.add_gridspec(3,1,height_ratios=[1.0,.8,5.0],hspace=.12)
    ax1=fig.add_subplot(gs[0]); ax2=fig.add_subplot(gs[1],sharex=ax1); ax3=fig.add_subplot(gs[2],sharex=ax1)
    for r in cons.itertuples(index=False):
        if xmin <= int(r.resi) <= xmax:
            ax1.add_patch(Rectangle((int(r.resi)-.5,0),1,1,linewidth=0,facecolor=hex_color(getattr(r,value_col))))
    ax1.set_xlim(xmin-.5,xmax+.5); ax1.set_ylim(0,1); ax1.set_yticks([]); ax1.tick_params(labelbottom=False); ax1.set_title(title)
    for name,s,e in (domains or PARP1_DOMAINS):
        s+=offset; e+=offset
        a=max(s,xmin); b=min(e,xmax)
        if b>=a:
            ax2.add_patch(Rectangle((a-.5,.2),b-a+1,.6,facecolor="white",edgecolor="black",linewidth=.8))
            ax2.text((a+b)/2,.5,name,ha="center",va="center",fontsize=8)
    ax2.set_ylim(0,1); ax2.set_yticks([]); ax2.tick_params(labelbottom=False)
    _draw_diffmap(ax3,df,value_col,offset,xmin,xmax,compact=True)
    ticks=list(range(((xmin+tick_step-1)//tick_step)*tick_step,xmax+1,tick_step))
    ax3.set_xticks(ticks); ax3.set_xlabel("Residue index"); ax3.legend(handles=_legend_handles(),title="Δ%D",loc="upper right")
    fig.tight_layout(); Path(output).parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=dpi,bbox_inches="tight",facecolor="white"); plt.close(fig)
    return cons
