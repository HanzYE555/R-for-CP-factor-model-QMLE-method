"""Describe the direction-dependent rho denominator; no model/search reruns."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/lambda_axis_geometry_20261030"
SOURCE = OUT / "primary_by_seed.csv"
PROJECT = ROOT / "outputs/Id_Prob_QMLE_CP_Model_JASA_S23"
KEY = ["p", "q", "regime", "seed", "T"]
CELL = ["p", "q", "regime", "T"]

d = pd.read_csv(SOURCE)
assert len(d) == 10796 and not d.duplicated(KEY).any()
lam = d[[f"rho_lambda_{j:02d}" for j in range(1,9)]].fillna(0).to_numpy(float)
i,j = np.triu_indices(8,k=1)
d["rho_lambda_norm"] = np.linalg.norm(lam,axis=1)
d["rho_k_norm"] = np.linalg.norm(lam[:,i]*lam[:,j],axis=1)
d["rho_H"] = d.rho_min*d.rho_k_norm
assert np.isfinite(d[["rho_H","rho_k_norm"]]).all().all()
assert (d.E1_at_rho > 0).all()
d["rho_H_over_E1"] = d.rho_H/d.E1_at_rho
r = d.E1_at_rho.to_numpy()**2
h2 = d.rho_H.to_numpy()**2
low = r*(1-r)
high = r*(1-r/2)
assert (h2 >= low-2e-10).all() and (h2 <= high+2e-10).all()
metrics = ["rho_min","E1_at_rho","rho_lambda_norm","rho_k_norm","rho_H","rho_H_over_E1","rho_max_loading_cosine"]
byseed = d[KEY+metrics]
summaries = []
for cell, g in d.groupby(CELL,sort=True):
    row = dict(zip(CELL,cell));row["n"] = len(g)
    for name in metrics:
        row[name+"_mean"] = g[name].mean()
        row[name+"_variance"] = g[name].var(ddof=1)
    summaries.append(row)
s = pd.DataFrame(summaries)
assert len(s) == 108 and s.n.sum() == 10796
vdir = OUT / "validation"
vdir.mkdir(exist_ok=True)
byseed.to_csv(vdir / "rho_denominator_by_seed.csv",index=False)
s.to_csv(vdir / "rho_denominator_mean_variance.csv",index=False)
s.to_csv(OUT / "tables/rho_denominator_mean_variance.csv",index=False)
small = s[(s.p==40)&(s.q==40)].copy()
small["_order"] = small.regime.map({"interior":0,"boundary":1,"exterior":2})
small = small.sort_values(["_order","T"]).drop(columns="_order")
small.to_csv(vdir / "rho_denominator_40x40_mean_variance.csv",index=False)
small.to_csv(OUT / "tables/rho_denominator_40x40_mean_variance.csv",index=False)

def cell(row, name, digits=6):
    m,v = row[name+"_mean"],row[name+"_variance"]
    return r"\shortstack{"+f"{m:.{digits}f}"+r"\\[-1pt]{\scriptsize["+f"{v:.2e}"+r"]}}"

tex = [r"% New coefficient-angle diagnostic; computed from saved unit-M lambda vectors.",
       r"\begin{table}[p]",r"\begingroup",r"\spacingset{1.0}",r"\small\centering",
       r"\caption{Direction-dependent scaling of the polynomial ratio on $\mathcal D^\lambda_{\delta_0}$ for $(p,q)=(40,40)$. Entries are mean [sample variance], with 100 seeds in every row. All quantities are evaluated at the lowest found $\rho_{\mathrm{ind}}$ candidate after imposing $\|m_r(\lambda)\|_F=1$.}",
       r"\label{s23:lambda:tab:denominator}",r"\footnotesize",r"\setlength{\tabcolsep}{4pt}",
       r"\renewcommand{\arraystretch}{1.10}",r"\begin{adjustbox}{max width=\linewidth}",
       r"\begin{tabular}{lrcccc}",r"\toprule",
       r"Design & $T$ & $\rho_{\mathrm{ind},\min}$ & $E_{1,\mathrm{center}}$ & $\|k(\lambda)\|_2$ & $H(m_r(\lambda))/E_{1,\mathrm{center}}$\\",
       r"\midrule"]
labels = {"interior":"Identifiable interior","boundary":"Identifiable boundary","exterior":"Non-identifiable"}
for regime,g in small.groupby("regime",sort=False):
    for k,(_,row) in enumerate(g.iterrows()):
        tex.append((labels[regime] if k==0 else "")+f" & {int(row['T'])} & "+" & ".join(
            cell(row,name,4 if name=="rho_k_norm" else 6) for name in ["rho_min","E1_at_rho","rho_k_norm","rho_H_over_E1"])+r"\\")
    if regime != "exterior":tex.append(r"\addlinespace")
tex += [r"\bottomrule",r"\end{tabular}",r"\end{adjustbox}",r"\par\smallskip",
        r"\begin{minipage}{0.98\linewidth}\footnotesize",
        r"Here $k(\lambda)=(\lambda_i\lambda_j)_{i<j}$, $H$ is the Euclidean norm of the distinct $2\times2$ minors, and $H=\rho_{\mathrm{ind}}\|k(\lambda)\|_2$. The ratio $H/E_1$ is computed within each seed before taking its mean. The matrix has unit Frobenius norm; the coefficient vector need not have unit Euclidean norm. No additional exclusion or residual threshold is applied. The much larger exterior denominator explains why a small polynomial ratio can coexist with a visibly larger rank-one error at the same direction. It does not by itself establish an additional exact rank-one direction. The full 108-cell mean/variance table and all per-seed denominators are retained in the accompanying data.",
        r"\end{minipage}",r"\endgroup",r"\end{table}",r"\clearpage"]
target = PROJECT / "additions/lambda/denominator_table.tex"
target.write_text("\n".join(tex)+"\n",encoding="utf-8")
audit = dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
             fits=len(d),cells=len(s),display_rows=len(small),normalization="matrix Frobenius norm 1",
             denominator="norm((lambda_i*lambda_j)_{i<j})",extra_thresholds_applied=False,
             maximum_lower_bound_violation=float(np.maximum(0,low-h2).max()),
             maximum_upper_bound_violation=float(np.maximum(0,h2-high).max()),
             H_over_E1_min=float(d.rho_H_over_E1.min()),H_over_E1_max=float(d.rho_H_over_E1.max()))
(vdir / "rho_denominator_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
print(small[["regime","T","rho_min_mean","E1_at_rho_mean","rho_k_norm_mean","rho_H_over_E1_mean"]].to_string(index=False))
print(json.dumps(audit,indent=2))
