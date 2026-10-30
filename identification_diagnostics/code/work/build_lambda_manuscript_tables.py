"""Generate primary lambda-angle-domain LaTeX tables from new geometry exports.

Never reads archived matrix-domain geometric scores.  The only reused numerical
source is the complete J spectrum, which depends on the saved QMLE fit and is
unchanged by the off-axis search domain.  All generated labels have the prefix
``s23:lambda:``.  This script never edits supplement.tex or existing additions.

Usage:
    python work/build_lambda_manuscript_tables.py --data-root outputs/lambda_axis_geometry_20261030
The default requires all 108 cells / 10,796 available fits; --wait-seconds can
wait for aggregation outputs without weakening the completeness checks.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import time
from itertools import product
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "outputs/Id_Prob_QMLE_CP_Model_JASA_S23"
TIMES = [100, 200, 400, 800]
DIMS = [20, 40, 80]
KEYS = ["p", "q", "regime", "T"]
REGIMES = {"interior": "Identifiable interior", "boundary": "Identifiable boundary",
           "exterior": "Non-identifiable"}
LETTERS = {"interior": "I", "boundary": "B", "exterior": "E"}
PAIRS = [("exterior", "interior"), ("exterior", "boundary"), ("boundary", "interior")]
REQUIRED = ["key_metrics_mean_var.csv", "geometry_summary.csv", "cap10_mean_variance.csv",
            "search_quality.csv", "geometry_auc.csv", "separation_auc_comparison_wide.csv"]
SPECTRA_DEFAULT = ROOT / "work/grid_final_spectra/all_J_singular_values_mean_variance.csv"
LABEL = "s23:lambda:"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def number(value):
    return pd.notna(value) and math.isfinite(float(value))


class Tables:
    def __init__(self, data_root, output, spectra):
        self.data_root, self.output = Path(data_root), Path(output)
        tables = self.data_root / "tables"
        self.paths = {name: tables / name for name in REQUIRED}
        self.paths["spectra"] = Path(spectra)
        self.key = pd.read_csv(self.paths[REQUIRED[0]])
        self.geom = pd.read_csv(self.paths[REQUIRED[1]])
        if "subset" in self.geom:
            self.geom = self.geom[self.geom["subset"].eq("all_valid")].copy()
        self.cap = pd.read_csv(self.paths[REQUIRED[2]])
        self.quality = pd.read_csv(self.paths[REQUIRED[3]])
        self.auc = pd.read_csv(self.paths[REQUIRED[4]])
        if "subset" in self.auc:
            self.auc = self.auc[self.auc["subset"].eq("all_valid")].copy()
        self.cap_auc = pd.read_csv(self.paths[REQUIRED[5]])
        self.J = pd.read_csv(spectra)
        self.mv_records, self.count_records, self.auc_records, self.tables = [], [], [], []
        self.outputs = {}
        expected = set(product(DIMS, DIMS, REGIMES, TIMES))
        for name, df in [("key", self.key), ("geometry", self.geom),
                         ("cap", self.cap), ("quality", self.quality)]:
            assert len(df) == 108, (name, len(df))
            assert not df.duplicated(KEYS).any(), name
            assert set(df[KEYS].itertuples(index=False, name=None)) == expected, name
        assert int(self.key.n.sum()) == 10796
        for z in self.key.itertuples():
            expected_n = 96 if (z.p, z.q, z.regime, z.T) == (80, 80, "interior", 100) else 100
            assert z.n == expected_n
        assert len(self.J) == 3456
        assert not self.J.duplicated(KEYS + ["ell", "svd_index"]).any()
        for name in ["cap10_axis_edge_count", "cap10_physical_cosine_gt090_count",
                     "cap10_max_loading_cosine_mean", "E1_max_loading_cosine_mean"]:
            assert name in self.key, f"Missing required new lambda-domain auxiliary: {name}"
        assert "cap10_loading_edge_count" not in self.key, "Ambiguous legacy physical constraint alias"

    def pick(self, name, p, q, reg, T):
        df = getattr(self, name)
        z = df[(df.p == p) & (df.q == q) & (df.regime == reg) & (df["T"] == T)]
        assert len(z) == 1, (name, p, q, reg, T, len(z))
        return z.iloc[0]

    def mv(self, mean, variance, context, decimals=6, allow_missing=False):
        if not number(mean):
            assert allow_missing, context
            self.mv_records.append(dict(context, source_mean=None, source_variance=None,
                                        display=None, missing=True))
            return r"\textemdash"
        mean = float(mean)
        if not number(variance):
            assert allow_missing, context
            out = r"\shortstack{" + f"{mean:.{decimals}f}" + r"\\[-1pt]{\scriptsize[\text{NA}]}}"
            self.mv_records.append(dict(context, source_mean=mean, source_variance=None,
                                        display=out, missing_variance=True))
            return out
        variance = float(variance)
        assert variance >= 0, (context, variance)
        m, v = f"{mean:.{decimals}f}", f"{variance:.2e}"
        out = r"\shortstack{" + m + r"\\[-1pt]{\scriptsize[" + v + "]}}"
        self.mv_records.append(dict(context, source_mean=mean, source_variance=variance,
                                    display=out, display_mean=m, display_variance=v))
        return out

    def cell(self, z, metric, context, decimals=6, conditional=False):
        return self.mv(z[metric + "_mean"], z[metric + "_variance"],
                       dict(context, metric=metric), decimals, conditional)

    @staticmethod
    def row(values):
        return " & ".join(map(str, values)) + r"\\" + "\n"

    def start(self, caption, label, clear=True, tabcolsep="4pt", arraystretch="1.03"):
        assert label.startswith(LABEL)
        self.tables.append(dict(label=label, caption=caption))
        return "\n".join(([r"\clearpage"] if clear else []) + [
            r"\begin{table}[p]", r"\begingroup", r"\spacingset{1.0}",
            r"\small\centering", r"\caption{" + caption + "}",
            r"\label{" + label + "}", r"\footnotesize",
            r"\setlength{\tabcolsep}{" + tabcolsep + "}",
            r"\renewcommand{\arraystretch}{" + arraystretch + "}", ""])

    @staticmethod
    def tabular(spec, headings):
        return (r"\begin{adjustbox}{max width=\linewidth}" + "\n" +
                r"\begin{tabular}{" + spec + "}\n" + r"\toprule" + "\n" +
                Tables.row(headings) + r"\midrule" + "\n")

    @staticmethod
    def endtab():
        return "\n".join([r"\bottomrule", r"\end{tabular}", r"\end{adjustbox}", ""])

    @staticmethod
    def finish(note):
        return "\n".join([r"\par\smallskip", r"\begin{minipage}{0.98\linewidth}\footnotesize",
                          note, r"\end{minipage}", r"\endgroup", r"\end{table}", ""])

    def save(self, filename, content, expected_count):
        assert content.count(r"\begin{table}") == content.count(r"\end{table}") == expected_count
        for env in ["tabular", "adjustbox", "minipage"]:
            assert content.count(r"\begin{" + env + "}") == content.count(r"\end{" + env + "}")
        assert content.count(r"\begingroup") == content.count(r"\endgroup")
        assert not re.search(r"\\begin\{(?:equation|align|gather|multline)\}", content)
        assert not re.search(r"[\u3400-\u9fff]", content)
        path = self.output / filename
        path.write_text(content, encoding="utf-8")
        self.outputs[filename] = dict(tables=expected_count, sha256=digest(path), bytes=path.stat().st_size)

    def oriented_point_auc(self, metric, first, second, T):
        filt = self.auc[(self.auc.metric == metric) & (self.auc["T"] == T)]
        z = filt[(filt["first"] == first) & (filt["second"] == second)].copy()
        reverse = False
        if len(z) == 0:
            z = filt[(filt["first"] == second) & (filt["second"] == first)].copy()
            z["auc_first_smaller"] = 1 - z.auc_first_smaller
            reverse = True
        assert len(z) == 9 and len(z[["p", "q"]].drop_duplicates()) == 9
        assert z.auc_first_smaller.between(0, 1).all()
        rec = dict(metric=metric, first=first, second=second, T=T,
                   minimum=float(z.auc_first_smaller.min()), maximum=float(z.auc_first_smaller.max()),
                   count_ge095=int((z.auc_first_smaller >= .95 - 1e-12).sum()),
                   archived_pair_complemented=reverse,
                   values=z[["p", "q", "auc_first_smaller"]].sort_values(["p", "q"]).to_dict("records"))
        self.auc_records.append(dict(family="finite", **rec))
        return rec

    def cap_auc_row(self, p, q, T, first, second):
        z = self.cap_auc[(self.cap_auc.p == p) & (self.cap_auc.q == q) & (self.cap_auc["T"] == T)]
        selected = z[(z["first"] == first) & (z["second"] == second)]
        reverse = False
        if len(selected) == 0:
            selected = z[(z["first"] == second) & (z["second"] == first)]
            reverse = True
        assert len(selected) == 1, (p, q, T, first, second)
        row = selected.iloc[0].copy()
        if reverse:
            # Reversing X<Y complements each AUC and reverses the paired change.
            for metric in ["E1_center", "E1_cap5_saved", "E1_cap10"]:
                point, low, high = [float(row[metric + s]) for s in ["_auc", "_ci_low", "_ci_high"]]
                row[metric + "_auc"], row[metric + "_ci_low"], row[metric + "_ci_high"] = 1-point, 1-high, 1-low
            for prefix in ["cap10_minus_center", "cap10_minus_cap5"]:
                point, low, high = [float(row[prefix + s]) for s in ["_delta_auc", "_ci_low", "_ci_high"]]
                row[prefix + "_delta_auc"], row[prefix + "_ci_low"], row[prefix + "_ci_high"] = -point, -high, -low
            row["n_first"], row["n_second"] = row.n_second, row.n_first
        row["first"], row["second"] = first, second
        self.auc_records.append(dict(family="cap", p=p, q=q, T=T, first=first, second=second,
                                    archived_pair_complemented=reverse))
        return row

    @staticmethod
    def ci(point, low, high):
        assert all(number(x) for x in [point, low, high])
        assert float(low) <= float(high) + 1e-12
        return r"\shortstack{$" + f"{point:.4f}" + r"$\\{\scriptsize$[" + f"{low:.4f},{high:.4f}" + r"]$}}"

    def finite(self):
        text = "% New lambda-axis-domain searches; generated by build_lambda_manuscript_tables.py.\n"
        text += self.start(
            r"Early finite-sample separation on the coefficient-angle domain $\mathcal D^\lambda_{\delta_0}$ across all nine ordered dimension pairs. Ranges are the minimum and maximum empirical AUC across pairs; counts indicate AUC at least $0.95$.",
            LABEL + "finite:auc", clear=False, arraystretch="1.10")
        text += self.tabular("llcccc", ["Metric", "Fixed direction", r"\multicolumn{2}{c}{$T=100$}", r"\multicolumn{2}{c}{$T=200$}"])
        text += self.row(["", "", "AUC range", r"$\ge0.95$", "AUC range", r"$\ge0.95$"])
        for metric, title in [("rho_min", r"$\rho_{\mathrm{ind},\min}$"),
                              ("E1_at_rho", r"$E_{1,\mathrm{center}}$"),
                              ("E1_cap10", r"$E_{1,10}$"),
                              ("E1_min", r"Independent $(E_1)_{\min}$"),
                              ("J_relative_min", r"$j_{\mathrm{nor}}$"),
                              ("center_to_E1_nearest_deg_090", r"$\Delta_T$")]:
            for first, second in PAIRS:
                vals = [title, "$" + LETTERS[first] + "<" + LETTERS[second] + "$"]
                for T in [100, 200]:
                    z = self.oriented_point_auc(metric, first, second, T)
                    vals += [f'{z["minimum"]:.4f}--{z["maximum"]:.4f}', f'{z["count_ge095"]}/9']
                text += self.row(vals)
            text += r"\addlinespace" + "\n"
        text += self.endtab() + self.finish(
            r"$I$, $B$, and $E$ denote the identifiable interior, identifiable boundary, and non-identifiable designs. For $X<Y$, AUC is the proportion of cross-group score pairs with the first score smaller, with half weight for ties. The same three directions $E<I$, $E<B$, and $B<I$ are used for every metric, sample size, and dimension pair, including weak or reversed orderings. The $0.95$ cutoff is descriptive rather than calibrated. Subscripts $\min$ refer to the lowest feasible numerical candidates found, without a global-optimality certificate. All geometric scores are recomputed on the new $\lambda$-angle domain; $j_{\mathrm{nor}}$ is unchanged because it depends only on the estimated loadings.")
        for reg, title in REGIMES.items():
            text += self.start(r"Early finite-sample levels on $\mathcal D^\lambda_{\delta_0}$: " + title.lower() + r". All nine ordered dimension pairs are retained; entries are mean [sample variance].", LABEL + "finite:" + reg)
            text += self.tabular("crrcccc", [r"$(p,q)$", "$T$", "$n$", r"$\rho_{\mathrm{ind},\min}$", r"Independent $(E_1)_{\min}$", r"$j_{\mathrm{nor}}$", r"\shortstack{$\Delta_T$\\(degrees)}"])
            for p, q in product(DIMS, repeat=2):
                for T in [100, 200]:
                    z = self.pick("key", p, q, reg, T)
                    context = dict(family="finite", p=p, q=q, regime=reg, T=T)
                    vals = [f"$({p},{q})$", T, int(z.n)]
                    for field in ["rho_min", "E1_min", "J_relative_min", "center_to_E1_nearest_deg_090"]:
                        vals.append(self.cell(z, field, context, 3 if "deg" in field else 6))
                    text += self.row(vals)
                text += r"\addlinespace[2pt]" + "\n"
            text += self.endtab() + self.finish(
                r"$\delta_0=\arccos(0.90)\simeq25.84^\circ$ excludes coordinate-axis directions in coefficient space. The loading-matrix cosine is an auxiliary quantity and is not constrained. $\rho_{\mathrm{ind}}=\rho/2$ uses each independent minor once; $\Delta_T$ is the projective $\lambda$ angle from the lowest $\rho_{\mathrm{ind}}$ candidate to the nearest independently found $E_1$ candidate. Values are not squared or multiplied by 100. Variances use the actual count minus one. Only interior $p=q=80,T=100$ has $n=96$; every other cell has $n=100$. QMLE stopping flags do not filter these summaries.")
        self.save("finite_tables.tex", text, 4)

    def count_value(self, k, g, qual, metric):
        n = int(k.n)
        if metric == "available":
            return n
        if metric == "finite_Q":
            count = int(g.pair_joint_max_n_finite)
            assert count == int(g.pair_exists_count)
            assert (count == 0 and pd.isna(k.pair_joint_max_mean)) or abs(k.pair_joint_max_mean-g.pair_joint_max_mean) < 1e-12
            return count
        if metric == "local":
            count = int(qual.local_candidate_exists10_count)
            assert count == int(k.E1_local_within10_found_n_finite)
            return count
        field = {"independent_le5": "center_to_independent_E1_le5_090_count",
                 "qmle_stopping": "qmle_converged_count", "angle_edge": "cap10_angle_edge_count",
                 "axis_edge": "cap10_axis_edge_count", "physical_gt090": "cap10_physical_cosine_gt090_count"}[metric]
        return int(g[field] if field in g else k[field] if field in k else qual[field])

    def complete(self):
        # Two additional physical-cosine rows document that they are not exclusion constraints.
        specs_a = [
            (r"$\rho_{\mathrm{ind},\min}$", "rho_min", 6, False),
            (r"$E_{1,\mathrm{center}}$", "E1_at_rho", 6, False),
            (r"$E_{1,5}$", "E1_cap5_saved", 6, False),
            (r"$E_{1,10}$", "E1_cap10", 6, False),
            (r"$\rho_{\mathrm{ind}}$ at the $10^\circ$ solution", "rho_cap10", 6, False),
            (r"Independent $(E_1)_{\min}$", "E1_min", 6, False),
            (r"$\rho_{\mathrm{ind}}$ at independent $(E_1)_{\min}$", "rho_at_E1", 6, False),
            (r"$E_{1,5}-E_{1,10}$", "E1_gain5_to10", 6, False)]
        specs_b = [
            (r"$j_{\mathrm{orig}}$", "J_min", 6, False),
            (r"$j_{\mathrm{nor}}$", "J_relative_min", 6, False),
            (r"$Q_T(5^\circ)$, finite values", "pair_joint_max", 6, True),
            (r"$\Delta_T$ (degrees)", "center_to_E1_nearest_deg_090", 3, False),
            (r"Center-to-$10^\circ$ solution angle (degrees)", "center_to_cap10_deg", 3, False),
            (r"Local $E_1$ within $10^\circ$, conditional", "E1_local_within10_found", 6, True),
            (r"Physical $c_{\mathrm{mat}}$ at independent $(E_1)_{\min}$", "E1_max_loading_cosine", 6, False),
            (r"Physical $c_{\mathrm{mat}}$ at $E_{1,10}$", "cap10_max_loading_cosine", 6, False),
            (r"$e_Z$", "subspace_error", 6, False),
            (r"$e_\Sigma$", "covariance_relative_error", 6, False)]
        counts = [("Available estimates", "available"), (r"Finite $Q_T$ / available", "finite_Q"),
                  (r"Found local candidate / available", "local"),
                  (r"Independent $\Delta_T\le5^\circ$ / available", "independent_le5"),
                  (r"Strict QMLE stopping / available", "qmle_stopping"),
                  (r"$10^\circ$ angular-cap contact / available", "angle_edge"),
                  (r"$\lambda$-axis exclusion contact / available", "axis_edge"),
                  (r"Auxiliary $c_{\mathrm{mat}}>0.90$ at $E_{1,10}$ / available", "physical_gt090")]
        text = "% Complete new lambda-domain geometry; archived fit-only Jacobians are reused unchanged.\n"
        for p, q, reg in product(DIMS, DIMS, REGIMES):
            text += self.start(r"Complete diagnostics on $\mathcal D^\lambda_{\delta_0}$: " + REGIMES[reg].lower() + f", $p={p}$, $q={q}$. " + r"Panels A--B report mean [sample variance]; Panel C gives the actual denominators and contacts.", LABEL + f"alltab:metrics:p{p}q{q}{reg}", tabcolsep="3.5pt", arraystretch="0.96")
            text += r"\scriptsize" + "\n"
            text += self.tabular("@{}lcccc@{}", ["Metric"] + [f"$T={T}$" for T in TIMES])
            for panel, specs in [("A: polynomial and rank-one residuals", specs_a), ("B: local geometry, correspondence and recovery", specs_b)]:
                text += r"\multicolumn{5}{@{}l}{\textit{Panel " + panel + r"}}\\" + "\n"
                for title, metric, decimals, conditional in specs:
                    vals = []
                    for T in TIMES:
                        z = self.pick("key", p, q, reg, T)
                        vals.append(self.cell(z, metric, dict(family="complete", p=p, q=q, regime=reg, T=T), decimals, conditional))
                    text += self.row([title] + vals)
                text += r"\addlinespace[2pt]" + "\n"
            text += r"\multicolumn{5}{@{}l}{\textit{Panel C: denominators and constraint contacts}}\\" + "\n"
            for title, metric in counts:
                vals = []
                for T in TIMES:
                    k, g, qual = [self.pick(name, p, q, reg, T) for name in ["key", "geom", "quality"]]
                    n = int(k.n)
                    assert n == int(g.n_seeds) == int(qual.n_seeds)
                    count = self.count_value(k, g, qual, metric)
                    assert 0 <= count <= n
                    vals.append(str(count) if metric == "available" else f"{count}/{n}")
                    self.count_records.append(dict(p=p, q=q, regime=reg, T=T, metric=metric, count=count, available=n))
                text += self.row([title] + vals)
            text += self.endtab() + self.finish(
                r"All geometric searches use $\max_j|\lambda_j|/\|\lambda\|\le0.90$ and projective coefficient angles. Physical $c_{\mathrm{mat}}$ values are auxiliary diagnostics and may exceed $0.90$. Subscripts $\min$ denote the lowest feasible candidates found; numerical stopping does not certify mathematical local or global optimality. $\rho_{\mathrm{ind}}=\rho/2$; $J$ uses the $\psi$-column convention. $Q_T$ and conditional local $E_1$ use only finite/found cases; unmatched $Q_T=+\infty$ is retained. A dash means no found value; NA means variance unavailable at fewer than two values. Variances use actual counts minus one. No imputation or strict-stopping filter is applied. QMLE estimates and budgets are unchanged; only geometric searches are new.")
        # These spectra are entirely fit-dependent and therefore independent of either off-axis domain.
        for p, q, reg in product(DIMS, DIMS, REGIMES):
            chunks = [(1, 3), (4, 6), (7, 8)] if reg == "interior" else [(1, 5)]
            nsv = 7 if reg == "interior" else 4
            for first, last in chunks:
                text += self.start(r"Complete $J_\ell$ spectra: " + REGIMES[reg].lower() + f", $p={p}$, $q={q}$, " + rf"$\ell={first},\ldots,{last}$. " + r"Entries are mean [sample variance]; singular values are in ascending order.", LABEL + f"alltab:spectra:p{p}q{q}{reg}:{first}-{last}", tabcolsep="2.2pt")
                text += r"\begin{adjustbox}{max width=\linewidth}" + "\n" + r"\begin{tabular}{@{}rr*{8}{c}@{}}" + "\n" + r"\toprule" + "\n"
                text += self.row([r"$\ell$", "$k$"] + [r"\multicolumn{2}{c}{$T=" + str(T) + "$}" for T in TIMES])
                text += r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}\cmidrule(lr){7-8}\cmidrule(lr){9-10}" + "\n"
                text += self.row(["", ""] + [r"\texttt{orig}", r"\texttt{nor}"] * 4) + r"\midrule" + "\n"
                for ell, k in product(range(first, last+1), range(1, nsv+1)):
                    if k == 1 and ell > first:
                        text += r"\addlinespace[2pt]" + "\n"
                    vals = [ell, k]
                    for T in TIMES:
                        df = self.J[(self.J.p == p) & (self.J.q == q) & (self.J.regime == reg) & (self.J["T"] == T) & (self.J.ell == ell) & (self.J.svd_index == k)]
                        assert len(df) == 1
                        z = df.iloc[0]
                        expected_n = 96 if (p, q, reg, T) == (80, 80, "interior", 100) else 100
                        assert int(z.n) == int(z.orig_n_finite) == int(z.nor_n_finite) == expected_n
                        for form in ["orig", "nor"]:
                            vals.append(self.cell(z, form, dict(family="spectra", p=p, q=q, regime=reg, T=T, ell=ell, svd_index=k)))
                    text += self.row(vals)
                text += self.endtab() + self.finish(
                    r"\texttt{orig} gives $s_{\ell,k}$; \texttt{nor} gives $s_{\ell,k}/s_{\ell,d-1}$, normalized within each seed before summarization. Values are unsquared and are not multiplied by 100. The $\psi$-column convention is used; the literal ordered-minor derivative doubles all \texttt{orig} values. These fit-only spectra are reused unchanged from the saved QMLE estimates. Each cell has 100 seeds except interior $p=q=80,T=100$, which has 96.")
        self.save("complete_tables.tex", text, 72)

    def auc_tables(self):
        text = "% Fixed-direction AUCs from new lambda-domain residuals; paired seed-block bootstrap intervals.\n"
        for p, q in product(DIMS, repeat=2):
            text += self.start(f"Complete separation on the coefficient-angle domain, $(p,q)=({p},{q})$. " + r"Score columns give AUC [95\% paired bootstrap interval]; the final column gives the paired AUC difference and its interval.", LABEL + f"allauc:p{p}q{q}", tabcolsep="3.5pt", arraystretch="1.10")
            text += self.tabular("lrc" + "c"*4, ["Comparison", "$T$", r"$n_1/n_2$", "Center", r"$5^\circ$", r"$10^\circ$", r"$\mathrm{AUC}_{10}-\mathrm{AUC}_{5}$"])
            for first, second in PAIRS:
                for T in TIMES:
                    z = self.cap_auc_row(p, q, T, first, second)
                    vals = ["$" + LETTERS[first] + "<" + LETTERS[second] + "$", T, f"{int(z.n_first)}/{int(z.n_second)}"]
                    for prefix in ["E1_center", "E1_cap5_saved", "E1_cap10"]:
                        vals.append(self.ci(z[prefix+"_auc"], z[prefix+"_ci_low"], z[prefix+"_ci_high"]))
                    vals.append(self.ci(z.cap10_minus_cap5_delta_auc, z.cap10_minus_cap5_ci_low, z.cap10_minus_cap5_ci_high))
                    text += self.row(vals)
                text += r"\addlinespace" + "\n"
            text += self.endtab() + self.finish(
                r"Directions are fixed as $E<I$, $E<B$, and $B<I$ throughout, with half weight for ties. Center, $5^\circ$, and $10^\circ$ denote $E_1$ at the lowest $\rho_{\mathrm{ind}}$ candidate and in its angular neighborhoods on the common $\mathcal D^\lambda_{\delta_0}$ domain. Reversing an archived comparison complements its AUC and interval and negates its paired change with reversed interval endpoints. Intervals use 1,000 shared-seed block bootstrap draws, without multiplicity adjustment. Empirical AUC equal to one does not establish zero population error. Actual group sizes are shown.")
        self.save("auc_tables.tex", text, 9)

    def selected(self):
        text = "% Selected new lambda-domain comparisons; no unrecomputed continuation-geometry table.\n"
        text += self.start(r"Rank-one errors at the lowest $\rho_{\mathrm{ind}}$ candidate and in its angular neighborhoods on $\mathcal D^\lambda_{\delta_0}$, $p=q=40$. Entries are mean [sample variance] over 100 seeds.", LABEL + "selectedtab:cap_comparison")
        text += self.tabular("lrccc", ["Design", "$T$", r"$E_{1,\mathrm{center}}$", "$E_{1,5}$", "$E_{1,10}$"])
        for reg in REGIMES:
            for T in TIMES:
                z = self.pick("cap", 40, 40, reg, T)
                vals = [LETTERS[reg], T] + [self.cell(z, field, dict(family="selected_cap", p=40, q=40, regime=reg, T=T)) for field in ["E1_center", "E1_cap5_saved", "E1_cap10"]]
                text += self.row(vals)
            text += r"\addlinespace" + "\n"
        text += self.endtab() + self.finish(
            r"The axis-exclusion angle is $\delta_0\simeq25.84^\circ$, distinct from the $5^\circ$ and $10^\circ$ neighborhood radii. The center and feasible $5^\circ$ witnesses are retained in the $10^\circ$ search, so $E_{1,10}\le E_{1,5}\le E_{1,\mathrm{center}}$ by construction. Each physical loading-matrix cosine is auxiliary and is not restricted by this domain.")

        text += self.start(r"Complementary local and off-axis geometry on $\mathcal D^\lambda_{\delta_0}$ at $p=q=40,T=800$. Entries are mean [sample variance]; counts refer to 100 seeds.", LABEL + "selectedtab:cap_correspondence", tabcolsep="3pt")
        text += r"\textit{Panel A: Jacobians and independent-pair correspondence.}\par\smallskip" + "\n"
        text += self.tabular("lccccc", ["Design", r"$j_{\mathrm{orig}}$", r"$j_{\mathrm{nor}}$", r"$Q_T(5^\circ)$", r"Finite $Q_T$", r"$\Delta_T$ (deg.)"])
        for reg in REGIMES:
            z = self.pick("key", 40, 40, reg, 800)
            context = dict(family="selected_correspondence", p=40, q=40, regime=reg, T=800)
            text += self.row([LETTERS[reg]] + [self.cell(z, field, context, conditional=(field == "pair_joint_max")) for field in ["J_min", "J_relative_min", "pair_joint_max"]] + [f"{int(z.pair_joint_max_n_finite)}/100", self.cell(z, "center_to_E1_nearest_deg_090", context, 3)])
        text += self.endtab() + r"\par\smallskip" + "\n"
        text += r"\textit{Panel B: cap contacts and uncapped local candidates.}\par\smallskip" + "\n"
        text += self.tabular("lccccc", ["Design", r"Cap angle (deg.)", r"Angle edge", r"$\lambda$-axis edge", "Local candidate", r"Conditional local $E_1$"])
        for reg in REGIMES:
            z = self.pick("key", 40, 40, reg, 800)
            context = dict(family="selected_correspondence", p=40, q=40, regime=reg, T=800)
            text += self.row([LETTERS[reg], self.cell(z, "center_to_cap10_deg", context, 3), f"{int(z.cap10_angle_edge_count)}/100", f"{int(z.cap10_axis_edge_count)}/100", f"{int(z.local_candidate_exists10_count)}/100", self.cell(z, "E1_local_within10_found", context, conditional=True)])
        text += self.endtab() + r"\par\smallskip" + "\n"
        text += r"\textit{Panel C: physical cosines, reported without exclusion.}\par\smallskip" + "\n"
        text += self.tabular("lccc", ["Design", r"$c_{\mathrm{mat}}$ at independent $(E_1)_{\min}$", r"$c_{\mathrm{mat}}$ at $E_{1,10}$", r"$c_{\mathrm{mat}}>0.90$ at $E_{1,10}$"])
        for reg in REGIMES:
            z = self.pick("key", 40, 40, reg, 800)
            context = dict(family="selected_physical", p=40, q=40, regime=reg, T=800)
            text += self.row([LETTERS[reg], self.cell(z, "E1_max_loading_cosine", context), self.cell(z, "cap10_max_loading_cosine", context), f"{int(z.cap10_physical_cosine_gt090_count)}/100"])
        text += self.endtab() + r"\par\smallskip" + "\n"
        text += r"\textit{Panel D: all four original singular values of boundary $J_1$.}\par\smallskip" + "\n"
        text += self.tabular("lcccc", ["$T$", "$s_{1,1}$", "$s_{1,2}$", "$s_{1,3}$", "$s_{1,4}$"])
        for T in [100, 800]:
            df = self.J[(self.J.p == 40) & (self.J.q == 40) & (self.J.regime == "boundary") & (self.J["T"] == T) & (self.J.ell == 1)].sort_values("svd_index")
            assert len(df) == 4
            text += self.row([T] + [self.cell(z, "orig", dict(family="selected_J", p=40, q=40, regime="boundary", T=T, ell=1, svd_index=int(z.svd_index))) for _, z in df.iterrows()])
        text += self.row(["Population"] + ["$0$" if int(z.svd_index) == 1 else f"${z.reference_orig_mean:.6f}$" for _, z in df.iterrows()])
        text += self.endtab() + self.finish(
            r"$Q_T$ uses independently searched $\rho_{\mathrm{ind}}$/$E_1$ pairs within $5^\circ$, with unmatched cases recorded as $+\infty$ and excluded only from its finite mean and variance. The cap solution has an imposed angular proximity; the local-candidate count additionally requires an uncapped feasible candidate within $10^\circ$. Numerical stopping is not a proof of mathematical local optimality. Physical cosine exceeding $0.90$ is allowed on the new domain. The boundary Jacobian spectra are unchanged fit-only quantities; all values shown are unsquared.")

        text += self.start(r"Interior--exterior separation of the new coefficient-angle-domain errors. AUC is fixed as the probability that an exterior score is smaller than an interior score, with half weight for ties.", LABEL + "selectedtab:cap_auc")
        text += r"\textit{Panel A: all nine dimension pairs at $T=800$.}\par\smallskip" + "\n"
        text += self.tabular("rrccc", ["$p$", "$q$", r"$E_{1,\mathrm{center}}$", "$E_{1,5}$", "$E_{1,10}$"])
        for p, q in product(DIMS, repeat=2):
            z = self.cap_auc_row(p, q, 800, "exterior", "interior")
            text += self.row([p, q] + [f"{z[field + '_auc']:.4f}" for field in ["E1_center", "E1_cap5_saved", "E1_cap10"]])
        text += self.endtab() + r"\par\smallskip" + "\n"
        text += r"\textit{Panel B: selected paired changes from $5^\circ$ to $10^\circ$.}\par\smallskip" + "\n"
        text += self.tabular("rrrcc", ["$p$", "$q$", "$T$", r"$\Delta\mathrm{AUC}_{10-5}$", r"95\% bootstrap interval"])
        for p, q, T in [(40, 40, 400), (20, 20, 800), (40, 40, 800)]:
            z = self.cap_auc_row(p, q, T, "exterior", "interior")
            text += self.row([p, q, T, f"{z.cap10_minus_cap5_delta_auc:.4f}", f"[{z.cap10_minus_cap5_ci_low:.4f},{z.cap10_minus_cap5_ci_high:.4f}]"])
        text += self.endtab() + self.finish(
            r"Every $T=800$ group has 100 seeds. Intervals use 1,000 shared-seed block bootstrap resamples and are exploratory, without multiplicity adjustment. The nested feasible sets guarantee nonincreasing pointwise errors, not an increase in AUC. AUC equal to one describes empirical ordering in these fixed designs. These AUC summaries use the original saved fits. The separately recomputed continued-fit sensitivity is reported in its own audit.")
        self.save("selected_tables.tex", text, 3)

    def index(self):
        text = "\n".join([r"\par\medskip\noindent\textit{Index to the complete coefficient-angle-domain grid.}", r"\begingroup", r"\spacingset{1.0}", r"\small\begin{center}", r"\begin{tabular}{crrrrrrr}", r"\toprule", r"& \multicolumn{3}{c}{Metric tables} & \multicolumn{3}{c}{Spectrum tables} & AUC\\", r"$(p,q)$ & $I$ & $B$ & $E$ & $I$ (first part) & $B$ & $E$ & All $T$\\", r"\midrule", ""])
        for p, q in product(DIMS, repeat=2):
            vals = [f"$({p},{q})$"]
            vals += [r"\ref{" + LABEL + f"alltab:metrics:p{p}q{q}{reg}" + "}" for reg in REGIMES]
            vals += [r"\ref{" + LABEL + f"alltab:spectra:p{p}q{q}{reg}:" + ("1-3" if reg == "interior" else "1-5") + "}" for reg in REGIMES]
            vals += [r"\ref{" + LABEL + f"allauc:p{p}q{q}" + "}"]
            text += self.row(vals)
        text += "\n".join([r"\bottomrule", r"\end{tabular}\end{center}\endgroup", r"\noindent All metric and spectrum groups report $T=100,200,400,800$. Interior spectra continue for $\ell=4$--$6$ and $7$--$8$ in the next two tables. $I$, $B$, and $E$ denote the three designs. The Jacobian spectra depend only on the unchanged fitted loadings; all off-axis geometric scores use $\mathcal D^\lambda_{\delta_0}$.", ""])
        self.save("index.tex", text, 0)

    def run(self):
        self.output.mkdir(parents=True, exist_ok=True)
        self.finite()
        self.complete()
        self.auc_tables()
        self.selected()
        self.index()
        labels = [x["label"] for x in self.tables]
        assert len(labels) == len(set(labels)) == 88
        alltext = "\n".join((self.output / name).read_text(encoding="utf-8") for name in self.outputs)
        assert alltext.count(r"\begin{table}") == 88
        parsed = re.findall(r"\\shortstack\{([-+eE.0-9]+)\\\\\[-1pt\]\{\\scriptsize\[([-+eE.0-9]+)\]\}\}", alltext)
        numeric = [x for x in self.mv_records if "display_mean" in x]
        assert parsed == [(x["display_mean"], x["display_variance"]) for x in numeric]
        audit = dict(status="passed", domain="lambda-projective-axis-angle", cosine_cutoff=.90,
                     axis_angle_degrees=math.degrees(math.acos(.90)), new_geometry=True,
                     original_qmle_fits_unchanged=True, J_spectra_reused_unchanged=True,
                     continued_fit_audit_included_in_this_table_set=False,
                     continued_fit_audit_file_exists=(self.output/"iteration_sensitivity.tex").is_file(),
                     available_fits=10796, cell_count=108,
                     table_count=len(labels), output_files=self.outputs,
                     source_hashes={str(p): digest(p) for p in self.paths.values()},
                     labels=labels, mean_variance_cells=len(self.mv_records),
                     mean_variance=self.mv_records, counts=self.count_records, auc=self.auc_records)
        (self.output / "tables_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
        print(json.dumps({k: audit[k] for k in ["status", "table_count", "cell_count", "mean_variance_cells", "domain"]}))
        print(str(self.output))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=ROOT / "outputs/lambda_axis_geometry_20261030")
    parser.add_argument("--output", type=Path, default=PROJECT / "additions/lambda")
    parser.add_argument("--spectra", type=Path, default=SPECTRA_DEFAULT)
    parser.add_argument("--wait-seconds", type=float, default=0,
                        help="Wait at most this long for aggregate CSV files, checking every 5 seconds.")
    args = parser.parse_args()
    deadline = time.monotonic() + args.wait_seconds
    while True:
        missing = [args.data_root / "tables" / name for name in REQUIRED if not (args.data_root / "tables" / name).is_file()]
        if not missing:
            break
        if time.monotonic() >= deadline:
            raise FileNotFoundError("New-domain aggregation outputs are not ready: " + ", ".join(str(p) for p in missing))
        time.sleep(min(5, max(0, deadline-time.monotonic())))
    Tables(args.data_root, args.output, args.spectra).run()


if __name__ == "__main__":
    main()
