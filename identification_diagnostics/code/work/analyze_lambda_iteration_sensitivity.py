"""Pair saved continued-fit lambda geometry with the full original lambda grid."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from build_lambda_manuscript_tables import Tables, LABEL, PROJECT

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "outputs/lambda_axis_geometry_20261030"
OUT = BASE / "iteration_sensitivity"
KEY = ["p", "q", "regime", "seed", "T"]
METRICS = ["rho_min", "E1_at_rho", "E1_min", "rho_at_E1", "E1_cap_090", "E1_cap10",
           "rho_cap10", "center_to_E1_nearest_deg_090", "pair_joint_max", "center_to_cap10_deg",
           "cap10_max_loading_cosine", "cap10_max_lambda_cosine", "E1_local_within10_found",
           "local_best_angle10", "J_min", "J_relative_min"]
PAIRS = [("exterior", "interior"), ("exterior", "boundary"), ("boundary", "interior")]


def build_iteration_table(summary, findings, direct):
    writer=Tables.__new__(Tables)
    writer.output=PROJECT/"additions/lambda"
    writer.mv_records=[];writer.tables=[];writer.outputs={}
    text=(r"\par\medskip\noindent\textit{Additional-iteration sensitivity on the new domain.} "
          r"We recomputed the $\mathcal D^\lambda_{\delta_0}$ searches on 28 previously saved continued QMLE fits, "
          r"after a common additional budget of 1,300 iterations, and on three separately saved adaptive stress fits. "
          r"No likelihood optimization was performed in this reassessment. The original and continued fits use the same new-domain search budgets and case-specific random seeds. "
          r"These configurations were selected purposively; the 15 $T=800$ cases comprise six boundary, four exterior, and five interior cases, pooled across their audited dimensions. They are not a fresh Monte Carlo sample or a higher-budget rerun of the full grid. "
          r"Table~\ref{"+LABEL+r"selectedtab:cap_iteration} reports their paired residual changes."+"\n\n")
    text+=(r"Replacing only these 15 records in the full $T=800$ grid leaves every center/$5^\circ$/$10^\circ$ AUC and its cross-regime ordering matrix unchanged. "
           r"The $j_{\mathrm{orig}}$, $j_{\mathrm{nor}}$, and $\Delta_T$ AUCs are also unchanged. "
           r"The largest AUC changes among $\rho_{\mathrm{ind},\min}$, independent $(E_1)_{\min}$, and $Q_T$ are $0.0006$, $0.0001$, and $0.0005$, respectively. "
           r"Thus the residual-based $T=800$ separation of these fixed designs survives the audited continuation, even though some boundary residuals move appreciably. "
           r"The smaller-sample stress results are less stable: for the 13 common $T=100$ cases the largest individual $E_{1,10}$ change is $0.01472$, "
           r"and the corresponding largest full-grid AUC change is $0.0029$ ($0.0099$ for some other scores). "
           r"The three adaptive exterior stress cases remain a separate comparison. No continued fits were available at $T=200$ or $400$, "
           r"so this audit supports the stated large-$T$ empirical separation without establishing optimization adequacy for every configuration."+"\n\n")
    text+=writer.start(r"Paired original and continued-fit diagnostics on $\mathcal D^\lambda_{\delta_0}$. Entries are mean [sample variance] over the same audited fits; the final column is the largest individual absolute change.",LABEL+"selectedtab:cap_iteration",tabcolsep="4pt")
    text+=r"\textit{Panel A: 15 common-budget $T=800$ cases.}\par\smallskip"+"\n"
    text+=writer.tabular("llrccc",["Design","Score","$n$","Original","Continued",r"$\max|\mathrm{change}|$"])
    labels={"interior":"Interior","boundary":"Boundary","exterior":"Exterior"}
    for reg in labels:
        for metric,title in [("E1_at_rho",r"$E_{1,\mathrm{center}}$"),("E1_cap_090",r"$E_{1,5}$"),("E1_cap10",r"$E_{1,10}$")]:
            z=summary[(summary.subset=="T800_common")&(summary.regime==reg)&(summary.metric==metric)]
            assert len(z)==1
            z=z.iloc[0]
            context=dict(subset="T800_common",regime=reg,metric=metric,n=int(z.n))
            text+=writer.row([labels[reg],title,int(z.n),writer.mv(z.old_mean,z.old_variance,dict(context,version="original")),writer.mv(z.new_mean,z.new_variance,dict(context,version="continued")),f"{z.max_absolute_difference:.6f}"])
        text+=r"\addlinespace"+"\n"
    text+=writer.endtab()+r"\par\smallskip"+"\n"
    text+=r"\textit{Panel B: $T=100$ common-budget and separate adaptive stress cases.}\par\smallskip"+"\n"
    text+=writer.tabular("llrccc",["Comparison","Design","$n$",r"Original $E_{1,10}$",r"Continued $E_{1,10}$",r"$\max|\mathrm{change}|$"])
    for subset,reg in [("T100_common","interior"),("T100_common","boundary"),("T100_common","exterior"),("adaptive3","exterior")]:
        z=summary[(summary.subset==subset)&(summary.regime==reg)&(summary.metric=="E1_cap10")]
        assert len(z)==1
        z=z.iloc[0]
        context=dict(subset=subset,regime=reg,metric="E1_cap10",n=int(z.n))
        text+=writer.row(["Common" if subset=="T100_common" else "Adaptive",labels[reg],int(z.n),writer.mv(z.old_mean,z.old_variance,dict(context,version="original")),writer.mv(z.new_mean,z.new_variance,dict(context,version="continued")),f"{z.max_absolute_difference:.6f}"])
    text+=writer.endtab()+writer.finish(
        r"Panel A uses the 15 common-budget configurations at $T=800$; Panel B uses 13 different common-budget configurations at $T=100$ and three additional adaptive continuation records for selected exterior stress cases. Adaptive and common rows are not pooled. The original full-grid means are retained in all main tables. Each diagnostic was newly optimized on the coefficient-angle domain; no historical matrix-domain geometric score is substituted. The same $E_1$ and singular-value scales are used throughout, without squaring or multiplication by 100. The common-budget continuation is additional to the original regime-specific QMLE caps. No continued-fit audit exists at $T=200,400$.")
    writer.save("iteration_sensitivity.tex",text,1)
    audit=dict(status="passed",tables=1,label=LABEL+"selectedtab:cap_iteration",source=str(OUT/"paired_lambda_iteration_mean_variance.csv"),
               mean_variance=writer.mv_records,output_files=writer.outputs,direct_validation=direct)
    (writer.output/"iteration_sensitivity_audit.json").write_text(json.dumps(clean(audit),indent=2),encoding="utf-8")


def finite(x):
    a = np.asarray(x, dtype=float)
    return a[np.isfinite(a)]


def stats(x):
    a = finite(x)
    return dict(n=len(a), mean=float(a.mean()) if len(a) else None,
                variance=float(a.var(ddof=1)) if len(a)>1 else None)


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        return float(x) if np.isfinite(x) else None
    return x


def main():
    original = pd.read_csv(BASE / "primary_by_seed_enriched.csv")
    paired = pd.read_csv(OUT / "paired_lambda_iteration_by_seed.csv")
    continued = pd.read_csv(OUT / "continued_lambda_primary.csv")
    assert len(original) == 10796 and len(paired) == len(continued) == 31
    common = paired[paired.comparison_group.eq("common_stage2")]
    adaptive = paired[paired.comparison_group.eq("adaptive_final")]
    assert len(common) == 28 and len(adaptive) == 3
    assert len(common[common["T"].eq(800)]) == 15
    assert common.groupby("regime").size().to_dict() == {"boundary":9,"exterior":10,"interior":9}
    assert paired.core_fit_gram_max_error.max() < 1e-10
    assert not paired.duplicated(KEY + ["comparison_group"]).any()
    # Baseline columns must agree exactly with the newly recomputed full grid.
    join = paired.merge(original[KEY+METRICS], on=KEY, validate="many_to_one")
    baseline_diffs = {}
    for metric in METRICS:
        a,b = join["old_"+metric].to_numpy(),join[metric].to_numpy()
        assert np.array_equal(np.isfinite(a),np.isfinite(b)), metric
        d = np.abs(a[np.isfinite(a)]-b[np.isfinite(b)])
        assert not len(d) or d.max() < 1e-11, (metric,d.max())
        baseline_diffs[metric] = float(d.max()) if len(d) else 0
    rows=[]
    for subset,frame in [("all28_common", common), ("T800_common",common[common["T"].eq(800)]),
                         ("T100_common",common[common["T"].eq(100)]), ("adaptive3",adaptive)]:
        for reg,g in frame.groupby("regime"):
            for metric in METRICS:
                a,b = g["old_"+metric].to_numpy(),g["new_"+metric].to_numpy()
                both = np.isfinite(a)&np.isfinite(b)
                d = b[both]-a[both]
                sa,sb = stats(a),stats(b)
                rows.append(dict(subset=subset,regime=reg,metric=metric,n=len(g),
                                 old_n=sa["n"],new_n=sb["n"],paired_n=int(both.sum()),
                                 old_mean=sa["mean"],old_variance=sa["variance"],
                                 new_mean=sb["mean"],new_variance=sb["variance"],
                                 difference_mean=float(d.mean()) if len(d) else np.nan,
                                 max_absolute_difference=float(np.abs(d).max()) if len(d) else np.nan,
                                 median_absolute_difference=float(np.median(np.abs(d))) if len(d) else np.nan,
                                 finite_indicator_changed=int((np.isfinite(a)!=np.isfinite(b)).sum())))
    summary=pd.DataFrame(rows)
    summary.to_csv(OUT/"paired_lambda_iteration_mean_variance.csv",index=False)
    orderrows=[]
    # Replace only audited fits; all other Monte Carlo estimates stay fixed.
    metrics=["rho_min","E1_at_rho","E1_cap_090","E1_cap10","E1_min",
             "J_min","J_relative_min","center_to_E1_nearest_deg_090","pair_joint_max"]
    for scenario in ["common_stage2","common_then_adaptive3"]:
        replacements=continued[continued.comparison_group.eq("common_stage2")].copy()
        if scenario=="common_then_adaptive3":
            replacements=pd.concat([replacements,continued[continued.comparison_group.eq("adaptive_final")]],ignore_index=True).drop_duplicates(KEY,keep="last")
        repl=replacements.set_index(KEY)
        updated=original.set_index(KEY).copy()
        updated.loc[repl.index,metrics]=repl[metrics]
        updated=updated.reset_index()
        for p,q,T in original[["p","q","T"]].drop_duplicates().itertuples(index=False,name=None):
            oa=original[(original.p==p)&(original.q==q)&(original["T"]==T)]
            na=updated[(updated.p==p)&(updated.q==q)&(updated["T"]==T)]
            for first,second in PAIRS:
                x=oa[oa.regime.eq(first)].sort_values("seed")
                y=oa[oa.regime.eq(second)].sort_values("seed")
                nx=na[na.regime.eq(first)].sort_values("seed")
                ny=na[na.regime.eq(second)].sort_values("seed")
                for metric in metrics:
                    a,b=x[metric].to_numpy(),y[metric].to_numpy()
                    c,d=nx[metric].to_numpy(),ny[metric].to_numpy()
                    # Q legitimately retains +Inf; conditional NaNs are not among these scores.
                    assert not any(np.isnan(z).any() for z in [a,b,c,d]), metric
                    old=(a[:,None]<b[None,:]).astype(float)+.5*(a[:,None]==b[None,:])
                    new=(c[:,None]<d[None,:]).astype(float)+.5*(c[:,None]==d[None,:])
                    orderrows.append(dict(scenario=scenario,p=p,q=q,T=T,first=first,second=second,
                                          metric=metric,n_first=len(a),n_second=len(b),
                                          old_auc=float(old.mean()),new_auc=float(new.mean()),
                                          change_auc=float(new.mean()-old.mean()),
                                          ordering_matrix_entries_changed=int((old!=new).sum()),
                                          ordering_matrix_size=old.size))
    auc=pd.DataFrame(orderrows)
    auc.to_csv(OUT/"lambda_auc_replacement_sensitivity.csv",index=False)
    flags=[]
    for subset,frame in [("all28_common",common),("T800_common",common[common["T"].eq(800)]),("adaptive3",adaptive)]:
        for reg,g in frame.groupby("regime"):
            for field in ["qmle_converged","pair_exists","cap10_angle_edge","cap10_axis_edge","local_candidate_exists10"]:
                flags.append(dict(subset=subset,regime=reg,flag=field,n=len(g),
                                  old_count=int(g["old_"+field].sum()),new_count=int(g["new_"+field].sum()),
                                  indicators_changed=int((g["old_"+field]!=g["new_"+field]).sum())))
    pd.DataFrame(flags).to_csv(OUT/"lambda_iteration_flag_changes.csv",index=False)
    findings=dict(status="passed",domain="lambda_axis_cosine_le_0.90",no_qmle_performed=True,
                  common_cases=28,common_T800_cases=15,adaptive_stress_cases=3,
                  original_lambda_baseline_max_errors=baseline_diffs,
                  core_fit_gram_max_error=float(paired.core_fit_gram_max_error.max()),
                  new_geometry_search_seconds=float(paired.seconds.sum()),
                  T800_summary=summary[summary.subset.eq("T800_common") & summary.metric.isin(["rho_min","E1_at_rho","E1_cap_090","E1_cap10","E1_min","J_min","J_relative_min","center_to_E1_nearest_deg_090"])].to_dict("records"),
                  T800_auc_changes=auc[(auc["T"]==800)&auc.scenario.eq("common_stage2")].groupby("metric").agg(max_absolute_auc_change=("change_auc",lambda x:abs(x).max()),changed_ordering_entries=("ordering_matrix_entries_changed","sum")).reset_index().to_dict("records"),
                  allT_auc_changes=auc[auc.scenario.eq("common_stage2")].groupby(["T","metric"]).agg(max_absolute_auc_change=("change_auc",lambda x:abs(x).max()),changed_ordering_entries=("ordering_matrix_entries_changed","sum")).reset_index().to_dict("records"),
                  selection_warning="Purposively selected audit configurations, not a fresh representative Monte Carlo sample; no continuation records at T=200 or T=400.")
    direct_points=pd.read_csv(OUT/"lambda_iteration_direct_point_validation.csv")
    direct_cases=pd.read_csv(OUT/"lambda_iteration_direct_case_validation.csv")
    assert len(direct_cases)==31 and direct_cases.source_unchanged.all()
    direct=dict(status="passed",continued_fit_records=len(direct_cases),
                unique_configuration_count=len(paired[KEY].drop_duplicates()),
                common_records=28,separate_adaptive_records=3,points=len(direct_points),
                maximum_Frobenius_error=float(direct_points.unit_Frobenius_error.max()),
                maximum_E1_direct_error=float(direct_points.E1_direct_error.max()),
                maximum_rho_direct_error=float(direct_points.rho_direct_error.max()),
                maximum_lambda_cosine=float(direct_points.max_lambda_cosine.max()),
                physical_cosine_gt090_points=int((direct_points.max_physical_cosine>.90).sum()),
                maximum_Q_direct_error=float(direct_cases.Q_direct_error.max()),
                maximum_Delta_direct_error=float(direct_cases.Delta_direct_error.max()),
                source_checksums_unchanged=True)
    assert direct["maximum_Frobenius_error"]<1e-5 and direct["maximum_E1_direct_error"]<2e-6 and direct["maximum_rho_direct_error"]<2e-6
    assert direct["maximum_lambda_cosine"]<=.90+2e-6
    (OUT/"lambda_iteration_direct_validation.json").write_text(json.dumps(clean(direct),indent=2),encoding="utf-8")
    findings["direct_validation"]=direct
    (OUT/"lambda_iteration_findings.json").write_text(json.dumps(clean(findings),indent=2),encoding="utf-8")
    build_iteration_table(summary,findings,direct)
    print(json.dumps({k:findings[k] for k in ["status","common_cases","common_T800_cases","adaptive_stress_cases","core_fit_gram_max_error","new_geometry_search_seconds"]}))
    print(summary[(summary.subset=="T800_common")&summary.metric.isin(["rho_min","E1_at_rho","E1_cap_090","E1_cap10","E1_min","J_relative_min"])][["regime","metric","n","old_mean","new_mean","max_absolute_difference"]].to_string(index=False))
    print(pd.DataFrame(findings["T800_auc_changes"]).to_string(index=False))


if __name__=="__main__":
    main()
