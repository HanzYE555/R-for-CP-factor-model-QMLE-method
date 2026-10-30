"""Independent validation of every new lambda-axis geometry export.

Reads source checkpoints and geometry exports only. Writes audit files under
the new experiment's validation folder; never changes fits or reported data.
Direct matrix-SVD/minor validation is separately performed by the R validator.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KEY = ["p", "q", "regime", "seed", "T"]
CELL = ["p", "q", "regime", "T"]
OUT = ROOT / "outputs/lambda_axis_geometry_20261030"


def flag(v):
    return v.astype(str).str.upper().isin(["TRUE", "1"])


def require(condition, message):
    if not bool(np.all(condition)):
        raise AssertionError(message)


def keyset(frame):
    return set(frame[KEY].itertuples(index=False, name=None))


def projective_angle(a, b):
    aa = np.abs((a*b).sum(axis=1))
    den = np.linalg.norm(a, axis=1)*np.linalg.norm(b, axis=1)
    require(den > 0, "Zero lambda vector in angle comparison")
    return np.degrees(np.arccos(np.clip(aa/den, 0, 1)))


def lambdas(frame, prefix):
    cols = [f"{prefix}_{i:02d}" for i in range(1, 9)]
    require(all(c in frame for c in cols), f"Missing lambda columns: {prefix}")
    raw = frame[cols].to_numpy(dtype=float)
    # Only d<8 padding may be absent; missing coordinates within the actual
    # model are never silently replaced with zero.
    width = np.where(frame.regime.eq("interior"), 8, 5)
    used = np.arange(8)[None, :] < width[:, None]
    require(np.isfinite(raw[used]), f"Missing/non-finite active lambda coordinate: {prefix}")
    require(np.isnan(raw[~used]) | (raw[~used] == 0), f"Nonzero lambda padding: {prefix}")
    ans = frame[cols].fillna(0).to_numpy(dtype=float)
    require(np.isfinite(ans), f"Non-finite lambda vector: {prefix}")
    return ans


def audit_lambda(frame, prefix, reported=None):
    a = lambdas(frame, prefix)
    norms = np.linalg.norm(a, axis=1)
    require(norms > 0, f"Zero lambda vector: {prefix}")
    cos = np.max(np.abs(a), axis=1)/norms
    require(cos <= .900002, f"Lambda-axis exclusion violated: {prefix}")
    error = 0.
    if reported:
        error = float(np.max(np.abs(cos-frame[reported].to_numpy(float))))
        require(error < 2e-7, f"Reported axis cosine mismatch: {reported}")
    return a, cos, error


def inputs(out, pilot, aggregate):
    result = {}
    for kind, aggregate_name in [("primary", "primary_by_seed.csv"),
                                 ("caps", "all_rho_by_seed.csv"),
                                 ("pairs", "minima_pair_angles_by_seed.csv")]:
        files = [out / aggregate_name] if aggregate else sorted((out / "parts").glob(
            ("pilot" if pilot else "shard*")+f"_{kind}.csv"))
        require(files and all(p.is_file() for p in files), f"Missing {kind} inputs")
        result[kind] = pd.concat([pd.read_csv(p) for p in files], ignore_index=True)
    return result


def validate(out=OUT, pilot=False, partial=False, aggregate=False):
    frames = inputs(out, pilot, aggregate)
    p, c, z = (frames[k] for k in ["primary", "caps", "pairs"])
    require(~p.duplicated(KEY), "Duplicate primary fit keys")
    require(~c.duplicated(KEY+["rho_order"]), "Duplicate cap center keys")
    require(~z.duplicated(KEY+["rho_id", "E1_id"]), "Duplicate independent pair keys")
    expected = pd.read_csv(ROOT / "work/grid_final_geometry_by_seed.csv")
    actualkeys, expectedkeys = keyset(p), keyset(expected)
    require(actualkeys <= expectedkeys, "Unexpected primary fit key")
    require(keyset(c) == actualkeys, "Cap case set differs from primary cases")
    require(keyset(z) == actualkeys, "Pair case set differs from primary cases")
    if not pilot and not partial:
        require(actualkeys == expectedkeys, "Full valid-fit coverage not attained")
        require(len(actualkeys) == 10796, "Expected 10,796 source fits")
        require(len(p.groupby(CELL)) == 108, "Expected all 108 dimension/regime/T cells")
    require(p.domain.eq("lambda_axis_cosine_le_0.90"), "Primary domain metadata mismatch")
    require(p.d == np.where(p.regime.eq("interior"),8,5), "Wrong model dimension")
    for col in ["rho_min", "E1_at_rho", "E1_min", "rho_at_E1", "E1_cap_090", "E1_cap10"]:
        require(np.isfinite(p[col]), f"Non-finite primary objective: {col}")
        require(p[col] >= -1e-12, f"Negative primary objective: {col}")
    for col in ["E1_at_rho", "E1_min", "E1_cap_090", "E1_cap10"]:
        require(p[col] <= 1+1e-10, f"E1 outside [0,1]: {col}")
    for col in ["rho_max_loading_cosine","E1_max_loading_cosine","cap10_max_loading_cosine"]:
        require(np.isfinite(p[col]) & (p[col] >= 0) & (p[col] <= 1.000002), f"Invalid physical auxiliary cosine: {col}")
    for col in ["E1_center", "E1_cap5_saved", "E1_cap10"]:
        require(np.isfinite(c[col]) & (c[col] >= 0) & (c[col] <= 1.000002), f"Invalid cap residual: {col}")
    ra, rc, err_r = audit_lambda(p, "rho_lambda", "rho_max_lambda_cosine")
    ea, ec, err_e = audit_lambda(p, "E_lambda", "E1_max_lambda_cosine")
    ca, cc, err_c = audit_lambda(c, "center_lambda", "center_max_axis_cosine")
    c5, cos5, err5 = audit_lambda(c, "cap5_lambda", "cap5_max_axis_cosine")
    c10, cos10, err10 = audit_lambda(c, "cap10_lambda", "cap10_max_axis_cosine")
    a5, a10 = projective_angle(ca, c5), projective_angle(ca, c10)
    angle5error = float(np.max(np.abs(a5-c.center_to_cap5_deg)))
    angle10error = float(np.max(np.abs(a10-c.center_to_cap10_deg)))
    require(a5 <= 5.0001, "5-degree cap constraint violated")
    require(a10 <= 10.0001, "10-degree cap constraint violated")
    require(angle5error < 5e-5, "Reported 5-degree lambda angle mismatch")
    require(angle10error < 5e-5, "Reported 10-degree lambda angle mismatch")
    require(c.E1_cap10 <= c.E1_cap5_saved+1e-8, "10-degree cap worse than nested 5-degree witness")
    require(c.E1_cap5_saved <= c.E1_center+1e-8, "5-degree cap worse than its center witness")
    require(p.E1_cap10 <= p.E1_cap_090+1e-8, "Primary cap nesting failed")
    require(p.E1_cap_090 <= p.E1_at_rho+1e-8, "Primary center nesting failed")
    # A local-refinement row can be infeasible. Only accepted refinements
    # carry a domain-feasible flag, and that flag is tested independently.
    refined = c[flag(c.refined_domain_feasible)]
    if len(refined):
        fr = lambdas(refined, "refined_lambda")
        fc = np.max(np.abs(fr), axis=1)/np.linalg.norm(fr, axis=1)
        require(fc <= .900002, "Accepted unconstrained refinement violates axis exclusion")
    # Every pair must come from the full Cartesian product of the independent
    # rho and E1 pools. The cap optimizer endpoints are not added to this pool.
    require(np.isfinite(z.lambda_angle_deg), "Non-finite independent-pair angle")
    require((z.lambda_angle_deg >= 0) & (z.lambda_angle_deg <= 90.00001), "Pair angle outside [0,90]")
    scores = z[["rho_at_rho", "E1_at_rho", "rho_at_E1", "E1_at_E1"]].max(axis=1)
    require(np.isfinite(z[["rho_at_rho", "E1_at_rho", "rho_at_E1", "E1_at_E1"]]), "Non-finite pair endpoint objective")
    require(np.max(np.abs(scores-z.joint_max)) < 1e-10, "Independent-pair joint score incorrect")
    zz = z.assign(_score=scores, _eligible=z.lambda_angle_deg <= 5+1e-8)
    counts = zz.groupby(KEY).size().rename("_pair_count")
    eligible = zz[zz._eligible].groupby(KEY).agg(
        _eligible_count=("_score", "size"), _Q=("_score", "min"))
    nearest = zz[flag(zz.is_rho_best)].groupby(KEY).lambda_angle_deg.min().rename("_nearest")
    cart = zz.groupby(KEY).agg(_nr=("rho_id", "nunique"), _ne=("E1_id", "nunique"))
    joined = p.set_index(KEY).join([counts, eligible, nearest, cart])
    require(joined._pair_count == joined.n_pairs, "Primary pair count differs from pair pool")
    require(joined._pair_count == joined._nr*joined._ne, "Independent pair pool is not Cartesian")
    require(joined._nr == joined.n_rho_local, "Rho pool candidate count mismatch")
    require(joined._ne == joined.n_E_local, "E1 pool candidate count mismatch")
    require(joined._eligible_count.fillna(0) == joined.n_pairs_angle_le5, "Eligible-pair count mismatch")
    qe = joined._Q.fillna(np.inf).to_numpy(float)
    qa = joined.pair_joint_max.to_numpy(float)
    require(np.isfinite(qe) == np.isfinite(qa), "Q missing-pair convention incorrect")
    finite = np.isfinite(qe)
    qerror = float(np.max(np.abs(qe[finite]-qa[finite]))) if finite.any() else 0.
    require(qerror < 1e-10, "Q differs from the minimum over all eligible independent pairs")
    require(flag(joined.pair_exists) == finite, "Pair-exists flag incorrect")
    derror = float(np.max(np.abs(joined._nearest-joined.center_to_E1_nearest_deg_090)))
    require(derror < 1e-9, "Nearest independent E1 lambda distance incorrect")
    # The first cap row is exactly the globally lowest searched rho center.
    order1 = c[c.rho_order == 1].set_index(KEY)
    check = p.set_index(KEY).join(order1, rsuffix="_cap")
    require(len(order1) == len(p), "Missing rho-order-one cap row")
    require(np.max(np.abs(check.E1_center-check.E1_at_rho)) < 1e-10, "Primary E1-at-rho mismatch")
    require(np.max(np.abs(check.rho_center-check.rho_min)) < 1e-10, "Primary rho minimum mismatch")
    require(np.max(np.abs(check.E1_cap10-check.E1_cap10_cap)) < 1e-10, "Primary cap10 mismatch")
    require(np.max(np.abs(check.E1_cap_090-check.E1_cap5_saved)) < 1e-10, "Primary cap5 mismatch")
    coverage = p.groupby(CELL).agg(n=("seed", "size"), seed_min=("seed", "min"), seed_max=("seed", "max")).reset_index()
    summary = dict(status="passed", fits=len(p), cells=len(coverage), cap_centers=len(c),
        independent_pairs=len(z), eligible_pairs=int(zz._eligible.sum()),
        missing_pairs=int((~finite).sum()), complete_source_coverage=actualkeys==expectedkeys,
        domain="max(abs(lambda))/||lambda|| <= 0.90", axis_exclusion_deg=float(np.degrees(np.arccos(.9))),
        maximum_axis_cosine=float(max(rc.max(),ec.max(),cc.max(),cos5.max(),cos10.max())),
        maximum_axis_report_error=max(err_r,err_e,err_c,err5,err10),
        maximum_angle5_report_error=angle5error, maximum_angle10_report_error=angle10error,
        maximum_Q_error=qerror, maximum_nearest_independent_angle_error=derror,
        physical_cosine_not_a_constraint=True,
        primary_rho_physical_cosine_above_old_cap=int((p.rho_max_loading_cosine>.900002).sum()),
        primary_E1_physical_cosine_above_old_cap=int((p.E1_max_loading_cosine>.900002).sum()),
        cap10_physical_cosine_above_old_cap=int((c.cap10_physical_max_cosine>.900002).sum()),
        QMLE_refitted=False, old_matrix_domain_results_not_renamed=True)
    auditdir = out / "validation"
    auditdir.mkdir(exist_ok=True)
    prefix = "pilot_" if pilot else "partial_" if partial else "full_"
    coverage.to_csv(auditdir / (prefix+"export_coverage.csv"), index=False)
    (auditdir / (prefix+"exports_validation.json")).write_text(json.dumps(summary,indent=2),encoding="utf-8")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--partial", action="store_true")
    ap.add_argument("--aggregate", action="store_true")
    args = ap.parse_args()
    print(json.dumps(validate(args.out,args.pilot,args.partial,args.aggregate),indent=2))
