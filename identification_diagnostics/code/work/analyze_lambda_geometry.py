"""Aggregate newly recomputed lambda-axis diagnostics, without fitting a model.

Inputs in --output-dir:
  primary_by_seed.csv: one newly searched record per saved QMLE fit.
  all_rho_by_seed.csv (or cap10/all_rho_by_seed.csv): all new rho centers.
The original by-seed table supplies only explicitly whitelisted quantities
that depend on the saved QMLE fit, never old geometric search results.

Primary input contract: p,q,regime,seed,T,rho_min,E1_at_rho,E1_min,
rho_at_E1,center_to_E1_nearest_deg_090,pair_joint_max,pair_exists,
E1_cap_090,rho_cap_090,E1_cap10,rho_cap10,center_to_cap10_deg,
cap10_axis_edge,cap10_max_loading_cosine,cap10_max_lambda_cosine,
local_candidate_exists10,E1_local_within10_found,local_best_angle10.
All-rho records additionally use rho_order. Aliases below accept the
center/cap field names used by the earlier ten-degree exporter.

Important: `cap10_max_loading_cosine` is a physical auxiliary diagnostic;
it is NOT constrained to .90. The new constraint applies to lambda only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from itertools import product
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KEY = ["p", "q", "regime", "seed", "T"]
CELL = ["p", "q", "regime", "T"]
REGIMES = ["interior", "boundary", "exterior"]
CAP_METRICS = ["E1_center", "E1_cap5_saved", "E1_cap10"]
GEOMETRY_PAIRS = [("exterior", "interior"), ("exterior", "boundary"),
                  ("boundary", "interior")]
CAP_PAIRS = [("exterior", "interior"), ("exterior", "boundary"),
             ("interior", "boundary")]
FIT_FIELDS = ["d", "noise_estimated", "qmle_converged", "nll", "projected_gradient",
              "gradient_norm", "maxit_per_start", "qmle_starts", "J_min",
              "J_relative_min", "subspace_error", "covariance_relative_error",
              "core_p", "core_q", "fitted_A_numerical_rank", "fitted_B_numerical_rank",
              "matching_to_original_mean", "matching_to_known_set_mean",
              "matching_to_known_set_min", "missing_true_direction_id", "extra_selected",
              "is_original_basis_set", "near_original_basis", "n_known_true_directions",
              "Mx_min_eigenvalue", "Mx_max_eigenvalue", "Mx_condition_number"]
CORE_NUMERIC = ["rho_min", "E1_at_rho", "E1_min", "rho_at_E1", "E1_cap_090",
                "rho_cap_090", "E1_cap10", "rho_cap10",
                "center_to_E1_nearest_deg_090", "pair_joint_max",
                "center_to_cap10_deg", "cap10_max_loading_cosine",
                "cap10_max_lambda_cosine", "E1_local_within10_found", "local_best_angle10"]
ALIASES = {"rho_min": ["rho_center"], "E1_at_rho": ["E1_center"],
           "E1_min": ["E1_min_saved"], "rho_at_E1": ["rho_at_E1_saved"],
           "E1_cap_090": ["E1_cap5", "E1_cap5_saved"],
           "rho_cap_090": ["rho_cap5", "rho_cap5_saved"],
           "center_to_cap10_deg": ["center_to_cap_deg_10"],
           "center_to_E1_nearest_deg_090": ["center_to_E1_nearest_deg", "center_to_nearest_saved_E1_deg", "Delta_T"],
           "rho_max_lambda_cosine": ["center_max_axis_cosine"],
           "rho_max_loading_cosine": ["center_physical_max_cosine"],
           "cap10_max_lambda_cosine": ["cap10_max_axis_cosine"],
           "cap10_axis_edge": ["cap10_axis_boundary"],
           "local_candidate_exists10": ["local_candidate_exists10_found"]}
AXIS_FIELDS = ["rho_max_lambda_cosine", "E1_max_lambda_cosine",
               "cap5_max_lambda_cosine", "cap10_max_lambda_cosine"]
PHYSICAL_FIELDS = ["rho_max_loading_cosine", "E1_max_loading_cosine",
                   "cap5_max_loading_cosine", "cap10_max_loading_cosine"]
BOOLEAN_FIELDS = ["pair_exists", "qmle_converged", "cap10_axis_edge", "cap10_angle_edge",
                  "cap5_angle_edge", "cap5_axis_edge",
                  "local_candidate_exists10", "cap10_selected_converged",
                  "cap10_fallback_selected", "refined_success", "refined_domain_feasible",
                  "refined_inside10", "saved_E1_within10_exists", "rho_best_success",
                  "E_best_success", "pair_joint_solvers_success"]
TOL = 2e-6
ANGLE_TOL = 1e-4


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    return None if x is pd.NA else x


def bools(series):
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)
    result = series.astype(str).str.strip().str.upper().map(
        {"TRUE": True, "FALSE": False, "T": True, "F": False,
         "1": True, "0": False, "1.0": True, "0.0": False})
    if result.isna().any():
        raise ValueError(f"Invalid/missing Boolean values: {series.name}")
    return result.astype(bool)


def finite_stat(values):
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    return dict(mean=float(a.mean()) if len(a) else np.nan,
                variance=float(a.var(ddof=1)) if len(a) > 1 else np.nan,
                n_finite=len(a), min=float(a.min()) if len(a) else np.nan,
                max=float(a.max()) if len(a) else np.nan,
                median=float(np.median(a)) if len(a) else np.nan)


def canonicalize(frame, primary_required=True):
    frame = frame.copy()
    for target, alternatives in ALIASES.items():
        if target not in frame:
            for old in alternatives:
                if old in frame:
                    frame[target] = frame[old]
                    break
    required_numeric = CORE_NUMERIC if primary_required else [
        "rho_min", "E1_at_rho", "E1_cap_090", "rho_cap_090", "E1_cap10", "rho_cap10",
        "center_to_cap10_deg", "cap10_max_loading_cosine", "cap10_max_lambda_cosine",
        "E1_local_within10_found", "local_best_angle10"]
    required = KEY + required_numeric + ["cap10_axis_edge", "local_candidate_exists10"]
    missing = sorted(set(required) - set(frame))
    if missing:
        raise ValueError(f"New lambda-geometry export lacks fields: {missing}")
    for col in ["p", "q", "T", "seed"]:
        values = pd.to_numeric(frame[col], errors="raise")
        if values.isna().any() or not np.equal(values, values.astype(int)).all():
            raise ValueError(f"Nonintegral case key {col}")
        frame[col] = values.astype(int)
    if not frame.regime.isin(REGIMES).all() or not frame.seed.between(1, 100).all():
        raise ValueError("Unexpected regime or seed index")
    for col in set(CORE_NUMERIC + AXIS_FIELDS + PHYSICAL_FIELDS) & set(frame):
        frame[col] = pd.to_numeric(frame[col], errors="raise")
    if "pair_exists" not in frame and "pair_joint_max" in frame:
        frame["pair_exists"] = np.isfinite(frame.pair_joint_max)
    if "cap10_angle_edge" not in frame:
        frame["cap10_angle_edge"] = frame.center_to_cap10_deg >= 10 - ANGLE_TOL
    for col in set(BOOLEAN_FIELDS) & set(frame):
        frame[col] = bools(frame[col])
    if "pair_joint_max" in frame and not np.equal(frame.pair_exists, np.isfinite(frame.pair_joint_max)).all():
        raise AssertionError("Q finite/matched indicator mismatch; preserve unmatched Q=+Inf")
    unconditional = [x for x in required_numeric if x not in
                     ["pair_joint_max", "E1_local_within10_found", "local_best_angle10"]]
    if not np.isfinite(frame[unconditional].to_numpy(dtype=float)).all():
        raise ValueError("A required unconditional search metric is nonfinite")
    for col in [x for x in ["E1_at_rho", "E1_min", "E1_cap_090", "E1_cap10"] if x in frame]:
        if not frame[col].between(-TOL, 1 + TOL).all():
            raise AssertionError(f"Invalid relative rank-one error {col}")
    if (frame.E1_cap10 > frame.E1_cap_090 + 1e-8).any() or (frame.E1_cap_090 > frame.E1_at_rho + 1e-8).any():
        raise AssertionError("New-domain retained-witness nesting E10<=E5<=center failed")
    for col in set(AXIS_FIELDS) & set(frame):
        if not frame[col].between(-TOL, .90 + TOL).all():
            raise AssertionError(f"New lambda-axis feasibility failed for {col}")
    for col in set(PHYSICAL_FIELDS) & set(frame):
        if not frame[col].between(-TOL, 1 + TOL).all():
            raise AssertionError(f"A physical cosine is outside [0,1]: {col}")
    if (frame.center_to_cap10_deg > 10 + ANGLE_TOL).any():
        raise AssertionError("A cap result lies outside its ten-degree lambda cap")
    found = frame.local_candidate_exists10
    if not np.isfinite(frame.loc[found, ["E1_local_within10_found", "local_best_angle10"]].to_numpy()).all():
        raise AssertionError("Found local candidates must have finite residual and angle")
    if (frame.loc[found, "local_best_angle10"] > 10 + ANGLE_TOL).any():
        raise AssertionError("A found local candidate lies outside ten degrees")
    if frame.loc[~found, "E1_local_within10_found"].notna().any():
        raise AssertionError("Unavailable local candidate residuals must be missing, never imputed")
    frame["E1_center"] = frame.E1_at_rho
    frame["E1_cap5_saved"] = frame.E1_cap_090
    frame["rho_center"] = frame.rho_min
    frame["E1_gain_center_to10"] = frame.E1_at_rho - frame.E1_cap10
    frame["E1_gain5_to10"] = frame.E1_cap_090 - frame.E1_cap10
    frame["cap10_angle_violation_deg"] = np.maximum(0, frame.center_to_cap10_deg - 10)
    frame["cap10_axis_cosine_violation"] = np.maximum(0, frame.cap10_max_lambda_cosine - .90)
    frame["cap10_feasible_export_check"] = True
    if "center_to_E1_nearest_deg_090" in frame:
        frame["center_to_independent_E1_le5_090"] = frame.center_to_E1_nearest_deg_090 <= 5
    if "E1_min" in frame:
        frame["search_gap_vs_independent_E1"] = np.maximum(0, frame.E1_min - frame.E1_cap10)
    frame["cap10_min_axis_angle_deg"] = np.degrees(np.arccos(np.clip(frame.cap10_max_lambda_cosine, 0, 1)))
    for field in [x for x in PHYSICAL_FIELDS if x in frame]:
        point = field.replace("_max_loading_cosine", "")
        for cutoff in [.90, .95, .99]:
            frame[f"{point}_physical_cosine_gt{int(cutoff*100):03d}"] = frame[field] > cutoff
    return frame


def collect_parts(out):
    """Collect completed shard exports; pilot parts are only a fallback.

    Repeated identical rows are harmless on restart. Conflicting rows for
    one semantic key are rejected, rather than selecting a favorable fit.
    """
    part_dir = out / "parts"
    shard_primary = sorted(part_dir.glob("shard*_primary.csv"))
    prefix = "shard*" if shard_primary else "pilot"
    result = {}
    for kind, destination, keys in [
        ("primary", "primary_by_seed.csv", KEY),
        ("caps", "all_rho_by_seed.csv", KEY+["rho_order"]),
        ("pairs", "minima_pair_angles_by_seed.csv", KEY+["rho_id", "E1_id"])]:
        paths = sorted(part_dir.glob(prefix + "_" + kind + ".csv"))
        if not paths and prefix == "pilot":
            paths = [p for p in [out/("pilot_"+kind+".csv")] if p.is_file()]
        if not paths:
            if kind == "primary":
                raise FileNotFoundError("No shard or pilot primary exports are available")
            continue
        frame = pd.concat([pd.read_csv(p) for p in paths], ignore_index=True, sort=False)
        if not set(keys) <= set(frame):
            raise ValueError(f"{kind} parts lack semantic keys {set(keys)-set(frame)}")
        duplicates = frame[frame.duplicated(keys, keep=False)]
        if len(duplicates):
            nonkeys = [x for x in frame if x not in keys]
            conflicts = duplicates.groupby(keys, dropna=False)[nonkeys].nunique(dropna=False).max(axis=1)
            if (conflicts > 1).any():
                raise AssertionError(f"Conflicting duplicate {kind} records: {conflicts[conflicts>1].index.tolist()[:3]}")
            frame = frame.drop_duplicates(keys, keep="last")
        frame = frame.sort_values(keys).reset_index(drop=True)
        frame.to_csv(out / destination, index=False)
        result[kind] = dict(rows=len(frame), source_files=[str(p.relative_to(out)) for p in paths],
                            identical_duplicate_rows_removed=len(duplicates)-len(duplicates.drop_duplicates(keys)))
    (out / "collection_metadata.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def crosscheck_pairs(primary, pairs_path):
    if not pairs_path.is_file():
        return dict(pair_export_present=False)
    pairs = pd.read_csv(pairs_path)
    required = KEY + ["lambda_angle_deg", "rho_at_rho", "E1_at_rho", "rho_at_E1", "E1_at_E1"]
    if not set(required) <= set(pairs):
        raise ValueError(f"Pair export lacks {set(required)-set(pairs)}")
    pairs["_joint"] = pairs[["rho_at_rho", "E1_at_rho", "rho_at_E1", "E1_at_E1"]].max(axis=1)
    mapped = {key: g for key, g in pairs.groupby(KEY, sort=False)}
    maximum = 0.
    for row in primary.itertuples(index=False):
        key = tuple(getattr(row, c) for c in KEY)
        g = mapped.get(key)
        if g is None:
            raise AssertionError(f"No pair rows for a completed primary case: {key}")
        near = g[g.lambda_angle_deg <= 5 + 1e-8]
        expected = near._joint.min() if len(near) else np.inf
        actual = row.pair_joint_max
        if np.isinf(expected):
            if not np.isposinf(actual):
                raise AssertionError(f"Unmatched independent pair must retain Q=+Inf: {key}")
        else:
            maximum = max(maximum, abs(float(expected)-float(actual)))
            if abs(float(expected)-float(actual)) > 1e-8:
                raise AssertionError(f"Q does not agree with independent pair export: {key}")
        if bool(row.pair_exists) != bool(len(near)):
            raise AssertionError(f"Pair existence disagrees with saved pair angles: {key}")
        if "n_pairs" in primary and row.n_pairs != len(g):
            raise AssertionError(f"Pair count mismatch: {key}")
        if "n_pairs_angle_le5" in primary and row.n_pairs_angle_le5 != len(near):
            raise AssertionError(f"Near-pair count mismatch: {key}")
    return dict(pair_export_present=True, independent_pair_rows=len(pairs),
                maximum_Q_crosscheck_error=maximum)


def merge_fit_fields(frame, original):
    original = original[KEY + [x for x in FIT_FIELDS if x in original]]
    merged = frame.merge(original, on=KEY, how="left", validate="many_to_one",
                         indicator=True, suffixes=("", "__saved_fit"))
    if not merged._merge.eq("both").all():
        raise AssertionError("New search contains a key absent from the saved QMLE table")
    merged = merged.drop(columns="_merge")
    for col in FIT_FIELDS:
        saved = col + "__saved_fit"
        if saved in merged:
            if col in ["qmle_converged", "noise_estimated"]:
                same = np.equal(bools(merged[col]), bools(merged[saved])).all()
            elif pd.api.types.is_numeric_dtype(merged[col]):
                same = np.allclose(merged[col], merged[saved], atol=1e-12, rtol=1e-12, equal_nan=True)
            else:
                same = merged[col].fillna("<NA>").astype(str).eq(merged[saved].fillna("<NA>").astype(str)).all()
            if not same:
                raise AssertionError(f"A supposedly unchanged saved-fit field changed: {col}")
            merged = merged.drop(columns=saved)
    merged["qmle_converged"] = bools(merged.qmle_converged)
    return merged


def verify_keys(primary, allow_partial):
    if primary.duplicated(KEY).any():
        raise AssertionError("Repeated primary fit key")
    expected = set(product([20, 40, 80], [20, 40, 80], REGIMES, range(1, 101), [100, 200, 400, 800]))
    unavailable = {(80, 80, "interior", seed, 100) for seed in range(1, 5)}
    valid = expected - unavailable
    actual = set(primary[KEY].itertuples(index=False, name=None))
    if actual - valid:
        raise AssertionError(f"Unexpected case keys: {list(actual-valid)[:5]}")
    if not allow_partial and actual != valid:
        raise AssertionError(f"Full export requires 10796 fits; found {len(actual)}")
    return dict(primary_fits=len(actual), expected_valid_fits=len(valid),
                complete=actual == valid, allow_partial=allow_partial,
                missing_saved_fits=len(valid-actual),
                known_original_missing=[dict(zip(KEY, x)) for x in sorted(unavailable)])


def summarize(frame, metrics, flags, groups=CELL, subsets=False):
    rows = []
    for subset in (["all_valid", "qmle_converged"] if subsets else [None]):
        selected = frame if subset != "qmle_converged" else frame[frame.qmle_converged]
        for key, g in selected.groupby(groups, sort=True):
            row = dict(zip(groups, key if isinstance(key, tuple) else (key,)))
            row.update(n=len(g), n_saved=len(g), n_seeds=g.seed.nunique())
            if subset:
                row["subset"] = subset
            for metric in metrics:
                s = finite_stat(g[metric])
                row.update({metric + "_" + k: v for k, v in s.items()})
                row[metric + "_n"] = s["n_finite"]
            for flag in flags:
                row[flag + "_count"] = int(g[flag].sum())
                row[flag + "_rate"] = float(g[flag].mean())
            rows.append(row)
    return pd.DataFrame(rows)


def point_auc(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if not len(x) or not len(y) or np.isnan(x).any() or np.isnan(y).any():
        return np.nan
    return float(np.mean((x[:, None] < y[None, :]) + .5 * (x[:, None] == y[None, :])))


def geometry_auc(frame, metrics):
    rows = []
    for subset in ["all_valid", "qmle_converged"]:
        selected = frame if subset == "all_valid" else frame[frame.qmle_converged]
        for key, g in selected.groupby(["p", "q", "T"], sort=True):
            for first, second in GEOMETRY_PAIRS:
                a, b = g[g.regime == first], g[g.regime == second]
                if not len(a) or not len(b):
                    continue
                for metric in metrics:
                    rows.append(dict(zip(["p", "q", "T"], key)) |
                                dict(subset=subset, first=first, second=second, metric=metric,
                                     n_first=len(a), n_second=len(b),
                                     auc_first_smaller=point_auc(a[metric], b[metric])))
    return pd.DataFrame(rows)


def cap_auc(frame, reps, seed, allow_partial):
    weights = np.random.default_rng(seed).multinomial(100, np.full(100, .01), size=reps).astype(float)
    long, changes, wide = [], [], []
    for key, g in frame.groupby(["p", "q", "T"], sort=True):
        for first, second in CAP_PAIRS:
            a = g[g.regime == first].sort_values("seed")
            b = g[g.regime == second].sort_values("seed")
            if not len(a) or not len(b):
                continue
            wa, wb = weights[:, a.seed.to_numpy()-1], weights[:, b.seed.to_numpy()-1]
            denom = wa.sum(axis=1) * wb.sum(axis=1)
            valid = denom > 0
            if reps and not allow_partial and not valid.all():
                raise AssertionError("A full-grid bootstrap draw has no available observations")
            base = dict(zip(["p", "q", "T"], key)) | dict(first=first, second=second,
                         n_first=len(a), n_second=len(b), bootstrap_reps=reps,
                         bootstrap_valid_draws=int(valid.sum()))
            wr, points, draws = dict(base), {}, {}
            for metric in CAP_METRICS:
                x, y = a[metric].to_numpy(), b[metric].to_numpy()
                wins = (x[:, None] < y[None, :]).astype(float) + .5 * (x[:, None] == y[None, :])
                point = float(wins.mean())
                boot = np.sum((wa[valid] @ wins) * wb[valid], axis=1) / denom[valid]
                low, high = np.quantile(boot, [.025, .975]) if len(boot) else (np.nan, np.nan)
                points[metric], draws[metric] = point, boot
                long.append(base | dict(metric=metric, auc_first_smaller=point,
                                        auc_ci_low=low, auc_ci_high=high))
                wr.update({metric + "_auc": point, metric + "_ci_low": low, metric + "_ci_high": high})
            for label, ref in [("cap10_minus_center", CAP_METRICS[0]), ("cap10_minus_cap5", CAP_METRICS[1])]:
                delta = points[CAP_METRICS[2]] - points[ref]
                boot = draws[CAP_METRICS[2]] - draws[ref]
                low, high = np.quantile(boot, [.025, .975]) if len(boot) else (np.nan, np.nan)
                changes.append(base | dict(comparison=label, delta_auc=delta,
                                           delta_auc_ci_low=low, delta_auc_ci_high=high))
                wr.update({label + "_delta_auc": delta, label + "_ci_low": low, label + "_ci_high": high})
            wide.append(wr)
    return pd.DataFrame(long), pd.DataFrame(changes), pd.DataFrame(wide)


def physical_checks(frame):
    rows = []
    for key, g in frame.groupby(CELL, sort=True):
        for physical in [x for x in PHYSICAL_FIELDS if x in frame]:
            s = finite_stat(g[physical])
            row = dict(zip(CELL, key)) | dict(point=physical.replace("_max_loading_cosine", ""),
                                               n_available=len(g)) | s
            for cutoff in [.90, .95, .99]:
                valid = g[physical].notna()
                row[f"physical_cosine_gt_{int(cutoff*100):03d}_count"] = int((g.loc[valid, physical] > cutoff).sum())
                row[f"physical_cosine_gt_{int(cutoff*100):03d}_rate"] = float((g.loc[valid, physical] > cutoff).mean()) if valid.any() else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def correlations(frame):
    rows = []
    for key, g in frame.groupby(CELL, sort=True):
        row = dict(zip(CELL, key)) | dict(n=len(g))
        for label, x, y in [("center_cap5", "E1_center", "E1_cap5_saved"),
                            ("center_cap10", "E1_center", "E1_cap10"),
                            ("cap5_cap10", "E1_cap5_saved", "E1_cap10")]:
            xy = g[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()
            okay = len(xy) >= 3 and xy[x].nunique() > 1 and xy[y].nunique() > 1
            row[label + "_pearson"] = float(xy[x].corr(xy[y])) if okay else np.nan
            row[label + "_spearman"] = float(xy[x].rank().corr(xy[y].rank())) if okay else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def scientific_findings(primary, key, auc, original, out, complete):
    """Describe the actual new results, including changed old-domain conclusions."""
    fields = ["rho_min", "E1_at_rho", "E1_min", "rho_at_E1", "E1_cap5_saved", "E1_cap10",
              "rho_cap10", "J_min", "J_relative_min", "pair_joint_max",
              "center_to_E1_nearest_deg_090", "subspace_error", "covariance_relative_error"]
    matched_old = original.merge(primary[KEY], on=KEY, validate="one_to_one").copy()
    matched_old["E1_cap5_saved"] = matched_old.E1_cap_090
    matched_old["E1_center"] = matched_old.E1_at_rho
    old_cap_path = ROOT / "outputs/cap10_supplement_20261030/best_rho_by_seed.csv"
    if old_cap_path.is_file():
        cap = pd.read_csv(old_cap_path, usecols=lambda c: c in KEY+[
            "E1_cap10", "rho_cap10", "cap10_max_loading_cosine", "local_candidate_exists10"])
        add = [x for x in cap if x not in KEY and x not in matched_old]
        matched_old = matched_old.merge(cap[KEY+add], on=KEY, validate="one_to_one")
    common = [x for x in fields if x in matched_old]
    old_summary = summarize(matched_old, common, [])
    comparison = []
    old_lookup = old_summary.set_index(CELL)
    for z in key.itertuples(index=False):
        lookup_key = tuple(getattr(z, col) for col in CELL)
        before = old_lookup.loc[lookup_key]
        for field in common:
            comparison.append(dict(zip(CELL, lookup_key)) | dict(metric=field,
                matrix_mean=before[field+"_mean"], matrix_variance=before[field+"_variance"],
                matrix_n_finite=int(before[field+"_n_finite"]),
                lambda_mean=getattr(z, field+"_mean"), lambda_variance=getattr(z, field+"_variance"),
                lambda_n_finite=int(getattr(z, field+"_n_finite")),
                mean_change=getattr(z, field+"_mean")-before[field+"_mean"]))
    comparison = pd.DataFrame(comparison)
    comparison.to_csv(out / "new_vs_matrix_domain_mean_variance.csv", index=False)
    new_auc = auc[auc.subset == "all_valid"].copy()
    matched_old["qmle_converged"] = bools(matched_old.qmle_converged)
    old_auc = geometry_auc(matched_old, [x for x in new_auc.metric.unique() if x in matched_old])
    old_auc = old_auc[old_auc.subset == "all_valid"].copy()
    auc_keys = ["p", "q", "T", "first", "second", "metric"]
    auc_comparison = new_auc.merge(old_auc[auc_keys+["auc_first_smaller"]], on=auc_keys,
                                   how="inner", validate="one_to_one", suffixes=("_lambda", "_matrix"))
    auc_comparison["auc_change"] = auc_comparison.auc_first_smaller_lambda-auc_comparison.auc_first_smaller_matrix
    auc_comparison.to_csv(out / "new_vs_matrix_domain_auc.csv", index=False)
    finite_ranges = []
    for (metric, first, second, T), g in new_auc.groupby(["metric", "first", "second", "T"], sort=True):
        values = g.auc_first_smaller
        finite_ranges.append(dict(metric=metric, first=first, second=second, T=T,
             dimension_pairs=len(g), minimum=values.min(), maximum=values.max(),
             n_at_least_095=int((values >= .95-1e-12).sum()),
             n_below_half=int((values < .5-1e-12).sum()),
             n_exact_one=int((values == 1).sum())))
    finite_ranges = pd.DataFrame(finite_ranges)
    finite_ranges.to_csv(out / "finite_sample_auc_ranges.csv", index=False)
    mean_ranges = []
    for (regime, T), g in key.groupby(["regime", "T"], sort=True):
        for field in fields + ["E1_max_loading_cosine", "rho_max_loading_cosine", "cap10_max_loading_cosine"]:
            if field+"_mean" not in g:
                continue
            values = g[field+"_mean"]
            mean_ranges.append(dict(regime=regime, T=T, metric=field,
                                    dimension_pairs=len(g), minimum_cell_mean=values.min(),
                                    maximum_cell_mean=values.max()))
    mean_ranges = pd.DataFrame(mean_ranges)
    mean_ranges.to_csv(out / "cell_mean_ranges.csv", index=False)
    physical = []
    for (regime, T), g in primary.groupby(["regime", "T"], sort=True):
        for point, field in [("E1_min", "E1_max_loading_cosine"), ("rho_center", "rho_max_loading_cosine"),
                             ("cap10", "cap10_max_loading_cosine")]:
            if field not in g:
                continue
            row = dict(regime=regime, T=T, point=point, n_available=len(g))
            row.update(finite_stat(g[field]))
            for cutoff in [.90, .95, .99]:
                count = int((g[field] > cutoff).sum())
                row[f"physical_gt{int(cutoff*100):03d}_count"] = count
                row[f"physical_gt{int(cutoff*100):03d}_rate"] = count/len(g)
            if point == "E1_min":
                for tol in [.001, .005, .01, .03]:
                    row[f"E1_le_{tol:g}_and_physical_gt095_count"] = int(((g.E1_min <= tol) & (g[field] > .95)).sum())
            physical.append(row)
    physical = pd.DataFrame(physical)
    physical.to_csv(out / "physical_proximity_by_regime_T.csv", index=False)
    coverage = primary.groupby(["regime", "T"], sort=True).agg(
        n_available=("seed", "size"), n_finite_Q=("pair_exists", "sum"),
        n_local_within10=("local_candidate_exists10", "sum")).reset_index()
    coverage.to_csv(out / "conditional_count_coverage.csv", index=False)
    invariant = comparison[comparison.metric.isin(["J_min", "J_relative_min", "subspace_error", "covariance_relative_error"])]
    invariant_change = float(invariant.mean_change.abs().max()) if len(invariant) else np.nan
    if invariant_change > 1e-12:
        raise AssertionError("A supposedly domain-independent fit summary changed")
    findings = dict(complete=complete, source="new lambda-axis geometry on saved QMLE fits",
        matched_case_count=len(primary), matrix_domain_reference_case_count=len(matched_old),
        domain_independent_cell_mean_maximum_change=invariant_change,
        finite_auc_ranges=finite_ranges.to_dict(orient="records"),
        cell_mean_ranges=mean_ranges.to_dict(orient="records"),
        physical_proximity=physical.to_dict(orient="records"),
        conditional_counts=coverage.to_dict(orient="records"),
        largest_absolute_mean_changes=comparison.sort_values("mean_change", key=lambda s:s.abs(), ascending=False).head(20).to_dict(orient="records"),
        strongest_early_fixed_direction_contrasts=finite_ranges[finite_ranges['T'].isin([100,200]) &
                 (finite_ranges.minimum >= .95-1e-12)].to_dict(orient="records"),
        interpretation=["The lambda-axis domain and matrix-angle domain are different feasible sets.",
            "Physical loading cosine above .90 is an observed auxiliary value, not a lambda-domain violation.",
            "A small E1 with physical cosine near one does not by itself exhibit a well-separated additional matrix direction.",
            "Close independent candidates and positive residuals do not certify a common zero.",
            "All AUC directions are fixed first-score-smaller contrasts; no cell-wise orientation maximization is used.",
            "Pooled physical counts across dimensions are descriptive and do not create independent additional seeds.",
            "Finite sample trajectories do not prove asymptotic limits or QMLE consistency."])
    (out / "scientific_findings.json").write_text(json.dumps(clean(findings), indent=2), encoding="utf-8")
    lines = ["New lambda-axis geometry: actual numerical findings", f"Complete dataset: {complete}; matched saved fits: {len(primary)}.",
             "Lambda-axis cosine is constrained to .90; physical loading cosine is only reported.",
             "", "Finite-sample AUC ranges across dimension pairs (fixed first-score-smaller direction):"]
    for z in finite_ranges[finite_ranges['T'].isin([100,200])].itertuples(index=False):
        lines.append(f"T={z.T}; {z.metric}; {z.first}<{z.second}: [{z.minimum:.6f},{z.maximum:.6f}], "
                     f"AUC>=.95 in {z.n_at_least_095}/{z.dimension_pairs} pairs, AUC<.5 in {z.n_below_half}.")
    lines += ["", "Physical closeness at the independently lowest E1 (descriptive counts; shared seeds across dimensions):"]
    for z in physical[physical.point == "E1_min"].itertuples(index=False):
        lines.append(f"{z.regime},T={z.T}: max-loading-cosine mean={z.mean:.6f}; "
                     f">.90: {z.physical_gt090_count}/{z.n_available}; >.95: {z.physical_gt095_count}/{z.n_available}; "
                     f">.99: {z.physical_gt099_count}/{z.n_available}.")
    lines += ["", f"Maximum change in domain-independent J/subspace/covariance cell means: {invariant_change:.3g}."]
    lines += findings["interpretation"]
    (out / "scientific_findings.txt").write_text("\n".join(lines)+"\n", encoding="utf-8")
    return findings


def self_test():
    assert point_auc([0, np.inf], [1, np.inf]) == .625
    s = finite_stat([1, np.inf, 3])
    assert s["mean"] == 2 and s["variance"] == 2 and s["n_finite"] == 2
    test = pd.DataFrame([dict(p=20, q=20, T=100, regime=reg, seed=seed,
                             E1_center=value, E1_cap5_saved=value, E1_cap10=value)
                         for reg, offset in [("interior", .1), ("exterior", 0)]
                         for seed, value in [(1, offset+.01), (2, offset+.02)]])
    _, delta, _ = cap_auc(test, 100, 20261030, True)
    assert (delta[["delta_auc", "delta_auc_ci_low", "delta_auc_ci_high"]].to_numpy() == 0).all()
    # A physical cosine above .90 is allowed; a lambda-axis cosine above
    # .90 is not. These are distinct geometric conditions, not aliases.
    one = pd.DataFrame([dict(p=40, q=40, regime="interior", seed=1, T=100,
        rho_min=.05, E1_at_rho=.04, E1_min=.03, rho_at_E1=.06,
        center_to_E1_nearest_deg_090=6., pair_joint_max=np.inf,
        E1_cap_090=.038, rho_cap_090=.052, E1_cap10=.037, rho_cap10=.053,
        center_to_cap10_deg=6., cap10_axis_edge=False,
        cap10_max_loading_cosine=.99, cap10_max_lambda_cosine=.85,
        local_candidate_exists10=False, E1_local_within10_found=np.nan,
        local_best_angle10=np.nan)])
    okay = canonicalize(one)
    assert okay.cap10_max_loading_cosine.iloc[0] == .99
    bad = one.copy()
    bad["cap10_max_lambda_cosine"] = .91
    try:
        canonicalize(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError("Self-test failed to reject a lambda-domain violation")
    cap_only = one.drop(columns=["E1_min", "rho_at_E1", "center_to_E1_nearest_deg_090", "pair_joint_max"])
    assert len(canonicalize(cap_only, primary_required=False)) == 1
    print("Self-test passed: Q infinity ties, conditional variance, paired zero increment, "
          "physical-vs-lambda feasibility, per-center cap-only schema.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs/lambda_axis_geometry_20261030")
    parser.add_argument("--original-fit-table", type=Path, default=ROOT / "work/grid_final_geometry_by_seed.csv")
    parser.add_argument("--bootstrap-reps", type=int, default=1000)
    parser.add_argument("--bootstrap-seed", type=int, default=20261030)
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--allow-missing-all-rho", action="store_true")
    parser.add_argument("--collect-parts", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.bootstrap_reps < 0 or (not args.allow_partial and args.bootstrap_reps < 100):
        raise ValueError("Use >=100 paired replicates for complete data; the planned analysis uses 1000")
    out = args.output_dir.resolve()
    primary_path = out / "primary_by_seed.csv"
    collection = None
    if args.collect_parts or not primary_path.is_file():
        collection = collect_parts(out)
    if not primary_path.is_file():
        raise FileNotFoundError("Export primary_by_seed.csv after completing the new lambda-domain geometry")
    original = pd.read_csv(args.original_fit_table)
    primary = merge_fit_fields(canonicalize(pd.read_csv(primary_path)), original)
    validation = verify_keys(primary, args.allow_partial)
    validation.update(crosscheck_pairs(primary, out / "minima_pair_angles_by_seed.csv"))
    if collection:
        validation["collection"] = collection
    all_path = next((p for p in [out/"all_rho_by_seed.csv", out/"cap10/all_rho_by_seed.csv"] if p.is_file()), None)
    all_rho = None
    if all_path:
        all_rho = merge_fit_fields(canonicalize(pd.read_csv(all_path), primary_required=False), original)
        if "rho_order" not in all_rho:
            raise ValueError("All-rho export needs rho_order")
        all_rho["rho_order"] = pd.to_numeric(all_rho.rho_order, errors="raise").astype(int)
        if all_rho.duplicated(KEY+["rho_order"]).any() or not (all_rho.rho_order >= 1).all():
            raise AssertionError("Invalid/repeated rho candidate order")
        pkeys = set(primary[KEY].itertuples(index=False, name=None))
        if set(all_rho[KEY].itertuples(index=False, name=None)) != pkeys:
            raise AssertionError("All-rho case set differs from primary export")
        order1 = all_rho[all_rho.rho_order == 1]
        check = primary.merge(order1[KEY+CAP_METRICS+["rho_min"]], on=KEY,
                              validate="one_to_one", suffixes=("", "__order1"))
        if len(check) != len(primary):
            raise AssertionError("Every primary case must have rho order one")
        for col in CAP_METRICS+["rho_min"]:
            if not np.allclose(check[col], check[col+"__order1"], atol=1e-8, rtol=0):
                raise AssertionError(f"Primary vs order-one mismatch: {col}")
        extras = [x for x in AXIS_FIELDS+PHYSICAL_FIELDS+[
            "cap5_axis_edge", "cap5_angle_edge", "cap10_optimizer_success_count", "cap10_optimizer_runs",
            "cap10_fallback_selected", "refined_success", "refined_domain_feasible", "refined_inside10",
            "saved_E1_within10_exists"] if x in order1 and x not in primary]
        if extras:
            primary = primary.merge(order1[KEY+extras], on=KEY, validate="one_to_one")
            primary = canonicalize(primary)
        validation["all_rho_candidates"] = len(all_rho)
    elif not args.allow_missing_all_rho:
        raise FileNotFoundError("Export all_rho_by_seed.csv, or explicitly allow omission for a pilot")

    metrics = list(dict.fromkeys(CORE_NUMERIC + CAP_METRICS + ["rho_center", "J_min", "J_relative_min",
        "subspace_error", "covariance_relative_error", "E1_gain_center_to10", "E1_gain5_to10",
        "cap10_min_axis_angle_deg", "search_gap_vs_independent_E1"] +
        [x for x in AXIS_FIELDS+PHYSICAL_FIELDS+["n_rho_local", "n_E_local", "n_pairs", "n_pairs_angle_le5"] if x in primary]))
    physical_flags = [x for x in primary if "_physical_cosine_gt" in x]
    flags = list(dict.fromkeys([x for x in BOOLEAN_FIELDS + ["center_to_independent_E1_le5_090", "cap10_feasible_export_check"] if x in primary] + physical_flags))
    summary = summarize(primary, metrics, flags, subsets=True)
    key = summary[summary.subset == "all_valid"].drop(columns="subset").copy()
    # Legacy report names refer to lambda-domain recomputations in this dataset.
    key["local_candidate_exists10_rate"] = key.local_candidate_exists10_rate
    key["pair_exists_rate"] = key.pair_exists_rate
    key["E1_local_within10_found_n"] = key.E1_local_within10_found_n_finite
    cap_fields = CAP_METRICS + ["E1_gain_center_to10", "E1_gain5_to10", "rho_center", "rho_cap10",
                 "E1_local_within10_found", "local_best_angle10", "center_to_cap10_deg"]
    cap_summary = summarize(primary, cap_fields, [], subsets=False)
    quality = summarize(primary, ["E1_local_within10_found", "local_best_angle10", "center_to_cap10_deg",
                        "cap10_max_loading_cosine", "cap10_max_lambda_cosine", "search_gap_vs_independent_E1"],
                        [x for x in flags if x not in ["qmle_converged"]])
    auclist = ["rho_min", "E1_min", "E1_at_rho", "E1_cap_090", "E1_cap10", "J_min", "J_relative_min",
               "pair_joint_max", "center_to_E1_nearest_deg_090"]
    auc = geometry_auc(primary, auclist)
    boot_long, boot_delta, boot_wide = cap_auc(primary, args.bootstrap_reps, args.bootstrap_seed, args.allow_partial)
    frames = {"primary_by_seed_enriched": primary, "key_metrics_mean_var": key,
              "geometry_summary": summary, "cap10_mean_variance": cap_summary,
              "search_quality": quality, "geometry_auc": auc,
              "separation_auc": boot_long, "separation_auc_increments": boot_delta,
              "separation_auc_comparison_wide": boot_wide,
              "physical_direction_checks": physical_checks(primary),
              "paired_correlations_and_gains": correlations(primary)}
    if all_rho is not None:
        frames["all_rho_mean_variance"] = summarize(all_rho, [x for x in cap_fields if x in all_rho], [], groups=CELL+["rho_order"])
        frames["all_rho_search_quality"] = summarize(all_rho,
            ["cap10_max_loading_cosine", "cap10_max_lambda_cosine"],
            [x for x in flags if x in all_rho and x not in ["qmle_converged"]], groups=CELL+["rho_order"])
    tables = out / "tables"
    tables.mkdir(exist_ok=True)
    for name, frame in frames.items():
        frame.to_csv(out / (name + ".csv"), index=False)
        if name != "primary_by_seed_enriched":
            frame.to_csv(tables / (name + ".csv"), index=False)
    # Common names consumed by existing publication/workbook generators.
    cap_summary.to_csv(tables / "best_rho_mean_variance.csv", index=False)
    key.to_csv(tables / "key_metrics_mean_var.csv", index=False)
    metadata = validation | {
        "domain": "lambda_axis_angle", "lambda_axis_cosine_cap": .90,
        "minimum_lambda_axis_angle_degrees": float(np.degrees(np.arccos(.90))),
        "matrix_Frobenius_norm": 1, "local_cap_radii_degrees": [5, 10],
        "physical_loading_cosine_is_auxiliary_not_constraint": True,
        "axis_feasibility_fields_checked": [x for x in AXIS_FIELDS if x in primary],
        "bootstrap_reps": args.bootstrap_reps, "bootstrap_seed": args.bootstrap_seed,
        "bootstrap_unit": "100 shared seed indices with common multinomial weights across dimensions/T/regimes",
        "variance_denominator": "actual finite/found count minus one",
        "Q_infinite_values_retained_in_AUC": True,
        "no_QMLE_refit_in_this_analysis": True,
        "primary_sha256": hashlib.sha256(primary_path.read_bytes()).hexdigest(),
        "all_rho_sha256": hashlib.sha256(all_path.read_bytes()).hexdigest() if all_path else None,
        "fit_source_sha256": hashlib.sha256(args.original_fit_table.read_bytes()).hexdigest(),
        "key_metric_cells": len(key), "primary_qmle_converged_count": int(primary.qmle_converged.sum()),
    }
    findings = scientific_findings(primary, key, auc, original, out, validation["complete"])
    metadata["scientific_findings_available"] = True
    metadata["maximum_domain_independent_cell_mean_change"] = findings["domain_independent_cell_mean_maximum_change"]
    (out / "analysis_metadata.json").write_text(json.dumps(clean(metadata), indent=2), encoding="utf-8")
    print(json.dumps(clean(metadata), ensure_ascii=True), flush=True)


if __name__ == "__main__":
    main()
