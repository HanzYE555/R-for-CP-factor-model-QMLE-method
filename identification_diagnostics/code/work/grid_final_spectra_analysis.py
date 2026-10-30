from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SOURCE=Path('C:/Users/Lenovo/AppData/Local/Temp/7a0a08b7-1844-4bc8-bfe8-652f5efa9228_CP_QMLE_grid_20_40_80_100seeds.zip.228/CP_QMLE_grid_results')
OUT=ROOT/'work/grid_final_spectra'
OUT.mkdir(exist_ok=True)
K=['p','q','regime','seed','T']
G=['p','q','regime','T']

def summarize(df, groups, fields):
    rows=[]
    for key,g in df.groupby(groups,dropna=False,sort=True):
        if not isinstance(key,tuple):key=(key,)
        r=dict(zip(groups,key));r['n']=len(g)
        for c in fields:
            a=pd.to_numeric(g[c],errors='coerce').replace([np.inf,-np.inf],np.nan).dropna()
            r.update({c+'_mean':a.mean(),c+'_variance':a.var(ddof=1),c+'_min':a.min(),c+'_max':a.max(),c+'_n_finite':len(a)})
        rows.append(r)
    return pd.DataFrame(rows)

m=pd.read_csv(SOURCE/'metrics_by_seed.csv')
fields=['J_min','J_relative_min','subspace_error','covariance_relative_error',
        'matching_to_known_set_mean','matching_to_known_set_min','Mx_min_eigenvalue',
        'Mx_max_eigenvalue','Mx_condition_number','fitted_A_numerical_rank','fitted_B_numerical_rank',
        'extra_selected','near_original_basis']
summarize(m,G,fields).to_csv(OUT/'metric_spectra_summary_by_pq_T.csv',index=False)
summarize(m,['regime','T'],fields).to_csv(OUT/'metric_spectra_summary_pooled.csv',index=False)

j=pd.read_csv(SOURCE/'J_svd_by_seed.csv')
assert np.max(np.abs(j.gram_eigen_orig-j.orig**2))<1e-12
assert np.max(np.abs(j.gram_eigen_nor-j.nor**2))<1e-12
j['orig_minus_reference']=j.orig-j.reference_orig
j['nor_minus_reference']=j.nor-j.reference_nor
j['orig_abs_error']=j.orig_minus_reference.abs()
j['nor_abs_error']=j.nor_minus_reference.abs()
summarize(j,G+['ell','svd_index'],['orig','nor','gram_eigen_orig','gram_eigen_nor','reference_orig','reference_nor',
          'orig_minus_reference','nor_minus_reference','orig_abs_error','nor_abs_error']).to_csv(
          OUT/'all_J_singular_values_mean_variance.csv',index=False)
summarize(j,['regime','T','population_basis_key','assigned_true_direction_id','svd_index'],
          ['orig','nor','reference_orig','reference_nor','orig_abs_error','nor_abs_error']).to_csv(
          OUT/'J_spectrum_by_actual_basis.csv',index=False)
perfit=j.groupby(K,as_index=False).agg(reference_J_min=('reference_orig','min'),
       reference_J_nor_min=('reference_nor','min'),J_min_check=('orig','min'),
       J_nor_min_check=('nor','min'),J_spectrum_orig_abs_error=('orig_abs_error','mean'),
       J_spectrum_nor_abs_error=('nor_abs_error','mean'))
perfit=perfit.merge(m[K+['J_min','J_relative_min','population_basis_key']] if 'population_basis_key' in m.columns else m[K+['J_min','J_relative_min']],on=K,validate='one_to_one')
assert np.max(np.abs(perfit.J_min_check-perfit.J_min))<1e-12
assert np.max(np.abs(perfit.J_nor_min_check-perfit.J_relative_min))<1e-12
perfit['J_min_to_reference']=perfit.J_min/perfit.reference_J_min
perfit['J_nor_min_to_reference']=perfit.J_relative_min/perfit.reference_J_nor_min
summarize(perfit,G,['J_min','J_relative_min','reference_J_min','reference_J_nor_min',
          'J_spectrum_orig_abs_error','J_spectrum_nor_abs_error']).to_csv(OUT/'J_errors_against_actual_basis_reference.csv',index=False)
ext=perfit[perfit.regime=='exterior'].copy()
ext.to_csv(OUT/'exterior_J_reference_by_seed.csv',index=False)
summarize(ext,G,['J_min','J_relative_min','reference_J_min','reference_J_nor_min','J_min_to_reference',
          'J_nor_min_to_reference','J_spectrum_orig_abs_error','J_spectrum_nor_abs_error']).to_csv(OUT/'exterior_J_reference_summary.csv',index=False)

# Each ell's minimum singular value has svd_index=1 in ascending order.
je=j[j.svd_index==1]
summarize(je,G+['ell'],['orig','nor','reference_orig','reference_nor','orig_abs_error','nor_abs_error']).to_csv(
          OUT/'Jell_minimum_singular_values.csv',index=False)
del j,je,perfit

dm=pd.read_csv(SOURCE/'direction_matching_by_seed.csv')
dmfit=dm.drop_duplicates(K).copy()
summarize(dmfit,G,['matching_to_original_mean','matching_to_known_set_mean','matching_to_known_set_min',
          'extra_selected','near_original_basis','is_original_basis_set']).to_csv(OUT/'direction_recovery_by_pq_T.csv',index=False)
counts=dmfit[dmfit.regime=='exterior'].groupby(['p','q','T','population_basis_key','missing_true_direction_id'],dropna=False).size().rename('n').reset_index()
counts.to_csv(OUT/'exterior_basis_counts.csv',index=False)
ext=ext.merge(dmfit[K+['population_basis_key']],on=K,validate='one_to_one')
summarize(ext,['T','population_basis_key'],['J_min','J_relative_min','reference_J_min','reference_J_nor_min',
          'J_min_to_reference','J_nor_min_to_reference']).to_csv(OUT/'exterior_J_by_actual_basis_summary.csv',index=False)
del dm,dmfit,ext

known=pd.read_csv(SOURCE/'known_rankone_directions_by_seed.csv')
kextra=known[known.is_extra].copy()
summarize(kextra,G,['projection_residual','E1_projected','rho_projected','feasible_primary','max_loading_cosine']).to_csv(
          OUT/'extra_true_direction_projection_summary.csv',index=False)
summarize(kextra,['T'],['projection_residual','E1_projected','rho_projected','feasible_primary','max_loading_cosine']).to_csv(
          OUT/'extra_true_direction_projection_pooled.csv',index=False)
del known,kextra

spec=pd.read_csv(SOURCE/'all_matrix_spectra_by_seed.csv',usecols=K+['matrix_name','svd_index','orig','nor','numerical_rank','n_svd'])
spec=spec[spec.matrix_name!='J'].copy()
summarize(spec,G+['matrix_name','svd_index'],['orig','nor','numerical_rank']).to_csv(
          OUT/'Ar_Br_Z_C_all_spectra_summary.csv',index=False)
true_ranks={('interior','Ar'):4,('interior','Br'):4,('interior','Z'):8,('interior','C'):26,
            ('boundary','Ar'):3,('boundary','Br'):3,('boundary','Z'):5,('boundary','C'):9,
            ('exterior','Ar'):3,('exterior','Br'):3,('exterior','Z'):5,('exterior','C'):9}
per_matrix=[]
for key,g in spec.groupby(K+['matrix_name']):
    p,q,regime,seed,T,name=key
    r=true_ranks[(regime,name)]
    g=g.sort_values('svd_index');tail=g.iloc[:len(g)-r]
    per_matrix.append(dict(zip(K+['matrix_name'],key))|dict(
        population_rank=r,estimated_numerical_rank=int(g.numerical_rank.iloc[0]),
        min_orig=g.orig.min(),min_nor=g.nor.min(),largest_orig=g.orig.max(),
        tail_singular_norm=np.sqrt(np.square(tail.orig).sum()),
        tail_squared_fraction=np.square(tail.orig).sum()/np.square(g.orig).sum()))
per_matrix=pd.DataFrame(per_matrix)
per_matrix.to_csv(OUT/'matrix_rank_and_tail_by_seed.csv',index=False)
summarize(per_matrix,G+['matrix_name'],['estimated_numerical_rank','min_orig','min_nor','largest_orig',
          'tail_singular_norm','tail_squared_fraction']).to_csv(OUT/'matrix_rank_and_tail_summary.csv',index=False)
summarize(per_matrix,['regime','T','matrix_name'],['estimated_numerical_rank','min_orig','min_nor',
          'tail_singular_norm','tail_squared_fraction']).to_csv(OUT/'matrix_rank_and_tail_pooled.csv',index=False)
del spec,per_matrix

fit=pd.read_csv(SOURCE/'fit_covariance_spectra_by_seed.csv')
summarize(fit,G+['spectrum_index'],['Mx_eigenvalue','signal_singular_value','signal_singular_value_nor']).to_csv(
          OUT/'Mx_eigenvalues_and_signal_svd_summary.csv',index=False)
fp=fit.groupby(K,as_index=False).agg(Mx_min=('Mx_eigenvalue','min'),Mx_max=('Mx_eigenvalue','max'),
         signal_min_nor=('signal_singular_value_nor','min'),signal_max_nor=('signal_singular_value_nor','max'))
summarize(fp,G,['Mx_min','Mx_max','signal_min_nor','signal_max_nor']).to_csv(OUT/'Mx_and_signal_extremes_summary.csv',index=False)
summarize(fp,['regime','T'],['Mx_min','Mx_max','signal_min_nor','signal_max_nor']).to_csv(OUT/'Mx_and_signal_extremes_pooled.csv',index=False)

summary={}
summary['by_regime_T']=m.groupby(['regime','T'])[['J_min','J_relative_min','matching_to_known_set_mean','matching_to_known_set_min','Mx_min_eigenvalue','subspace_error']].mean().reset_index().to_dict('records')
summary['exterior_extra_selected_by_T']=m[m.regime=='exterior'].groupby('T').extra_selected.agg(['sum','mean']).reset_index().to_dict('records')
summary['J_min_ranges']=m.groupby(['regime','T'])[['J_min','J_relative_min']].agg(['min','max']).to_json()
(OUT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print('Wrote',len(list(OUT.glob('*.csv'))),'CSV tables to',OUT)
print('T800 J ranges:')
print(m[m['T']==800].groupby('regime')[['J_min','J_relative_min']].agg(['min','max']).to_string())
print('T800 per pq metrics:')
print(m[m['T']==800].groupby(['regime','p','q'])[['J_min','J_relative_min','subspace_error','matching_to_known_set_min']].mean().to_string())
