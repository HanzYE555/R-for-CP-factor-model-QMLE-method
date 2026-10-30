"""Plot the recomputed lambda-axis diagnostic grid and unchanged fit spectra.

This creates new figures only. Archived matrix-angle figures are preserved.
"""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, PercentFormatter
import numpy as np
import pandas as pd
from PIL import Image, ImageOps, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'outputs/lambda_axis_geometry_20261030/tables/key_metrics_mean_var.csv'
PROJECT = ROOT / 'outputs/Id_Prob_QMLE_CP_Model_JASA_S23'
OUT = PROJECT / 'images'
TEX = PROJECT / 'additions/lambda/figures.tex'
DATA = pd.read_csv(SOURCE)
DIMS = [20, 40, 80]
TS = [100, 200, 400, 800]
REGIMES = ['interior', 'boundary', 'exterior']
STYLES = {
    'interior': dict(color='#2878B5', marker='o', linestyle='-', label='Identifiable interior'),
    'boundary': dict(color='#D88B20', marker='s', linestyle='--', label='Identifiable boundary'),
    'exterior': dict(color='#168572', marker='^', linestyle='-.', label='Non-identifiable'),
}

assert len(DATA) == 108
assert not DATA.duplicated(['p', 'q', 'T', 'regime']).any()
assert DATA['n'].sum() == 10796
assert set(DATA.regime) == set(REGIMES)
assert set(DATA['T']) == set(TS)
special = DATA[DATA.n != 100]
assert len(special) == 1
assert special[['p', 'q', 'T', 'n']].iloc[0].tolist() == [80, 80, 100, 96]
assert special.regime.iloc[0] == 'interior'

SPECS = [
    dict(key='rho', column='rho_min_mean', title='Polynomial residual over the complete dimension grid',
         ylabel=r'Mean $\rho_{\mathrm{ind},\min}$', multiplier=1,
         caption=r'Lowest found polynomial ratio in the recomputed coefficient-angle experiment, over all nine ordered dimension pairs. Each panel fixes $(p,q)$ and shows the three regimes at $T=100,200,400,800$. The plotted ratio is $\rho_{\mathrm{ind}}=\rho/2$, using each independent minor once. Searches impose $\max_j|\lambda_j|/\|\lambda\|_2\le0.90$, so each direction is at least $\arccos(0.90)\approx25.84^\circ$ from every coordinate axis in the fitted unit-column basis. This is a fixed exclusion from the axes; the polynomial score alone does not detect an axis Jacobian singularity.'),
    dict(key='cap10', column='E1_cap10_mean', title='Rank-one error in the 10-degree cap over the complete dimension grid',
         ylabel=r'Mean $E_{1,10}$', multiplier=1,
         caption=r'Lowest found relative rank-one error within the $10^\circ$ projective coefficient-angle cap around the lowest $\rho_{\mathrm{ind}}$ candidate, for all nine ordered dimension pairs in the recomputed coefficient-angle experiment. The common exclusion is $\max_j|\lambda_j|/\|\lambda\|_2\le0.90$. The matrix cosine to an original loading is reported separately as an auxiliary quantity and is not constrained. Cap proximity is imposed, so a small cap residual and correspondence of independently searched local candidates have different interpretations.'),
    dict(key='Jorig', column='J_min_mean', title='Absolute axis-Jacobian margin over the complete dimension grid',
         ylabel=r'Mean $j_{\mathrm{orig}}$', multiplier=1,
         caption=r'Absolute axis-Jacobian margin $j_{\mathrm{orig}}=\min_\ell s_{\ell,1}$ over all nine ordered dimension pairs. These singular values are unchanged by replacing the matrix-angle exclusion with the coefficient-angle exclusion: they depend only on the saved fitted loadings. Singular values are unsquared and use the $\psi$-column convention of the main article. Non-identifiable fits are retained across their alternative bases; a mixed-basis mean is not interpreted as a single-basis population reference.'),
    dict(key='Jnor', column='J_relative_min_mean', title='Relative axis-Jacobian margin over the complete dimension grid',
         ylabel=r'Mean $j_{\mathrm{nor}}$', multiplier=1,
         caption=r'Relative axis-Jacobian margin $j_{\mathrm{nor}}=\min_\ell s_{\ell,1}/s_{\ell,d-1}$ over all nine ordered dimension pairs. These fitted-loading spectra are unchanged by the new coefficient-angle exclusion. Each Jacobian is normalized within each seed before taking the mean. This unsquared ratio is distinct from the $s_{\ell,2}/s_{\ell,1}$ spectral gap used in the real-data table.'),
    dict(key='angle', column='center_to_E1_nearest_deg_090_mean', title='Independent-minimum correspondence over the complete dimension grid',
         ylabel=r'Mean $\Delta_T$ (degrees)', multiplier=1,
         caption=r'Projective coefficient angle $\Delta_T$ from the lowest $\rho_{\mathrm{ind}}$ candidate to its nearest independently searched $E_1$ candidate, within the same fitted unit-column basis and the common exclusion $\max_j|\lambda_j|/\|\lambda\|_2\le0.90$. This is distinct from the automatically bounded angle to a cap solution. The full observed range of cell means is shown; proximity of positive minima alone does not identify a common zero.'),
    dict(key='localrate', column='local_candidate_exists10_rate', title='Local-candidate detection within 10 degrees over the complete dimension grid',
         ylabel=r'Found within $10^\circ$ (%)', multiplier=100,
         caption=r'Fraction of estimates with a found local $E_1$ candidate within a $10^\circ$ projective coefficient angle of the lowest $\rho_{\mathrm{ind}}$ center, under the coefficient-axis exclusion $\max_j|\lambda_j|/\|\lambda\|_2\le0.90$. An eligible candidate is a successful saved independent candidate there or a successful feasible refinement that remains there after removing only the angular-cap constraint. The fixed coordinate-axis exclusion remains active. The denominator includes every available estimate, and no missing candidate is imputed as a zero residual. Numerical success is not an analytic local-optimality certificate.'),
    dict(key='physicalcos', column='E1_max_loading_cosine_mean', title='Physical loading proximity at the independent rank-one-error minimum',
         ylabel=r'Mean maximum loading cosine', multiplier=1,
         caption=r'Mean maximum absolute Frobenius cosine between the independently selected $E_1$-minimum matrix and the fitted original rank-one loadings, in the recomputed coefficient-angle experiment. This is an auxiliary physical-space proximity measure: it is not constrained by the new coefficient-axis exclusion. The dotted line at $0.90$ is the exclusion threshold used by the archived matrix-angle experiment, shown only for comparison. The full interval from zero to one is displayed, so no cell means are hidden by a zoom. A direction far from the coefficient axes can still produce a matrix close to an original loading when the loading atoms are correlated.'),
]

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 8.8,
    'axes.labelsize': 9.0,
    'axes.titlesize': 9.3,
    'xtick.labelsize': 7.8,
    'ytick.labelsize': 8.0,
    'legend.fontsize': 8.4,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.edgecolor': '#A6AEB4',
    'axes.linewidth': .55,
    'grid.color': '#D8DEE3',
    'grid.linewidth': .45,
    'mathtext.fontset': 'dejavusans',
    'savefig.facecolor': 'white',
})

OUT.mkdir(exist_ok=True)
TEX.parent.mkdir(exist_ok=True, parents=True)
tex = ['% Complete dimension-grid figures; generated from the new coefficient-angle 108-cell table.',
       '% No simulations or geometric optimizations are performed by the plot script.']
audit = {
    'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'cells': len(DATA), 'available_estimates': int(DATA.n.sum()),
    'ordered_dimension_pairs': [[p, q] for p in DIMS for q in DIMS],
    'sample_sizes': TS,
    'regimes': REGIMES,
    'missing_record_policy': 'Available records retained; only 80x80 interior T100 has n=96; no imputation.',
    'figures': [],
}

for spec in SPECS:
    col = spec['column']
    assert np.isfinite(DATA[col]).all(), col
    vals = DATA[col].to_numpy() * spec['multiplier']
    upper = 104 if spec['key'] == 'localrate' else 1.03 if spec['key'] == 'physicalcos' else vals.max() * 1.10
    fig, axes = plt.subplots(3, 3, figsize=(7.0, 6.5), sharex=True, sharey=True)
    handles = []
    points = []
    for i, p in enumerate(DIMS):
        for j, q in enumerate(DIMS):
            ax = axes[i, j]
            for reg in REGIMES:
                d = DATA[(DATA.p == p) & (DATA.q == q) & (DATA.regime == reg)].sort_values('T')
                assert d['T'].tolist() == TS
                x = d['T'].to_numpy(dtype=float)
                y = d[col].to_numpy(dtype=float) * spec['multiplier']
                line, = ax.plot(x, y, linewidth=1.35, markersize=4.2, markeredgewidth=.5, **STYLES[reg])
                if i == 0 and j == 0:
                    handles.append(line)
                for (_, z), displayed in zip(d.iterrows(), y):
                    assert displayed == z[col] * spec['multiplier']
                    points.append({'p': int(p), 'q': int(q), 'regime': reg,
                                   'T': int(z['T']), 'n': int(z.n),
                                   'source_value': float(z[col]), 'displayed_value': float(displayed)})
            ax.set_title(fr'$p={p},\ q={q}$', pad=5)
            ax.set_xlim(55, 845)
            ax.set_ylim(0, upper)
            ax.set_xticks(TS)
            ax.set_xticklabels([str(t) for t in TS])
            if spec['key'] == 'localrate':
                ax.set_yticks([0, 25, 50, 75, 100])
            else:
                ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
            ax.grid(axis='both', alpha=.72)
            if spec['key'] == 'physicalcos':
                ax.axhline(.90, color='#747D85', linestyle=':', linewidth=.85, zorder=1)
            ax.set_axisbelow(True)
            if i == 2:
                ax.set_xlabel(r'Sample size $T$', labelpad=5)
            if j == 0:
                ax.set_ylabel(spec['ylabel'], labelpad=6)
    fig.suptitle(spec['title'], fontsize=10.4, y=.985)
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .951),
               ncol=3, frameon=False, handlelength=2.5, columnspacing=1.4)
    fig.subplots_adjust(left=.112, right=.987, bottom=.087, top=.872, hspace=.34, wspace=.16)
    f = OUT / f"lambda_grid_{spec['key']}.png"
    fig.savefig(f, dpi=300)
    plt.close(fig)
    assert len(points) == 108
    assert len({(z['p'], z['q'], z['regime'], z['T']) for z in points}) == 108
    with Image.open(f) as im:
        dims = list(im.size)
    audit['figures'].append({
        'key': spec['key'], 'source_column': col, 'source_min': float(DATA[col].min()),
        'source_max': float(DATA[col].max()), 'display_multiplier': spec['multiplier'],
        'shared_y_limits': [0, float(upper)], 'pixels': dims,
        'sha256': hashlib.sha256(f.read_bytes()).hexdigest(), 'points': points,
    })
    common = (r' Points are cell means (or detection fractions), without uncertainty bands; '
              r'mean [sample variance] and actual denominators are tabulated separately. '
              r'All cells have 100 estimates except the interior $80\times80$, $T=100$ cell, which has 96. '
              r'Shared seed indices across cells do not form 900 independent replications. '
              r'The axes are shared across panels within each figure.')
    tex.extend([
        r'\clearpage', r'\begin{figure}[p]', r'\begingroup', r'\spacingset{1.08}', r'\centering',
        r'\setkeys{Gin}{width=!,height=!}',
        r'\includegraphics[width=0.98\linewidth,height=0.75\textheight,keepaspectratio]{images/' + f.name + '}',
        r'\caption{' + spec['caption'] + common + '}',
        r'\label{s23:lambda:allfig:' + spec['key'] + '}', r'\endgroup', r'\end{figure}', '',
    ])
tex.append(r'\clearpage')
TEX.write_text('\n'.join(tex) + '\n', encoding='utf-8')
(ROOT / 'work/lambda_geometry_figures_audit.json').write_text(
    json.dumps(audit, indent=2), encoding='utf-8')

# Small contact sheet for layout review; publication PNGs remain full resolution.
thumbs = []
for spec in SPECS:
    f = OUT / f"lambda_grid_{spec['key']}.png"
    with Image.open(f).convert('RGB') as im:
        thumb = ImageOps.contain(im, (700, 650))
        thumbs.append(thumb.copy())
sheet = Image.new('RGB', (1400, 650*((len(SPECS)+1)//2)), 'white')
for k, im in enumerate(thumbs):
    sheet.paste(im, ((k % 2) * 700, (k // 2) * 650))
sheet.save(ROOT / 'work/lambda_geometry_figures_contact.png')
print(json.dumps({'figures': len(SPECS), 'cells_per_figure': 108,
                  'points_verified': sum(len(x['points']) for x in audit['figures']),
                  'available_estimates': int(DATA.n.sum()),
                  'files': [f"lambda_grid_{x['key']}.png" for x in SPECS]}, ensure_ascii=True))
