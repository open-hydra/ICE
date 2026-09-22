"""V&V figures for the Berthon Riemann problems.

Run after `ctest --test-dir build -L 1D` (it reads each case's OUTPUT/part-field.tec) with
    cd test/Berthon && python3 plot_vv.py
Writes transparent SVGs into docs/vv/images/ and prints the errors quoted in the docs.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from berthon import QUANTITIES, read_exact, read_field, sample   # noqa: E402

# Transparent, theme-aware SVGs, as in Doisneau/plot_vv.py
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
    "axes.edgecolor": "#808080",
    "xtick.color": "#808080",
    "ytick.color": "#808080",
    "font.size": 11,
})
EXACT = dict(color='#808080', lw=2.2, label='exact')
ICE = dict(color='#d1495b', lw=1.1, marker='o', ms=2.6, markevery=20, label='ICE')

HERE = Path(__file__).resolve().parent
IMG_DIR = HERE.parents[1] / "docs" / "vv" / "images"

TITLES = {'rho': r'$\rho_p$', 'u': r'$u_p$', 'v': r'$v_p$',
          'P11': r'$P_{11}$', 'P22': r'$P_{22}$', 'det(P)': r'$\det P$'}


def figure(case):
    exact = read_exact(HERE / 'Results' / case / ('%s_exact.tec' % case))
    xc, var = read_field(HERE / case / 'OUTPUT/part-field.tec')

    fig, axes = plt.subplots(2, 3, figsize=(11, 5.5), sharex=True)
    for ax, (zone, label, get) in zip(axes.ravel(), QUANTITIES):
        xe, ve = exact[zone]
        num = [get([var[v][c] for v in range(12)]) for c in range(len(xc))]
        ref = [sample(xe, ve, x) for x in xc]
        ax.plot(xe, ve, **EXACT)
        ax.plot(xc, num, **ICE)
        ax.set_title(TITLES[label], pad=4)
        ax.grid(alpha=0.15)
        err = np.abs(np.array(num) - np.array(ref)).mean() / (max(ve) - min(ve))
        print('   %-4s %-8s %.3e' % (case, label, err))
    for ax in axes[1]:
        ax.set_xlabel('x [m]')
    axes[0][0].legend(frameon=False, fontsize=9)
    fig.tight_layout()
    out = IMG_DIR / ('berthon-%s.svg' % case.lower())
    fig.savefig(out, transparent=True)
    plt.close(fig)
    print('   wrote', out.relative_to(HERE.parents[1]))


if __name__ == '__main__':
    print('L1 error, relative to the range of each exact profile:')
    for case in ('SCS', 'RCS', 'RCR'):
        figure(case)
