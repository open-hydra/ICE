"""V&V figures and metrics for the Doisneau crossing-jets cases.

Run after `ctest --test-dir build` (it reads each case's OUTPUT/part-field.tec) with
    cd test/Doisneau && python3 plot_vv.py
Writes transparent SVGs into docs/vv/images/ and prints the numbers quoted in the docs.
"""
import re
from pathlib import Path

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Transparent, theme-aware SVGs for the docs
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
    # Mid-grey frame: the docs recolour SVG text for the theme, but not lines
    "axes.edgecolor": "#808080",
    "xtick.color": "#808080",
    "ytick.color": "#808080",
})
REF = dict(color='#808080', lw=2.2)                 # reference curves, readable on both themes
RHO = 'ρ [kg/m³]'
BIG = {'font.size': 15}                           # field maps and chimera panels

HERE = Path(__file__).resolve().parent
IMG_DIR = HERE.parents[1] / "docs" / "vv" / "images"

# Inlets (face 1): y in [-0.3,-0.2] at +45 deg and y in [0.2,0.3] at -45 deg, rho = g/|u| = 0.1
RHO_JET, SPEED = 0.1, 10.0
FLUX_NOMINAL = 2 * RHO_JET * SPEED / np.sqrt(2) * 0.1       # int(rho u dy) of both jets


def read_zones(path):
    """[(x_nodes, y_nodes, fields[var, j, i])] for each zone of a 2-D BLOCK Tecplot file."""
    zones = re.split(r'^\s*ZONE', Path(path).read_text(), flags=re.IGNORECASE | re.MULTILINE)[1:]
    out = []
    for z in zones:
        header, body = z.split('\n', 1)
        I, J, K = (int(re.search(rf'\b{k}\s*=\s*(\d+)', header, re.IGNORECASE).group(1)) for k in 'IJK')
        v = np.array([float(t) for t in body.split() if re.fullmatch(r'[-+0-9.Ee]+', t)])
        nn, nc = I * J * K, (I - 1) * (J - 1)
        x = v[:nn].reshape(K, J, I)[0, 0, :]
        y = v[nn:2 * nn].reshape(K, J, I)[0, :, 0]
        nvar = (len(v) - 3 * nn) // nc
        f = v[3 * nn:3 * nn + nvar * nc].reshape(nvar, J - 1, I - 1)
        out.append((x, y, f))
    return out


def exact_density(x, y, sub=16):
    """Cell average of the free-streaming solution: two 45-degree bands of density 0.1."""
    xc = x[:-1, None] + (np.arange(sub) + 0.5)[None, :] / sub * np.diff(x)[:, None]
    yc = y[:-1, None] + (np.arange(sub) + 0.5)[None, :] / sub * np.diff(y)[:, None]
    X = xc.reshape(-1)[None, :]
    Y = yc.reshape(-1)[:, None]
    rho = RHO_JET * (((Y - X) >= -0.3) & ((Y - X) <= -0.2)) + RHO_JET * (((Y + X) >= 0.2) & ((Y + X) <= 0.3))
    return rho.reshape(len(y) - 1, sub, len(x) - 1, sub).mean(axis=(1, 3))


def column_flux(x, y, f):
    """int(rho u dy) through each cell column, and the column centres."""
    return (f[0] * f[1] * np.diff(y)[:, None]).sum(axis=0), 0.5 * (x[1:] + x[:-1])


def crossing_jets():
    models = {m: read_zones(HERE / m / 'OUTPUT' / 'part-field.tec')[0] for m in ('MK', 'IG', 'AG')}
    x, y, _ = models['IG']
    ex = exact_density(x, y)
    xc, yc = 0.5 * (x[1:] + x[:-1]), 0.5 * (y[1:] + y[:-1])

    print('Crossing jets (100x100):')
    print(f'  nominal inlet mass flux  {FLUX_NOMINAL:.5f}')
    for m, (_, _, f) in models.items():
        F, _ = column_flux(x, y, f)
        l1 = np.abs(f[0] - ex).sum() / ex.sum()
        print(f'  {m}: max rho {f[0].max():.4f} | column flux x=0.105 {F[10]:.5f} ({F[10] / FLUX_NOMINAL:.1%})'
              f' | rel L1 vs free streaming {l1:.3f}')

    # Density fields
    with plt.rc_context(BIG):
        density_fields(x, y, ex, models)

    # Density profiles across the crossing and downstream of it
    density_profiles(x, y, xc, yc, ex, models)


def density_fields(x, y, ex, models):
    fig, axs = plt.subplots(1, 4, figsize=(18, 5.2), sharey=True)
    panels = [('Free streaming (exact)', ex)] + [(m, f[0]) for m, (_, _, f) in models.items()]
    for ax, (title, rho) in zip(axs, panels):
        pc = ax.pcolormesh(x, y, rho, shading='flat', cmap='OrRd', vmin=0, vmax=0.2, rasterized=True)
        ax.set_title(title)
        ax.set_aspect('equal')
        ax.set_xlabel('x [m]')
    axs[0].set_ylabel('y [m]')
    fig.colorbar(pc, ax=axs.tolist(), label=RHO, fraction=0.012, pad=0.01)
    plt.savefig(IMG_DIR / 'doisneau-fields.svg', bbox_inches='tight', transparent=True, dpi=100)
    plt.close(fig)


def density_profiles(x, y, xc, yc, ex, models):
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, xs in zip(axs, (0.25, 0.60)):
        i = np.argmin(np.abs(xc - xs))
        ax.plot(yc, ex[:, i], label='Free streaming (exact)', **REF)
        for m, (_, _, f) in models.items():
            ax.plot(yc, f[0][:, i], lw=1.8, label=m)
        # MK collapses the crossing jets into a delta on the axis: keep it off-scale
        mk = models['MK'][2][0][:, i].max()
        ax.annotate(f'MK peak {mk:.2f} ↑', xy=(0.0, 0.245), xytext=(0.03, 0.225), color='C0', fontsize=9)
        ax.set_ylim(0, 0.25)
        ax.set_title(f'x = {xc[i]:.3f} m')
        ax.set_xlabel('y [m]')
        ax.grid(alpha=0.3)
    axs[0].set_ylabel(RHO)
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 0.0))
    plt.savefig(IMG_DIR / 'doisneau-profiles.svg', bbox_inches='tight', transparent=True, dpi=100)
    plt.close(fig)


def chimera():
    xr, yr, ref = read_zones(HERE / 'IG' / 'OUTPUT' / 'part-field.tec')[0]
    blocks = read_zones(HERE / 'IG-chimera' / 'OUTPUT' / 'part-field.tec')
    xrc = 0.5 * (xr[1:] + xr[:-1])
    Fref, _ = column_flux(xr, yr, ref)

    print('\nChimera overset (block 1 50x100 on [0,0.5], block 2 45x73 on [0.4,1]):')
    for b, (x, y, f) in enumerate(blocks):
        F, xc = column_flux(x, y, f)
        dev = F / np.interp(xc, xrc, Fref) - 1
        print(f'  block {b + 1}: column flux vs single block  min {dev.min():+.3%}  max {dev.max():+.3%}')

    with plt.rc_context(BIG):
        chimera_figure(blocks, xr, yr, ref, xrc, Fref)


def chimera_figure(blocks, xr, yr, ref, xrc, Fref):
    fig, axs = plt.subplots(1, 3, figsize=(20, 5.6), gridspec_kw={'width_ratios': [1, 1, 1.25]})
    # Overset density field with block outlines
    ax = axs[0]
    for b in (1, 0):                                      # block 1 drawn on top of block 2
        x, y, f = blocks[b]
        pc = ax.pcolormesh(x, y, f[0], shading='flat', cmap='OrRd', vmin=0, vmax=ref[0].max(), rasterized=True)
    for (x, y, _), c, ls in zip(blocks, ('C0', 'C2'), ('-', '--')):
        ax.plot([x[0], x[-1], x[-1], x[0], x[0]], [y[0], y[0], y[-1], y[-1], y[0]], color=c, ls=ls, lw=1.8)
    ax.set_aspect('equal')
    ax.set_title('Overset IG: block 1 (solid), block 2 (dashed)')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    fig.colorbar(pc, ax=ax, label=RHO, fraction=0.046, pad=0.02)

    # Profile in the overlap: both blocks and the single-block solution at the same x
    ax = axs[1]
    xs = 0.445
    for (x, y, f), c, ls, lab in zip(blocks, ('C0', 'C2'), ('-', '--'), ('block 1', 'block 2')):
        xc, yc = 0.5 * (x[1:] + x[:-1]), 0.5 * (y[1:] + y[:-1])
        i = np.argmin(np.abs(xc - xs))
        ax.plot(yc, f[0][:, i], color=c, ls=ls, lw=1.8, label=f'{lab}, x = {xc[i]:.3f}')
    yrc = 0.5 * (yr[1:] + yr[:-1])
    ax.plot(yrc, [np.interp(xs, xrc, row) for row in ref[0]], ls=':', label='single block', **REF)
    ax.set_title('Density across the overlap region')
    ax.set_xlabel('y [m]')
    ax.set_ylabel(RHO)
    ax.legend(frameon=False, loc='upper left')
    ax.grid(alpha=0.3)

    # Column mass flux along x
    ax = axs[2]
    ax.plot(xrc, Fref, label='single block 100x100', **REF)
    for (x, y, f), c, m, lab in zip(blocks, ('C0', 'C2'), ('o', 's'), ('block 1', 'block 2')):
        F, xc = column_flux(x, y, f)
        ax.plot(xc, F, m, color=c, ms=3.5, label=lab)
    ax.axvspan(0.4, 0.5, color='0.5', alpha=0.15, lw=0, label='overlap')
    ax.set_title('Mass flux ∫ρu dy through vertical lines')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('kg/(m s)')
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(IMG_DIR / 'chimera-overset.svg', bbox_inches='tight', transparent=True, dpi=100)
    plt.close(fig)


if __name__ == '__main__':
    IMG_DIR.mkdir(parents=True, exist_ok=True)
    crossing_jets()
    chimera()
