"""V&V figures and metrics for the Doisneau crossing-jets cases.

Run after `ctest --test-dir build` (it reads each case's OUTPUT/part-field.tec and
OUTPUT/part-residual-history.dat) with
    cd test/Doisneau && python3 plot_vv.py
Writes transparent SVGs into docs/vv/images/ and prints the numbers quoted in the docs.
The SVGs carry fixed ids and no date, so the same outputs give the same files with the same
Matplotlib; those in docs/vv/images were written by Python 3.13 with Matplotlib 3.10.1 and
NumPy 2.2.4 (another Matplotlib writes other bytes for the same figures).
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
    "svg.hashsalt": "doisneau",                     # ids from the content, not from a random salt
    # Mid-grey frame: the docs recolour SVG text for the theme, but not lines
    "axes.edgecolor": "#808080",
    "xtick.color": "#808080",
    "ytick.color": "#808080",
})
REF = dict(color='#808080', lw=2.2)                 # reference curves, readable on both themes
RHO = 'ρ [kg/m³]'
BIG = {'font.size': 15}                           # field maps and chimera panels
SVG = dict(bbox_inches='tight', transparent=True, dpi=100, metadata={'Date': None})

HERE = Path(__file__).resolve().parent
IMG_DIR = HERE.parents[1] / "docs" / "vv" / "images"

# Inlets (face 1): y in [-0.3,-0.2] at +45 deg and y in [0.2,0.3] at -45 deg, rho = g/|u| = 0.1.
# MK's upper inlet is one cell higher, y in [0.21,0.31]: its jets are mirror images about the
# centre of row 51, y = 0.005
RHO_JET, SPEED, WIDTH = 0.1, 10.0, 0.1
LOWER, UPPER, UPPER_MK, Y_MK = (-0.3, -0.2), (0.2, 0.3), (0.21, 0.31), 0.005
FLUX_NOMINAL = 2 * RHO_JET * SPEED / np.sqrt(2) * WIDTH     # int(rho u dy) of both jets
FLOOR = 1.0000001e-25                                       # P11, P22 at the floor of the AG primitives


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


def cell_variables(path):
    """Names of the cell variables of a Tecplot file, in the order of read_zones' fields."""
    head = re.split(r'^\s*ZONE', Path(path).read_text(), flags=re.IGNORECASE | re.MULTILINE)[0]
    return head.split('=', 1)[1].split()[3:]


def ramp_integral(t, h):
    """Integral of clip(s, 0, h) over s < t."""
    t = np.maximum(t, 0.0)
    return np.where(t < h, 0.5 * t * t, h * (t - 0.5 * h))


def band_area(x, y0, y1, band, s):
    """Area of each cell [x_i, x_i+1] x [y0_j, y1_j] (y1 >= y0) in the band
    band[0] <= y + s x <= band[1], s = +-1."""
    lo, hi = (x[None, :-1], x[None, 1:]) if s > 0 else (-x[None, 1:], -x[None, :-1])     # s x over the cell
    y0, h = y0[:, None], (y1 - y0)[:, None]
    # The area below the line y + s x = c: the height clip(c - s x - y0, 0, h) integrated over the cell
    below = [ramp_integral(c - y0 - lo, h) - ramp_integral(c - y0 - hi, h) for c in band]
    return np.maximum(below[1] - below[0], 0.0)                 # no round-off below zero outside the band


def exact_free(x, y):
    """Cell average of the free-streaming solution: two 45-degree bands of density 0.1 that cross unchanged."""
    area = np.diff(y)[:, None] * np.diff(x)[None, :]
    return RHO_JET * (band_area(x, y[:-1], y[1:], LOWER, -1) + band_area(x, y[:-1], y[1:], UPPER, 1)) / area


def exact_sticky(x, y):
    """Cell average of the MK solution, a pressureless gas: each jet streams freely until it reaches
    the mirror line y = Y_MK, where the two join a single jet moving along the line at their common
    x velocity, 7.07 m/s. The column mass per unit x is therefore conserved: the merged jet carries
    0.1 kg/m3 times the width of the parts of both jets that have reached the line, 0.02 kg/m2 once
    they are merged, a cell value of 2.0 in the 0.01 m row that holds the line."""
    area = np.diff(y)[:, None] * np.diff(x)[None, :]
    y_up, y_lo = np.maximum(y, Y_MK), np.minimum(y, Y_MK)     # each row's part above and below the line
    rho = RHO_JET * (band_area(x, y_up[:-1], y_up[1:], UPPER_MK, 1) +
                     band_area(x, y_lo[:-1], y_lo[1:], LOWER, -1)) / area
    # The merged mass of each column: per jet, the width that has reached the line, clip(x - x0, 0, WIDTH)
    merged = sum(ramp_integral(x[1:] - x0, WIDTH) - ramp_integral(x[:-1] - x0, WIDTH)
                 for x0 in (UPPER_MK[0] - Y_MK, Y_MK - LOWER[1]))
    r = int(np.searchsorted(y, Y_MK)) - 1
    rho[r] += RHO_JET * merged / area[r]
    return rho


def column_flux(x, y, f):
    """int(rho u dy) through each cell column, and the column centres."""
    return (f[0] * f[1] * np.diff(y)[:, None]).sum(axis=0), 0.5 * (x[1:] + x[:-1])


def residuals(m):
    """The last density residual of a case and its range over the last 1000 iterations."""
    h = np.loadtxt(HERE / m / 'OUTPUT' / 'part-residual-history.dat')
    last = h[h[:, 0] > h[-1, 0] - 1000, 2]
    return (f'residual at {h[-1, 0]:.0f} {h[-1, 2]:.1e}'
            f' (last 1000 iterations {last.min():.1e}..{last.max():.1e})')


def where(mask, xc, yc):
    """The extent of a set of cells."""
    jj, ii = np.nonzero(mask)
    if not len(ii):
        return 'none'
    xs, ys = xc[ii], np.abs(yc[jj])
    return f'x {xs.min():.3f}..{xs.max():.3f}, |y| {ys.min():.3f}..{ys.max():.3f}'


def crossing_jets():
    models = {m: read_zones(HERE / m / 'OUTPUT' / 'part-field.tec')[0] for m in ('MK', 'IG', 'AG')}
    x, y, _ = models['IG']
    free = exact_free(x, y)
    exact = {'MK': exact_sticky(x, y), 'IG': free, 'AG': free}
    xc, yc = 0.5 * (x[1:] + x[:-1]), 0.5 * (y[1:] + y[:-1])

    print('Crossing jets (100x100), distance to free streaming (IG, AG) or to the sticky merger (MK):')
    print(f'  nominal inlet mass flux  {FLUX_NOMINAL:.5f}')
    for m, (_, _, f) in models.items():
        F, _ = column_flux(x, y, f)
        d = f[0] - exact[m]
        print(f'  {m}: max rho {f[0].max():#.3g} | column flux x=0.105 {F[10]:.5f}'
              f' ({F[10] / FLUX_NOMINAL:.1%}) | RMS distance {np.sqrt(np.mean(d * d)):.1e} | {residuals(m)}')

    # MK: the merged jet's row upstream of the meeting point, and the mirror image about its centre
    rho = models['MK'][2][0]
    r, i = int(np.searchsorted(y, Y_MK)) - 1, int(np.argmin(np.abs(xc - 0.145)))
    n = min(r, len(yc) - 1 - r)                           # rows r-n..r+n mirror about row r
    a = np.abs(rho[r - n:r + n + 1] - rho[r - n:r + n + 1][::-1]).max(axis=1)
    k = int(a.argmax())
    print(f'  MK: row {r + 1} at x={xc[i]:.3f} {rho[r, i]:.1e} (exact {exact["MK"][r, i]:.0f})'
          f' | largest mirror asymmetry {a[k]:.1e}, rows {r - n + k + 1} and {r + n - k + 1}')

    # AG: the crossing, the two bands beside it, and the populated cells whose pressure tensor is at the floor
    names = cell_variables(HERE / 'AG' / 'OUTPUT' / 'part-field.tec')
    f = models['AG'][2]
    rho, dia = f[0], free > 0.199
    band = (rho > 0.15) & ~dia
    j, i = np.unravel_index(int(np.argmax(np.where(band, rho, 0))), rho.shape)
    print(f'  AG: crossing (exact 0.2) mean {rho[dia].mean():.2f} | {band.sum()} cells above 0.15 beside it'
          + (f' ({where(band, xc, yc)}): max {rho[j, i]:#.3g} at x={xc[i]:.3f}, |y|={abs(yc[j]):.3f}, exact'
             f' {free[j, i]:.1f}; up to {np.max(rho[band & (free == 0)], initial=0):.2f} where the exact'
             f' solution is empty' if band.any() else ''))
    fl = (rho > 1e-3) & ((f[names.index('P11_p1')] <= FLOOR) | (f[names.index('P22_p1')] <= FLOOR))
    between = fl & (np.abs(yc)[:, None] < UPPER[0] - xc[None, :])
    print(f'  AG: {fl.sum()} populated cells (rho > 1e-3) with P11 or P22 at the floor'
          + (f', rho {rho[fl].min():.1e}..{rho[fl].max():.1e}, exact at most {free[fl].max():.1f}:'
             f' {between.sum()} between the jets ({where(between, xc, yc)}),'
             f' {(fl & ~between).sum()} beside the crossing ({where(fl & ~between, xc, yc)})'
             if fl.any() else ''))

    # Density fields
    with plt.rc_context(BIG):
        density_fields(x, y, exact, models)

    # Density profiles across the crossing and downstream of it
    density_profiles(x, y, xc, yc, exact, models)


def density_fields(x, y, exact, models):
    fig, axs = plt.subplots(1, 5, figsize=(22.5, 5.2), sharey=True)
    panels = [('Sticky merger (exact)', exact['MK']), ('MK', models['MK'][2][0]),
              ('Free streaming (exact)', exact['IG']), ('IG', models['IG'][2][0]), ('AG', models['AG'][2][0])]
    for ax, (title, rho) in zip(axs, panels):
        pc = ax.pcolormesh(x, y, rho, shading='flat', cmap='OrRd', vmin=0, vmax=0.2, rasterized=True)
        ax.set_title(title)
        ax.set_aspect('equal')
        ax.set_xlabel('x [m]')
    axs[0].set_ylabel('y [m]')
    fig.colorbar(pc, ax=axs.tolist(), label=RHO, fraction=0.01, pad=0.01)
    plt.savefig(IMG_DIR / 'doisneau-fields.svg', **SVG)
    plt.close(fig)


def density_profiles(x, y, xc, yc, exact, models):
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, xs in zip(axs, (0.25, 0.60)):
        i = np.argmin(np.abs(xc - xs))
        ax.plot(yc, exact['IG'][:, i], label='Free streaming (exact)', **REF)
        ax.plot(yc, exact['MK'][:, i], color='C0', ls=':', lw=1.8, label='Sticky merger (exact)')
        for m, (_, _, f) in models.items():
            ax.plot(yc, f[0][:, i], lw=1.8, label=m)
        # MK and its exact solution hold the merged jet in one row: keep it off-scale
        mk, mkex = models['MK'][2][0][:, i].max(), exact['MK'][:, i].max()
        ax.annotate(f'MK peak {mk:.2f}, exact {mkex:.2f} ↑', xy=(0.0, 0.245), xytext=(0.03, 0.225), color='C0',
                    fontsize=9)
        ax.set_ylim(0, 0.25)
        ax.set_title(f'x = {xc[i]:.3f} m')
        ax.set_xlabel('y [m]')
        ax.grid(alpha=0.3)
    axs[0].set_ylabel(RHO)
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc='upper center', ncol=5, bbox_to_anchor=(0.5, 0.0))
    plt.savefig(IMG_DIR / 'doisneau-profiles.svg', **SVG)
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
    plt.savefig(IMG_DIR / 'chimera-overset.svg', **SVG)
    plt.close(fig)


if __name__ == '__main__':
    IMG_DIR.mkdir(parents=True, exist_ok=True)
    crossing_jets()
    chimera()
