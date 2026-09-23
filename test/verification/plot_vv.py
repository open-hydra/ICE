"""V&V figures for the code-verification cases.

Unlike the Berthon and Doisneau scripts, this one does not read the output of a ctest
run: the verification cases build their inputs on the fly and delete them again, so
there is nothing left behind to plot. It re-runs what each figure needs instead, which
takes a few minutes, mostly in the vortex.

    cd test/verification && python3 plot_vv.py [figure ...]

Writes transparent SVGs into docs/vv/images/ and prints the numbers quoted in
docs/vv/verification.md. With no arguments it draws all of them.
"""
import cmath
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
IMG_DIR = HERE.parents[1] / 'docs' / 'vv' / 'images'
WORK = HERE / 'figure-work'

# Transparent, theme-aware SVGs, as in Berthon/plot_vv.py and Doisneau/plot_vv.py
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
    "axes.edgecolor": "#808080",
    "xtick.color": "#808080",
    "ytick.color": "#808080",
    # A title 1.2x the body, the matplotlib default, is wider than a panel of a
    # three-across figure once the type is scaled up for the page.
    "axes.titlesize": "medium",
})
EXACT = dict(color='#808080', lw=2.2, zorder=2, label='exact')
ICE = dict(color='#d1495b', lw=0.0, marker='o', ms=4.5, zorder=3, label='ICE')
ICE_LINE = dict(color='#d1495b', lw=1.4, zorder=3, label='ICE')
GUIDE = dict(color='#808080', lw=0.9, ls=':', zorder=1)

# Font sizes are chosen per figure rather than once, because the docs scale an SVG down
# to the width of the content column: a wide multi-panel figure loses a third of its
# type on the way in, and 10 pt drawn at 15 inches arrives as about 7 px. `font_for`
# inverts that, so every figure on the page ends up with text the same size on screen.
PAGE_PX = 750.0          # Material's content column, near enough
TARGET_PX = 13.0         # what the reader should end up with


def font_for(width_in):
    """Point size that renders as TARGET_PX once the page has scaled the figure."""
    return TARGET_PX * max(1.0, 72.0 * width_in / PAGE_PX)


def figure(width_in, height_in, ncols=1, **kwargs):
    """subplots() with the type pre-scaled for how wide the figure will be drawn."""
    mpl.rcParams['font.size'] = font_for(width_in)
    return plt.subplots(1, ncols, figsize=(width_in, height_in), **kwargs)


def case_module(directory, name):
    """Import a run.py that lives in a directory whose name is not an identifier."""
    spec = importlib.util.spec_from_file_location(name, str(HERE / directory / 'run.py'))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def save(fig, stem):
    out = IMG_DIR / ('verification-%s.svg' % stem)
    fig.savefig(str(out), bbox_inches='tight', transparent=True)
    plt.close(fig)
    print('   wrote %s' % out.relative_to(HERE.parents[1]))


def tidy(ax, xlabel=None, ylabel=None, title=None):
    ax.grid(alpha=0.15)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, pad=5)


# ---------------------------------------------------------------------------
#  I. The vortex
# ---------------------------------------------------------------------------

STOKES = (0.01, 0.1, 1.0, 10.0)


def vortex_runs(vx, n=96):
    """One run per Stokes number, plus the exact map at the time each one reached."""
    out = []
    for St in STOKES:
        tau = St / vx.OMEGA
        case, sol = vx.vortex_case(vx.Physics(), n, tau, vx.QUARTER, 'fig-st%g' % St)
        C, Cdot = vx.stretch(tau, sol['time'])
        out.append((St, tau, case, sol, C, Cdot))
    return out


def vortex_fields(vx, runs):
    """Four density maps: where the cloud went, against where it should have gone."""
    fig, axes = figure(16, 5.0, 4, sharey=True)
    theta = np.linspace(0.0, 2.0 * math.pi, 241)
    ring = np.exp(1j * theta)
    view = 1.65

    for k, (ax, (St, tau, case, sol, C, Cdot)) in enumerate(zip(axes, runs)):
        x, y = np.array(case.xn), np.array(case.yn)
        rho = np.array(sol['var'][0]).reshape(len(y) - 1, len(x) - 1)
        # A common scale would hide the low-St clouds, which dilate least: normalise
        # each panel by its own exact peak instead, which is rho0 / |C|^2.
        peak = vx.RHO0 / abs(C) ** 2
        pc = ax.pcolormesh(x, y, rho / peak, shading='flat', cmap='OrRd',
                           vmin=0.0, vmax=1.0, rasterized=True)

        gas = vx.X0 * ring
        ax.plot(gas.real, gas.imag, label='gas orbit', **GUIDE)
        start = vx.X0 + 2.0 * vx.SIG * ring                            # cloud at t = 0
        ax.plot(start.real, start.imag, color='#808080', lw=1.0, alpha=0.6,
                label='cloud at t = 0')
        end = C * start                                                # exact image of it
        ax.plot(end.real, end.imag, color='#808080', lw=2.0, label='exact at t_end')
        ax.plot([0.0], [0.0], marker='+', color='#808080', ms=7)

        # The numbers go inside the panel: at this type size a one-line subtitle is
        # wider than a quarter of the figure and would run into its neighbour.
        lag = vx.OMEGA * sol['time'] - cmath.phase(C)
        ax.set_title('St = %g' % St, pad=8)
        ax.annotate('lag %.3f rad\ndilation %.3f' % (lag, abs(C)),
                    xy=(0.04, 0.03), xycoords='axes fraction', va='bottom',
                    color='#808080', fontsize='small')
        ax.set_aspect('equal')
        ax.set_xlim(-view, view)
        ax.set_ylim(-view, view)
        ax.set_xlabel('x [m]')
    axes[0].set_ylabel('y [m]')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize='small', ncol=3,
               loc='upper center', bbox_to_anchor=(0.5, 0.06))
    fig.colorbar(pc, ax=axes.tolist(), label='density / exact peak',
                 fraction=0.012, pad=0.01)
    save(fig, 'vortex')


def vortex_transition(vx, runs):
    """The spiral each cloud follows, and the lag and ejection against a continuum."""
    fig, axes = figure(11.5, 4.6, 2)
    colours = ('#2e6f9e', '#4f9d69', '#e09f3e', '#9e2a2b')

    def centroid(case, sol):
        mass = sum(sol['var'][0])
        return sum(r * complex(px, py) for r, (px, py)
                   in zip(sol['var'][0], case.centres)) / mass

    ax = axes[0]
    theta = np.linspace(0.0, 2.0 * math.pi, 241)
    gas = vx.X0 * np.exp(1j * theta)
    ax.plot(gas.real, gas.imag, label='gas orbit', **GUIDE)
    for (St, tau, case, sol, C, Cdot), colour in zip(runs, colours):
        ts = np.linspace(0.0, sol['time'], 160)
        path = np.array([vx.stretch(tau, t)[0] for t in ts]) * vx.X0
        ax.plot(path.real, path.imag, color=colour, lw=1.8, label='St = %g' % St)
        centre = centroid(case, sol)
        ax.plot([centre.real], [centre.imag], marker='o', ms=6, mfc='none',
                mec=colour, mew=1.6)
    ax.plot([vx.X0], [0.0], marker='s', ms=5, color='#808080')
    ax.annotate('release', xy=(vx.X0, 0.0), xytext=(vx.X0 - 0.30, -0.16),
                color='#808080', fontsize='small')
    ax.plot([0.0], [0.0], marker='+', color='#808080', ms=8)
    ax.set_aspect('equal')
    tidy(ax, 'x [m]', 'y [m]',
         'Centroid over a quarter turn')
    ax.legend(frameon=False, fontsize='small', loc='lower left')

    ax = axes[1]
    sweep = np.logspace(-2.6, 1.6, 200)
    exact = [vx.stretch(St / vx.OMEGA, vx.QUARTER)[0] for St in sweep]
    ax.semilogx(sweep, [vx.OMEGA * vx.QUARTER - cmath.phase(C) for C in exact],
                color='#808080', lw=2.2, label='phase lag [rad]')
    ax.semilogx(sweep, [abs(C) - 1.0 for C in exact],
                color='#808080', lw=2.2, ls='--', label='dilation - 1')
    for (St, tau, case, sol, C, Cdot), colour in zip(runs, colours):
        centre = centroid(case, sol)
        ax.semilogx([St], [vx.OMEGA * sol['time'] - cmath.phase(centre)],
                    marker='o', ms=6, color=colour)
        ax.semilogx([St], [abs(centre) / vx.X0 - 1.0], marker='s', ms=5.5, mfc='none',
                    mec=colour, mew=1.6)
    ax.set_ylim(0.0, 1.0)
    tidy(ax, 'Stokes number', None, 'Lag and ejection')
    ax.legend(frameon=False, fontsize='small', loc='upper left')

    fig.tight_layout()
    save(fig, 'vortex-transition')


def figure_vortex():
    vx = case_module('I-vortex-cloud', 'ice_vx')
    vx.WORK = WORK / 'vortex'
    print('I. Vortex: four Stokes numbers on a 96^2 mesh over a quarter turn')
    runs = vortex_runs(vx)
    for St, tau, case, sol, C, Cdot in runs:
        mass = sum(sol['var'][0])
        centre = sum(r * complex(px, py) for r, (px, py)
                     in zip(sol['var'][0], case.centres)) / mass
        print('   St=%-6g |C| ICE %.6f exact %.6f   lag ICE %.6f exact %.6f'
              % (St, abs(centre) / vx.X0, abs(C),
                 vx.OMEGA * sol['time'] - cmath.phase(centre),
                 vx.OMEGA * sol['time'] - cmath.phase(C)))
    vortex_fields(vx, runs)
    vortex_transition(vx, runs)


# ---------------------------------------------------------------------------
#  A and B. Drag and thermal relaxation
# ---------------------------------------------------------------------------

def figure_relaxation():
    """The two exponentials the source terms have to reproduce, sampled in time."""
    dr = case_module('A-drag-relaxation', 'ice_drag')
    th = case_module('B-thermal-relaxation', 'ice_heat')
    dr.WORK, th.WORK = WORK / 'drag', WORK / 'heat'
    ph = dr.Physics()
    print('A/B. Relaxation: tau_drag = %.4e s, tau_thermal = %.4e s'
          % (ph.tau_stokes, ph.tau_thermal))

    fig, axes = figure(11, 4.2, 2)
    fractions = (0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0)
    curve = np.linspace(0.0, 4.2, 300)

    ax = axes[0]
    ax.plot(curve, 1.0 - np.exp(-curve), **EXACT)
    pts = [(f, dr.relax_case(ph, f * ph.tau_stokes, name='fig-drag-%g' % f))
           for f in fractions]
    ax.plot([s['time'] / ph.tau_stokes for _, s in pts],
            [s['var'][1][0] / ph.ug for _, s in pts], **ICE)
    worst = max(abs(s['var'][1][0] / ph.ug - (1.0 - math.exp(-s['time'] / ph.tau_stokes)))
                for _, s in pts)
    print('   drag:    worst departure from the exponential %.2e of u_g' % worst)
    tidy(ax, 'time / relaxation time', 'particle velocity / gas velocity',
         'A. Momentum, Stokes drag')
    ax.legend(frameon=False, fontsize='small', loc='lower right')

    ax = axes[1]
    ax.plot(curve, np.exp(-curve), **EXACT)
    pts = [(f, th.thermal_case(ph, f * ph.tau_thermal, 'fig-heat-%g' % f))
           for f in fractions]
    ax.plot([s['time'] / ph.tau_thermal for _, s in pts],
            [(s['var'][4][0] - ph.Tg) / (th.TP0 - ph.Tg) for _, s in pts], **ICE)
    worst = max(abs((s['var'][4][0] - ph.Tg) / (th.TP0 - ph.Tg)
                    - math.exp(-s['time'] / ph.tau_thermal)) for _, s in pts)
    print('   thermal: worst departure from the exponential %.2e of the initial excess'
          % worst)
    tidy(ax, 'time / relaxation time', 'temperature excess / initial excess',
         'B. Energy, Nu = 2')
    ax.legend(frameon=False, fontsize='small', loc='upper right')

    fig.tight_layout()
    save(fig, 'relaxation')


# ---------------------------------------------------------------------------
#  G. A cloud released into a uniform gas
# ---------------------------------------------------------------------------

def figure_translation():
    """Initial cloud, accelerating, translating - and the velocity curve behind it."""
    tr = case_module('G-cloud-translation', 'ice_tr')
    tr.WORK = WORK / 'translation'
    ph = tr.Physics(ug=tr.UG)
    tau, nx = ph.tau_stokes, 200
    dx = tr.L / nx
    xc = np.array([(i + 0.5) * dx for i in range(nx)])
    print('G. Translation: tau = %.4e s' % tau)

    fig, axes = figure(15, 4.6, 3)

    for ax, (shape, label) in zip(axes[:2], ((tr.wrapped_gaussian, 'Gaussian'),
                                             (tr.top_hat, 'top hat'))):
        top = 0.0
        for k, frac in enumerate((1.0, 2.0, 3.0)):
            sol = tr.cloud_case(ph, nx, shape, frac * tau, tr.DT_REF,
                                'fig-%s-%g' % (label.replace(' ', '-'), frac))
            X = tr.displacement(sol['time'], tau)
            ref = [tr.cell_average(shape, i * dx, (i + 1) * dx, X) for i in range(nx)]
            ax.plot(xc, ref, **dict(EXACT, label='exact' if k == 0 else None))
            ax.plot(xc, sol['var'][0], label='ICE' if k == 0 else None,
                    color='#d1495b', lw=1.3, zorder=3)
            top = max(top, max(ref))
            ax.annotate('%g tau' % frac, xy=(tr.XC + X, 1.06 * max(ref)),
                        ha='center', color='#808080', fontsize='small')
        initial = [tr.cell_average(shape, i * dx, (i + 1) * dx, 0.0) for i in range(nx)]
        ax.plot(xc, initial, color='#808080', lw=1.0, ls='--', label='t = 0')
        # Headroom for the time labels, so they clear the curves and the title alike
        ax.set_ylim(-0.03 * top, 1.22 * top)
        tidy(ax, 'x [m]', 'particle bulk density [kg/m\u00b3]' if label == 'Gaussian' else None,
             'Cloud profile, %s' % label)

    ax = axes[2]
    curve = np.linspace(0.0, 3.2, 300)
    ax.plot(curve, 1.0 - np.exp(-curve), **EXACT)
    fracs = (0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
    sols = [tr.cloud_case(ph, 100, tr.wrapped_gaussian, f * tau, tr.DT_REF,
                          'fig-u-%g' % f) for f in fracs]
    ax.plot([s['time'] / tau for s in sols],
            [sum(r * u for r, u in zip(s['var'][0], s['var'][1])) / sum(s['var'][0]) / tr.UG
             for s in sols], **ICE)
    tidy(ax, 'time / relaxation time', 'cloud velocity / gas velocity',
         'The velocity it moves at')

    # One legend under the row: inside the profile panels it collides with the
    # elapsed-time labels, and all three panels draw the same two things anyway.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize='small', ncol=3,
               loc='upper center', bbox_to_anchor=(0.5, 0.06))
    fig.tight_layout()
    save(fig, 'translation')


# ---------------------------------------------------------------------------
#  H. A cloud in a straining gas
# ---------------------------------------------------------------------------

def figure_strain():
    """The cloud stretching, the velocity it lags at, and the small-St asymptote."""
    st = case_module('H-linear-strain', 'ice_st')
    st.WORK = WORK / 'strain'
    ph = st.Physics()
    tau = ph.rho_al * ph.dp ** 2 / (18.0 * ph.mu)
    nx, dx = 400, st.L / 400
    xc = np.array([(i + 0.5) * dx for i in range(nx)])
    print('H. Strain: a = %.1f 1/s, tau = %.4e s, St = %.4f'
          % (st.A_STRAIN, tau, st.A_STRAIN * tau))

    fig, axes = figure(15, 4.6, 3)
    sol = st.strain_case(ph, nx, ph.dp, 'fig-exact')
    C, Cdot = st.stretch(tau, sol['time'])

    ax = axes[0]
    initial = [st.cell_average(st.profile, i * dx, (i + 1) * dx) for i in range(nx)]
    ref = [st.cell_average(lambda x: st.profile(st.xi(x) / C - st.U0 / st.A_STRAIN) / C,
                           i * dx, (i + 1) * dx) for i in range(nx)]
    ax.plot(xc, initial, color='#808080', lw=1.0, ls='--', label='t = 0')
    ax.plot(xc, ref, **EXACT)
    ax.plot(xc, sol['var'][0], color='#d1495b', lw=1.3, zorder=3, label='ICE')
    ax.annotate('stretched %.2fx,\ndiluted the same' % C,
                xy=(0.97, 0.80), xycoords='axes fraction', ha='right',
                color='#808080', fontsize='small')
    tidy(ax, 'x [m]', 'particle bulk density [kg/m\u00b3]', 'Cloud profile at t = 0.1 s')
    ax.legend(frameon=False, fontsize='small', loc='upper left')

    # Only where there are particles: in the near-vacuum outside the cloud the
    # velocity ICE reports is not a meaningful quantity.
    ax = axes[1]
    heavy = [i for i, r in enumerate(sol['var'][0]) if r > 0.01 * max(sol['var'][0])]
    xs = xc[heavy[0]:heavy[-1] + 1]
    gas = st.U0 + st.A_STRAIN * xs
    part = np.array([Cdot / C * st.xi(x) for x in xs])
    ax.fill_between(xs, part, gas, color='#d1495b', alpha=0.10, lw=0,
                    label='slip')
    ax.plot(xs, gas, color='#808080', lw=1.0, ls='--', label='gas velocity')
    ax.plot(xs, part, **EXACT)
    ax.plot(xs, sol['var'][1][heavy[0]:heavy[-1] + 1], color='#d1495b', lw=1.3,
            zorder=3, label='ICE')
    tidy(ax, 'x [m]', 'velocity [m/s]', 'The particles lag the gas')
    ax.legend(frameon=False, fontsize='small', loc='upper left')

    ax = axes[2]
    sweep = np.logspace(-2.4, -0.4, 200)
    ratio = []
    for St in sweep:
        t_s = St / st.A_STRAIN
        C_s, Cd_s = st.stretch(t_s, st.T_END)
        ratio.append((st.A_STRAIN - Cd_s / C_s) / (st.A_STRAIN ** 2 * t_s))
    ax.semilogx(sweep, ratio, **EXACT)
    ax.semilogx(sweep, [1.0 - 2.0 * St for St in sweep], color='#808080', lw=1.0,
                ls='--', label='1 - 2 St')
    for St in (0.20, 0.10, 0.05, 0.025):
        t_s = St / st.A_STRAIN
        dp_s = math.sqrt(18.0 * ph.mu * t_s / ph.rho_al)
        sol_s = st.strain_case(ph, 200, dp_s, 'fig-st%g' % St)
        d2 = st.L / 200
        xs = [(i + 0.5) * d2 for i in range(200)]
        r, u = sol_s['var'][0], sol_s['var'][1]
        slope = (sum(a * b * st.xi(x) for a, b, x in zip(r, u, xs))
                 / sum(a * st.xi(x) ** 2 for a, x in zip(r, xs)))
        ax.semilogx([St], [(st.A_STRAIN - slope) / (st.A_STRAIN ** 2 * t_s)], **ICE)
    ax.axhline(1.0, color='#808080', lw=0.7, ls=':')
    tidy(ax, 'Stokes number', 'slip / (a\u00b2 tau)', 'The small-Stokes limit')
    handles, labels = ax.get_legend_handles_labels()
    keep = dict(zip(labels, handles))
    ax.legend(keep.values(), keep.keys(), frameon=False, fontsize='small', loc='lower left')

    fig.tight_layout()
    save(fig, 'strain')


# ---------------------------------------------------------------------------
#  C, G, H, I. Where the spatial error goes under refinement
# ---------------------------------------------------------------------------

def figure_orders():
    """Every refinement study on one pair of axes."""
    print('Convergence: L1 density error against the exact solution')
    fig, ax = figure(6.4, 5.0)
    series = []

    ad = case_module('C-advection-sine', 'ice_ad')
    ad.WORK = WORK / 'advect'
    for rec, lim, tag in (('MUSCL', 'vanleer', 'C. sine, MUSCL'),
                          ('first-order', 'none', 'C. sine, first order')):
        grids, errs = (25, 50, 100, 200), []
        for nx in grids:
            sol = ad.advect_case(ad.Physics(), nx, rec, lim, 'fig-%s-%d' % (rec, nx))
            dx = ad.L / nx
            ref = [ad.cell_average(i * dx, (i + 1) * dx, sol['time']) for i in range(nx)]
            errs.append(sum(abs(a - b) for a, b in zip(sol['var'][0], ref)) / nx / ad.AMP)
        series.append((tag, grids, errs))

    tr = case_module('G-cloud-translation', 'ice_tr2')
    tr.WORK = WORK / 'orders-translation'
    ph = tr.Physics(ug=tr.UG)
    grids, errs = (50, 100, 200), []
    for nx in grids:
        sol = tr.cloud_case(ph, nx, tr.wrapped_gaussian, 3.0 * ph.tau_stokes,
                            tr.DT_REF, 'fig-order-%d' % nx)
        dx = tr.L / nx
        X = tr.displacement(sol['time'], ph.tau_stokes)
        ref = [tr.cell_average(tr.wrapped_gaussian, i * dx, (i + 1) * dx, X)
               for i in range(nx)]
        errs.append(sum(abs(a - b) for a, b in zip(sol['var'][0], ref)) / sum(ref))
    series.append(('G. translation', grids, errs))

    st = case_module('H-linear-strain', 'ice_st2')
    st.WORK = WORK / 'orders-strain'
    ph = st.Physics()
    tau = ph.rho_al * ph.dp ** 2 / (18.0 * ph.mu)
    grids, errs = (100, 200, 400), []
    for nx in grids:
        sol = st.strain_case(ph, nx, ph.dp, 'fig-order-%d' % nx)
        dx = st.L / nx
        C, _ = st.stretch(tau, sol['time'])
        ref = [st.cell_average(lambda x: st.profile(st.xi(x) / C - st.U0 / st.A_STRAIN) / C,
                               i * dx, (i + 1) * dx) for i in range(nx)]
        errs.append(sum(abs(a - b) for a, b in zip(sol['var'][0], ref)) / sum(ref))
    series.append(('H. strain', grids, errs))

    vx = case_module('I-vortex-cloud', 'ice_vx2')
    vx.WORK = WORK / 'orders-vortex'
    grids, errs = (48, 96, 192), []
    for n in grids:
        case, sol = vx.vortex_case(vx.Physics(), n, 0.1 / vx.OMEGA, vx.QUARTER,
                                   'fig-order-%d' % n)
        C, _ = vx.stretch(0.1 / vx.OMEGA, sol['time'])
        errs.append(vx.l1_error(case, sol, C))
    series.append(('I. vortex', grids, errs))

    colours = ('#2e6f9e', '#7b9eb5', '#4f9d69', '#e09f3e', '#9e2a2b')
    for (tag, grids, errs), colour in zip(series, colours):
        ax.loglog(grids, errs, marker='o', ms=5, lw=1.6, color=colour, label=tag)
        orders = [math.log(errs[i] / errs[i + 1]) / math.log(grids[i + 1] / float(grids[i]))
                  for i in range(len(errs) - 1)]
        print('   %-22s %s   orders %s'
              % (tag, ' '.join('%.2e' % e for e in errs),
                 ' '.join('%.2f' % o for o in orders)))

    guide = np.array([25.0, 400.0])
    ax.loglog(guide, 3.0e-1 * (guide / 25.0) ** -1.0, color='#808080', lw=0.9, ls=':')
    ax.loglog(guide, 3.0e-1 * (guide / 25.0) ** -2.0, color='#808080', lw=0.9, ls='--')
    ax.annotate('slope -1', xy=(330, 3.0e-1 * (330 / 25.0) ** -1.0), color='#808080',
                fontsize='small', va='bottom')
    ax.annotate('slope -2', xy=(330, 3.0e-1 * (330 / 25.0) ** -2.0), color='#808080',
                fontsize='small', va='bottom')
    ticks = (25, 50, 100, 200, 400)
    ax.set_xticks(ticks)
    ax.set_xticks([], minor=True)
    ax.set_xticklabels([str(t) for t in ticks])
    tidy(ax, 'cells across the domain', 'L1 density error',
         'Spatial convergence, every refinement study')
    ax.legend(frameon=False, fontsize='small', loc='lower left')
    fig.tight_layout()
    save(fig, 'orders')


# ---------------------------------------------------------------------------
#  D and E. Every correlation against an independent integration
# ---------------------------------------------------------------------------

def figure_correlations():
    """How far each of the nineteen correlations sits from its RK4 reference."""
    dg = case_module('D-drag-matrix', 'ice_dg')
    ht = case_module('E-heat-matrix', 'ice_ht')
    dg.WORK, ht.WORK = WORK / 'dragmatrix', WORK / 'heatmatrix'
    print('D/E. Correlations against an independent RK4 integration')

    fig, axes = figure(12, 5.0, 2)

    ph = dg.Physics()
    errs = []
    for law in dg.LAWS:
        sol = dg.drag_case(ph, law, dg.T_END, 'fig-%s' % law)
        ref, _ = ph.reference(sol['time'], drag=law, heat='Stokes')
        errs.append(abs(sol['var'][1][0] - ref) / ph.ug)
    _lollipop(axes[0], dg.LAWS, errs, 1e-3,
              'D. Twelve drag laws', 'velocity error / gas velocity')
    print('   drag:    worst %-18s %.2e' % (dg.LAWS[int(np.argmax(errs))], max(errs)))

    ph = ht.Physics()
    errs = []
    for law in ht.LAWS:
        sol = ht.heat_case(ph, law, ht.T_END, 'fig-%s' % law)
        _, tref = ph.reference(sol['time'], drag='Stokes', heat=law, Tp0=ht.TP0)
        errs.append(abs(sol['var'][4][0] - tref) / (ht.TP0 - ph.Tg))
    _lollipop(axes[1], ht.LAWS, errs, 1e-3,
              'E. Seven Nusselt laws', 'temperature error / initial gap')
    print('   thermal: worst %-18s %.2e' % (ht.LAWS[int(np.argmax(errs))], max(errs)))

    fig.tight_layout()
    save(fig, 'correlations')


def _lollipop(ax, names, errs, tol, title, xlabel):
    y = np.arange(len(names))
    floor = 10.0 ** math.floor(math.log10(min(e for e in errs if e > 0.0)) - 0.5)
    ax.hlines(y, floor, [max(e, floor) for e in errs], color='#d1495b', lw=1.4,
              alpha=0.55)
    ax.plot([max(e, floor) for e in errs], y, 'o', ms=6, color='#d1495b')
    ax.axvline(tol, color='#808080', lw=1.4, ls='--')
    ax.annotate('tolerance', xy=(tol, len(names) - 0.6), color='#808080', fontsize='small',
                ha='right', rotation=90, va='top')
    ax.set_xscale('log')
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.set_ylim(-0.7, len(names) - 0.3)
    ax.set_xlim(floor, tol * 4.0)
    tidy(ax, xlabel, None, title)


# ---------------------------------------------------------------------------
#  F. Evaporation
# ---------------------------------------------------------------------------

def figure_evaporation():
    """The d-squared law each model obeys, and the temperature two of them release."""
    ev = case_module('F-evaporation', 'ice_ev')
    ev.WORK = WORK / 'evaporation'
    print('F. Evaporation: water at %.0f K in air at %.0f K' % (ev.TP0, ev.TG))

    fig, axes = figure(15, 4.6, 3)
    colours = ('#2e6f9e', '#4f9d69', '#e09f3e', '#9e2a2b', '#7b4b94')

    # --- the d-squared law, temperature frozen ------------------------------------
    ax = axes[0]
    times = (0.01, 0.02, 0.03, 0.04)
    for model, colour in zip(ev.MODELS, colours):
        K = ev.slope(model)
        span = np.linspace(0.0, 0.04, 50)
        # At Re = 0 the ASM film correction is inactive and its rate IS CEM's, so it
        # would sit exactly on top of it: dash it to leave both visible.
        style = dict(ls='--', lw=2.4) if model == 'ASM' else dict(lw=1.6)
        ax.plot(span, 1.0 - K * span / ev.D0 ** 2, color=colour, label=model, **style)
        pts = [ev.run('fig-d2-%s-%g' % (model, t), model, t, 2.0e-5, ev.CS_FROZEN)
               for t in times]
        ax.plot([s['time'] for s in pts],
                [ev.diameter(s['var'][0][0]) ** 2 / ev.D0 ** 2 for s in pts],
                'o', ms=5, mfc='none', mec=colour, mew=1.4)
        print('   %-8s K = %.6e m2/s' % (model, K))
    tidy(ax, 'time [s]', 'squared diameter / initial', 'The d-squared law')
    ax.legend(frameon=False, fontsize='small', loc='lower left')

    # --- the coupled system, temperature released ---------------------------------
    ax = axes[1]
    times = np.linspace(0.005, 0.05, 6)
    for model, colour in zip(('d2-law', 'ASM', 'TC'), colours):
        ref = [ev.coupled(model, t, ev.CS_L)[1] for t in np.linspace(0.0, 0.05, 40)]
        ax.plot(np.linspace(0.0, 0.05, 40), ref, color=colour, lw=1.6, label=model)
        pts = [ev.run('fig-hot-%s-%g' % (model, t), model, t, 5.0e-6, ev.CS_L)
               for t in times]
        ax.plot([s['time'] for s in pts], [s['var'][4][0] for s in pts],
                'o', ms=5, mfc='none', mec=colour, mew=1.4)
    tidy(ax, 'time [s]', 'droplet temperature [K]', 'Wet-bulb approach')
    ax.legend(frameon=False, fontsize='small', loc='lower right')

    # --- the isothermal identity ---------------------------------------------------
    ax = axes[2]
    for blowing, colour, label in (('none', '#9e2a2b', 'no blowing factor'),
                                   ('LK', '#2e6f9e', 'with blowing factor')):
        pts = [ev.run('fig-blow-%s-%g' % (blowing, t), 'd2-law', t, 5.0e-6, ev.CS_L,
                      blowing=blowing) for t in times]
        ax.plot([s['time'] for s in pts], [s['var'][4][0] - ev.TP0 for s in pts],
                marker='o', ms=5, lw=1.6, color=colour, label=label)
    ax.axhline(0.0, color='#808080', lw=1.4, ls='--')
    ax.annotate('exactly isothermal', xy=(0.97, 0.08), xycoords='axes fraction',
                ha='right', color='#808080', fontsize='small')
    tidy(ax, 'time [s]', 'temperature change [K]', 'With and without blowing')
    ax.legend(frameon=False, fontsize='small', loc='upper left')

    fig.tight_layout()
    save(fig, 'evaporation')


FIGURES = {'relaxation': figure_relaxation,
           'translation': figure_translation,
           'correlations': figure_correlations,
           'evaporation': figure_evaporation,
           'strain': figure_strain,
           'orders': figure_orders,
           'vortex': figure_vortex}


if __name__ == '__main__':
    wanted = sys.argv[1:] or list(FIGURES)
    for name in wanted:
        if name not in FIGURES:
            raise SystemExit('unknown figure %r; choose from %s'
                             % (name, ', '.join(FIGURES)))
        FIGURES[name]()
