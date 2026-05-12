"""
Optics and emittance analysis functions matching plot_g4bl_emittance.py,
working directly on a BDSIM-style pandas DataFrame.

Expected DataFrame columns:
    samplerName  - string station label
    x, y, z      - position [metres]
    xp, yp, zp   - direction cosines (px/p, py/p, pz/p)
    p            - total momentum [GeV/c]
    energy       - total energy [GeV]
    T            - time of flight [ns]
    weight       - statistical weight

All physics calculations use mm and MeV/c internally (matching xboa conventions).
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

C_LIGHT = 299.792   # mm/ns
M_MU    = 105.658   # MeV/c²

# ---------------------------------------------------------------------------
# Publication style — applied at import time
# ---------------------------------------------------------------------------

plt.rcParams.update({
    'font.family':               'serif',
    'font.serif':                ['cmr10', 'DejaVu Serif'],
    'mathtext.fontset':          'cm',
    'axes.formatter.use_mathtext': True,
    'font.size':                 12,
    'axes.labelsize':            13,
    'axes.titlesize':            14,
    'axes.titleweight':          'normal',
    'axes.grid':                 True,
    'grid.alpha':                0.3,
    'grid.linestyle':            '--',
    'legend.frameon':            True,
    'legend.framealpha':         0.9,
    'legend.edgecolor':          'gray',
    'legend.fancybox':           True,
    'figure.dpi':                110,
    'axes.unicode_minus':        False,
})

_COLOR_X   = '#d62728'
_COLOR_Y   = '#1f77b4'
_COLOR_AUX = 'black'
_FIGSIZE   = (8, 5)


def _style_axes(ax):
    for spine in ax.spines.values():
        spine.set_visible(True)
    ax.margins(x=0)


# ---------------------------------------------------------------------------
# Unit helpers
# ---------------------------------------------------------------------------

def _momenta(df):
    """Return px, py, pz in MeV/c and p in MeV/c from direction cosines."""
    p  = df["p"].to_numpy()  * 1000.0
    px = df["xp"].to_numpy() * p
    py = df["yp"].to_numpy() * p
    pz = df["zp"].to_numpy() * p
    return px, py, pz, p

def _positions(df):
    """Return x, y, z in mm."""
    x = df["x"].to_numpy() * 1000.0
    y = df["y"].to_numpy() * 1000.0
    z = df["z"].to_numpy() * 1000.0
    return x, y, z

def _ct(df):
    """Return c*t in mm."""
    return df["T"].to_numpy() * C_LIGHT

def _energy(df):
    """Return total energy in MeV."""
    return df["energy"].to_numpy() * 1000.0

def _w(df):
    return df["weight"].to_numpy()


# ---------------------------------------------------------------------------
# Per-station physics
# ---------------------------------------------------------------------------

def transverse_emittance(df):
    """
    Normalised 4D transverse emittance [mm].

    Formula: det(Cov(x, px, y, py))^(1/4) / m_mu
    """
    x, y, _ = _positions(df)
    px, py, _, _ = _momenta(df)
    w = _w(df)
    cov = np.cov([x, px, y, py], aweights=w)
    return np.linalg.det(cov)**0.25 / M_MU


def longitudinal_emittance(df):
    """
    Normalised longitudinal emittance [mm].

    Formula: det(Cov(ct, E))^(1/2) / m_mu
    """
    ct = _ct(df)
    E  = _energy(df)
    w  = _w(df)
    cov = np.cov([ct, E], aweights=w)
    return np.linalg.det(cov)**0.5 / M_MU


def emittance_6d(df):
    """
    Normalised 6D emittance [mm^3].

    Formula: det(Cov(ct, E, x, px, y, py))^(1/2) / m_mu^3
    """
    ct = _ct(df)
    E  = _energy(df)
    x, y, _ = _positions(df)
    px, py, _, _ = _momenta(df)
    w = _w(df)
    cov = np.cov([ct, E, x, px, y, py], aweights=w)
    return np.linalg.det(cov)**0.5 / M_MU**3


def beta_function(df):
    """
    Geometric beta functions [mm].

    beta_x = Var(x) / eps_x,  eps_x = sqrt(det(Cov(x, xp)))
    beta_y = Var(y) / eps_y,  eps_y = sqrt(det(Cov(y, yp)))

    Returns (beta_x, beta_y) in mm.
    """
    x, y, _ = _positions(df)
    xp = df["xp"].to_numpy()
    yp = df["yp"].to_numpy()
    w  = _w(df)

    cov_x = np.cov([x, xp], aweights=w)
    cov_y = np.cov([y, yp], aweights=w)
    eps_x = np.linalg.det(cov_x)**0.5
    eps_y = np.linalg.det(cov_y)**0.5
    var_x = np.cov(x, aweights=w)
    var_y = np.cov(y, aweights=w)
    return var_x / eps_x, var_y / eps_y


def amplitude_pz_correlation(df):
    """
    Pearson correlation between transverse amplitude A = sqrt(x²+y²) and pz.
    """
    x, y, _ = _positions(df)
    _, _, pz, _ = _momenta(df)
    A = np.sqrt(x**2 + y**2)
    w = _w(df)
    cov = np.cov([A, pz], aweights=w)
    return cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])


# ---------------------------------------------------------------------------
# Cuts
# ---------------------------------------------------------------------------

def apply_cuts(df,
               r_max=None,
               s_max=None,
               transmitted=False,
               x="x", y="y", s="S",
               simName="SimName", trackID="trackID"):
    """
    Apply analysis cuts to a DataFrame.

    Parameters
    ----------
    r_max       : float or None — transverse aperture cut [same units as x, y]
    s_max       : float or None — max arc length S to include
    transmitted : bool — if True, keep only particles present at the last S
                  value within each SimName
    x, y, s     : column names for transverse positions and arc length
    simName     : column name for simulation identifier
    trackID     : column name for particle track ID

    Returns a filtered copy of df.
    """
    df = df.copy()

    if r_max is not None:
        df = df[df[x]**2 + df[y]**2 <= r_max**2]

    if s_max is not None:
        df = df[df[s] <= s_max]

    if transmitted:
        keep = []
        for sim, sub in df.groupby(simName, sort=False):
            s_end = sub[s].max()
            surviving_ids = sub.loc[sub[s] == s_end, trackID].unique()
            keep.append(sub[sub[trackID].isin(surviving_ids)])
        df = pd.concat(keep, ignore_index=True) if keep else df.iloc[0:0]

    return df


# ---------------------------------------------------------------------------
# Per-station sweep
# ---------------------------------------------------------------------------

def _stations(df):
    """Yield (samplerName, group) sorted by mean z."""
    mean_z = df.groupby("samplerName")["z"].mean()
    for name in mean_z.sort_values().index:
        yield name, df[df["samplerName"] == name]


def compute_optics(df):
    """
    Compute all optics quantities at each station.

    Returns a DataFrame with columns:
        samplerName, mean_z,
        emittance_trans, emittance_long, emittance_6d,
        beta_x, beta_y, mean_p, transmission, a_pz_corr
    """
    w0 = None
    rows = []
    for name, group in _stations(df):
        w_sum = group["weight"].sum()
        if w0 is None:
            w0 = w_sum
        _, _, _, p = _momenta(group)
        w = _w(group)
        rows.append({
            "samplerName":     name,
            "mean_z":          group["z"].mean() * 1000.0,
            "emittance_trans": transverse_emittance(group),
            "emittance_long":  longitudinal_emittance(group),
            "emittance_6d":    emittance_6d(group),
            "beta_x":          beta_function(group)[0],
            "beta_y":          beta_function(group)[1],
            "mean_p":          np.average(p, weights=w),
            "transmission":    100.0 * w_sum / w0,
            "a_pz_corr":       amplitude_pz_correlation(group),
        })
    return pd.DataFrame(rows).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_beta(optics, title="Beta function", save_path=None):
    """Beta functions beta_x and beta_y vs z."""
    z = optics["mean_z"] / 1000.0

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.plot(z, optics["beta_x"], label=r'$\beta_x$', color=_COLOR_X, linewidth=2)
    ax.plot(z, optics["beta_y"], label=r'$\beta_y$', color=_COLOR_Y, linewidth=2, linestyle='--')
    _style_axes(ax)
    ax.set_xlabel(r'$z$ (m)')
    ax.set_ylabel(r'Beta function (mm)')
    ax.set_title(title)
    ax.legend(loc='best', fontsize=12)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path)
    return fig


def plot_emittance(optics, title="Emittance", save_path=None):
    """
    Transverse and longitudinal emittance on left axis, 6D on right axis.
    Exponential fits overlaid on the last 60% of data points.
    """
    z       = optics["mean_z"].values / 1000.0
    e_trans = optics["emittance_trans"].values
    e_long  = optics["emittance_long"].values
    e_6d    = optics["emittance_6d"].values

    fig, ax1 = plt.subplots(figsize=_FIGSIZE)
    ax2 = ax1.twinx()

    l1, = ax1.plot(z, e_trans, label=r'$\varepsilon_\perp$ (transverse)',
                   color=_COLOR_X, linewidth=2)
    l2, = ax1.plot(z, e_long,  label=r'$\varepsilon_\parallel$ (longitudinal)',
                   color=_COLOR_Y, linewidth=2, linestyle='--')

    start = int(np.floor(0.40 * len(z)))
    tail  = np.zeros(len(z), dtype=bool)
    tail[start:] = True

    def _fit_exp(z_fit, y_fit):
        m = y_fit > 0
        k, logA = np.polyfit(z_fit[m], np.log(y_fit[m]), 1)
        return np.exp(logA), k

    At, kt = _fit_exp(z[tail], e_trans[tail])
    Al, kl = _fit_exp(z[tail], e_long[tail])
    ax1.plot(z, At * np.exp(kt * z), color=_COLOR_X, linewidth=1.5,
             linestyle=':', label='_nolegend_')
    ax1.plot(z, Al * np.exp(kl * z), color=_COLOR_Y, linewidth=1.5,
             linestyle=':', label='_nolegend_')

    _style_axes(ax1)
    ax1.set_xlabel(r'$z$ (m)')
    ax1.set_ylabel(r'Transverse / longitudinal emittance (mm$\,$rad)')
    ax1.set_title(title)

    l3, = ax2.plot(z, e_6d, label=r'$\varepsilon_{6\mathrm{D}}$',
                   color=_COLOR_AUX, linewidth=2, linestyle=':')
    ax2.set_ylabel(r'6D emittance (mm$^3\,$rad$^3$)', color=_COLOR_AUX)
    ax2.tick_params(axis='y', labelcolor=_COLOR_AUX)
    ax2.grid(False)
    for spine in ax2.spines.values():
        spine.set_visible(True)

    ax1.legend(handles=[l1, l2, l3], loc='best', fontsize=11)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path)
    return fig


def plot_transmission(df, r_max=None, s_max=None, n_initial=None,
                      title="Transmission", save_path=None):
    """
    Transmission vs z computed from raw (unfiltered) df.

    Groups by the 'z' column (set equal to S before calling) and counts
    particles at each station. n_initial defaults to the count at the
    first station.
    """
    df_cut = apply_cuts(df, r_max=r_max, s_max=s_max, transmitted=False)
    counts = df_cut.groupby("z").size()
    if n_initial is None:
        n_initial = counts.iloc[0]
    transmission = counts / n_initial

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.plot(transmission.index / 1000, transmission.values,
            color=_COLOR_Y, linewidth=2)
    _style_axes(ax)
    ax.set_xlabel(r'$z$ (m)')
    ax.set_ylabel(r'Transmission')
    ax.set_title(title)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path)
    return fig


def plot_momentum(optics, title="Mean momentum", save_path=None):
    """Mean momentum vs z."""
    z = optics["mean_z"] / 1000.0

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.plot(z, optics["mean_p"], color=_COLOR_Y, linewidth=2)
    _style_axes(ax)
    ax.set_xlabel(r'$z$ (m)')
    ax.set_ylabel(r'Mean momentum $\langle p \rangle$ (MeV/$c$)')
    ax.set_title(title)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path)
    return fig


def plot_amplitude_pz_corr(optics, title=r'Amplitude--$p_z$ correlation',
                            save_path=None):
    """Amplitude–pz correlation vs z."""
    z = optics["mean_z"] / 1000.0

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.plot(z, optics["a_pz_corr"], color=_COLOR_X, linewidth=2)
    _style_axes(ax)
    ax.set_xlabel(r'$z$ (m)')
    ax.set_ylabel(r'$a_{p_z}$ correlation')
    ax.set_title(title)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path)
    return fig


def plot_all(optics, df_raw, r_max=None, s_max=None, cell_length=None,
             show=True, save=None):
    """
    Render all five figures, optionally showing and/or saving them.

    Parameters
    ----------
    optics      : DataFrame returned by compute_optics()
    df_raw      : unfiltered DataFrame (used for the transmission plot)
    r_max       : aperture cut [m]
    s_max       : zoom all plots to z <= s_max [m]; None shows full range
    cell_length : if given, beta/momentum/a_pz_corr plots show only stations
                  at multiples of cell_length [m] (e.g. 2.0 for 2 m cells)
    show        : bool — display figures interactively
    save        : str or None — directory to save figures into; None skips saving
    """
    import os

    # Zoom optics to s_max [m] → mean_z is stored in mm
    optics_cut = (optics[optics["mean_z"] <= s_max * 1000]
                  if s_max is not None else optics)

    # Cell-boundary filter: keep only stations at multiples of cell_length
    if cell_length is not None:
        cl_mm = round(cell_length * 1000)
        optics_abs = optics_cut[
            optics_cut["mean_z"].round().astype(int) % cl_mm == 0
        ]
    else:
        optics_abs = optics_cut

    # Transmission plot uses apply_cuts internally which expects mm
    s_max_mm = s_max * 1000 if s_max is not None else None

    figs = {
        "beta.png":         plot_beta(optics_abs),
        "emittance.png":    plot_emittance(optics_cut),
        "transmission.png": plot_transmission(df_raw, r_max=r_max, s_max=s_max_mm),
        "momentum.png":     plot_momentum(optics_abs),
        "a_pz_corr.png":    plot_amplitude_pz_corr(optics_abs),
    }

    if save:
        os.makedirs(save, exist_ok=True)

    for fname, fig in figs.items():
        if save:
            fig.savefig(os.path.join(save, fname))
        if show:
            plt.show()
        plt.close(fig)

    if save:
        print(f"Saved {len(figs)} figures to {save}/")


# ---------------------------------------------------------------------------
# Phase-space movie
# ---------------------------------------------------------------------------

def make_phase_space_movie(df, out_path="animation_ps.gif", fps=5, figsize=(20, 10)):
    """
    Render a phase space animation over stations, saved as a GIF.

    One frame per station (sorted by mean z), showing:
        top-left:     x vs px
        top-right:    t vs p
        bottom-left:  x vs p
        bottom-right: y vs p

    Requires Pillow (pip install Pillow).
    """
    import matplotlib.animation as animation

    writer = animation.PillowWriter(fps=fps)
    mean_z = df.groupby("samplerName")["z"].mean()
    stations = [(name, df[df["samplerName"] == name])
                for name in mean_z.sort_values().index]

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    with writer.saving(fig, out_path, dpi=100):
        for name, group in stations:
            x, y, _ = _positions(group)
            px, _, _, p = _momenta(group)
            t = group["T"].to_numpy()
            t = t - t[0]
            mean_z_val = group["z"].mean()

            for ax in axes.flat:
                ax.cla()

            axes[0, 0].scatter(x,  px, s=2)
            axes[0, 0].set_xlabel("x [mm]")
            axes[0, 0].set_ylabel("px [MeV/c]")

            axes[0, 1].scatter(t,  p,  s=2)
            axes[0, 1].set_xlabel("t - t0 [ns]")
            axes[0, 1].set_ylabel("p [MeV/c]")

            axes[1, 0].scatter(x,  p,  s=2)
            axes[1, 0].set_xlabel("x [mm]")
            axes[1, 0].set_ylabel("p [MeV/c]")

            axes[1, 1].scatter(y,  p,  s=2)
            axes[1, 1].set_xlabel("y [mm]")
            axes[1, 1].set_ylabel("p [MeV/c]")

            fig.suptitle(f"{name}  z = {mean_z_val:.3f} m  N = {len(group)}")
            fig.tight_layout()
            writer.grab_frame()

    plt.close(fig)
    print(f"Movie saved to {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(description="Compute and plot optics from a CSV file.")
    parser.add_argument("csv",            help="Path to input CSV file")
    parser.add_argument("--out",          default="plots",   help="Output directory (default: plots/)")
    parser.add_argument("--r-max",        type=float, default=None, help="Aperture cut [same units as x, y]")
    parser.add_argument("--s-max",        type=float, default=None, help="Max arc length S")
    parser.add_argument("--transmitted",  action="store_true",      help="Keep only particles reaching last station")
    parser.add_argument("--movie",        action="store_true",      help="Also render phase-space GIF")
    args = parser.parse_args()

    matplotlib.use("Agg")   # non-interactive backend for saving

    df_raw = pd.read_csv(args.csv)

    df = apply_cuts(df_raw, r_max=args.r_max, s_max=args.s_max,
                    transmitted=args.transmitted)
    optics = compute_optics(df)

    print(optics[["samplerName", "mean_z", "emittance_trans",
                  "emittance_long", "transmission"]].to_string(index=False))

    plot_all(optics, df_raw, args.out, r_max=args.r_max, s_max=args.s_max)

    if args.movie:
        make_phase_space_movie(df, out_path=os.path.join(args.out, "animation_ps.gif"))
