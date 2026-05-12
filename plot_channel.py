"""
CLI tool: load sampler files from a directory, compute optics, and show/save plots.

Usage
-----
python plot_channel.py <directory> [options]

Examples
--------
python plot_channel.py /path/to/sims
python plot_channel.py /path/to/sims --s-offset 21.5 --r-max 0.0819 --s-max 100 --save plots/
python plot_channel.py /path/to/sims --show --no-save
python plot_channel.py /path/to/sims --cell-length 2.0 --movie
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from analyse.loader import find_files, load_sims
from analyse.optics import apply_cuts, compute_optics, plot_all, make_phase_space_movie


def run(directory, pattern="*.txt", s_offset=0.0, r_max=None, s_max=None,
        transmitted=False, cell_length=None, show=False, save="plots", movie=False):
    """
    Load sampler files, compute optics, and show/save plots.

    Parameters
    ----------
    directory   : str        - directory to search for files
    pattern     : str        - glob pattern (default "*.txt")
    s_offset    : float      - subtracted from S column [mm]
    r_max       : float      - transverse aperture cut [m]
    s_max       : float      - zoom all plots to z <= s_max [m]; also cuts data
    transmitted : bool       - keep only particles reaching the last station
    cell_length : float      - if given, beta/momentum/a_pz_corr use only stations
                               at multiples of cell_length [m] (e.g. 2.0)
    show        : bool       - display figures interactively
    save        : str or None - directory to save figures into; None skips saving
    movie       : bool       - also render a phase-space GIF

    Returns
    -------
    optics : DataFrame from compute_optics()
    df     : filtered DataFrame used for optics calculation
    dfAll  : raw concatenated DataFrame
    """
    files = find_files(directory, pattern=pattern)
    if not files:
        raise FileNotFoundError(f"No files matching '{pattern}' found in {directory}")
    print(f"Found {len(files)} file(s)")

    dfAll = load_sims(files)

    if s_offset:
        dfAll["S"] = dfAll["S"] - s_offset
    dfAll["z"] = dfAll["S"]

    # apply_cuts uses mm internally; s_max is in metres
    s_max_mm = s_max * 1000 if s_max is not None else None
    df = apply_cuts(dfAll, r_max=r_max, s_max=s_max_mm, transmitted=transmitted)
    optics = compute_optics(df)

    plot_all(optics, dfAll, r_max=r_max, s_max=s_max, cell_length=cell_length,
             show=show, save=save)

    if movie:
        movie_path = os.path.join(save, "phase_space.gif") if save else "phase_space.gif"
        make_phase_space_movie(df, out_path=movie_path)

    return optics, df, dfAll


def main():
    parser = argparse.ArgumentParser(
        description="Load sampler files, compute optics, and show/save plots."
    )
    parser.add_argument("directory",
                        help="Directory to search for files (recursive)")
    parser.add_argument("--pattern",      default="*.txt",
                        help="Glob pattern for input files (default: *.txt)")
    parser.add_argument("--s-offset",     type=float, default=0.0,
                        help="Offset subtracted from S column [mm] (default: 0)")
    parser.add_argument("--r-max",        type=float, default=None,
                        help="Transverse aperture cut [m] (default: no cut)")
    parser.add_argument("--s-max",        type=float, default=None,
                        help="Zoom all plots to z <= s-max [m] (default: no cut)")
    parser.add_argument("--transmitted",  action="store_true",
                        help="Keep only particles that reach the last station")
    parser.add_argument("--cell-length",  type=float, default=None, metavar="M",
                        help="Show beta/momentum/a_pz_corr at cell boundaries only [m]")
    parser.add_argument("--show",         action="store_true",
                        help="Display figures interactively")
    parser.add_argument("--save",         default="plots", metavar="DIR",
                        help="Directory to save figures into (default: plots/)")
    parser.add_argument("--no-save",      action="store_true",
                        help="Disable saving figures to disk")
    parser.add_argument("--movie",        action="store_true",
                        help="Also render a phase-space GIF")
    args = parser.parse_args()

    run(
        directory=args.directory,
        pattern=args.pattern,
        s_offset=args.s_offset,
        r_max=args.r_max,
        s_max=args.s_max,
        transmitted=args.transmitted,
        cell_length=args.cell_length,
        show=args.show,
        save=None if args.no_save else args.save,
        movie=args.movie,
    )


if __name__ == "__main__":
    main()
