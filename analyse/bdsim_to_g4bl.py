"""
Convert a BDSIM-style pandas DataFrame to G4BL/ICOOL track format (output.txt).

Expected DataFrame columns:
    samplerName  - string station label
    trackID      - combined with SimName to form a unique eventNumber (col 1)
    partID       - PDG particle code (-13 = mu+, 13 = mu-, etc.)
    x, y, z      - position in metres
    xp, yp, zp   - direction cosines (px/p, py/p, pz/p)
    p            - total momentum in GeV/c
    T            - time in ns
    weight       - statistical weight

Dummy columns written as 0: bx, by, bz, ex, ey, ez, SARC, sx, sy, sz, status

Station integers are assigned by sorting samplerName by mean z position.
"""

import pandas as pd


# PDG -> ICOOL PID mapping (only codes xboa knows about)
_PDG_TO_ICOOL = {
    -11:  1,   # e+
     11: -1,   # e-
    -13:  2,   # mu+
     13: -2,   # mu-
    211:  3,   # pi+
   -211: -3,   # pi-
    321:  4,   # K+
   -321: -4,   # K-
   2212:  5,   # proton
  -2212: -5,   # antiproton
}

_HEADER = """\
#NTuple/Z0
#Units are ns, meters, GeV/c, Tesla, and V/m
#IEVT IPNUM IPTYP IPFLG JSRG T X Y Z Px Py Pz Bx By Bz Weight Ex Ey Ez SARC POLx POLy POLz
"""


def bdsim_to_g4bl(df: pd.DataFrame, out_path: str, station_map: dict = None) -> dict:
    """
    Write df to a G4BL/ICOOL track text file at out_path.

    Parameters
    ----------
    df : DataFrame with the columns described at the top of this file.
    out_path : path to write (e.g. "tmp/optimisation/cooling_a-pz=0.9/output.txt")
    station_map : optional dict {samplerName: int} overriding the default
                  alphabetical sort. If None, names are sorted and assigned 1, 2, 3...

    Returns
    -------
    dict mapping samplerName -> station integer (so you can inspect the assignment).
    """
    if station_map is None:
        mean_z = df.groupby("samplerName")["z"].mean()
        names = mean_z.sort_values().index.tolist()
        station_map = {name: i + 1 for i, name in enumerate(names)}

    # unique integer per (trackID, SimName) pair — used as IEVT for duplicate removal
    event_ids, _ = pd.factorize(list(zip(df["trackID"], df["SimName"])))

    df = df.copy()
    df["station_int"] = df["samplerName"].map(station_map)
    df["event_id"] = event_ids
    df = df.sort_values(["station_int", "event_id"]).reset_index(drop=True)

    import os
    os.makedirs(os.path.dirname(out_path) if os.path.dirname(out_path) else ".", exist_ok=True)

    with open(out_path, "w") as fout:
        fout.write(_HEADER)
        for row in df.itertuples(index=False):
            icool_pid = _PDG_TO_ICOOL.get(int(row.partID), 0)
            station   = row.station_int
            t_s       = row.T / 1e9          # ns -> seconds
            px        = row.xp * row.p       # direction cosine * |p|  -> GeV/c
            py        = row.yp * row.p
            pz        = row.zp * row.p

            fout.write(
                f"{row.event_id} 0 {icool_pid} 0 {station} "
                f"{t_s:.10g} {row.x:.10g} {row.y:.10g} {row.z:.10g} "
                f"{px:.10g} {py:.10g} {pz:.10g} "
                f"0 0 0 "
                f"{row.weight:.6g} "
                f"0 0 0 0 "
                f"0 0 0\n"
            )

    return station_map


if __name__ == "__main__":
    # quick smoke test with a tiny synthetic DataFrame
    import numpy as np

    rng = np.random.default_rng(0)
    n = 100
    p = 0.2 + rng.normal(0, 0.01, n)
    xp = rng.normal(0, 0.05, n)
    yp = rng.normal(0, 0.05, n)
    zp = np.sqrt(np.clip(1 - xp**2 - yp**2, 0, None))

    df_test = pd.DataFrame({
        "samplerName": rng.choice(["Sampler0", "Sampler1", "Sampler2"], n),
        "turnNumber":  np.zeros(n, int),
        "trackID":     np.arange(n),
        "parentID":    np.zeros(n, int),
        "partID":      np.full(n, -13),   # mu+
        "x":           rng.normal(0, 0.01, n),
        "y":           rng.normal(0, 0.01, n),
        "xp":          xp,
        "yp":          yp,
        "z":           np.zeros(n),
        "zp":          zp,
        "energy":      np.sqrt((105.658e-3)**2 + p**2),
        "p":           p,
        "T":           rng.uniform(0, 1, n),
        "weight":      np.ones(n),
        "S":           np.zeros(n),
        "SimName":     ["test"] * n,
    })

    mapping = bdsim_to_g4bl(df_test, "/tmp/test_bdsim_to_g4bl.txt")
    print("Station map:", mapping)
    print("First 5 lines of output:")
    with open("/tmp/test_bdsim_to_g4bl.txt") as f:
        for i, line in enumerate(f):
            if i >= 8:
                break
            print(line, end="")
