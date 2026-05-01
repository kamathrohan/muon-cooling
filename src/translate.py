import numpy as np

MUON_MASS_MEV = 105.66


def g4bl_to_beamgen(input_path: str, output_path: str,
                    mass_MeV: float = MUON_MASS_MEV,
                    z_override: float = None) -> int:
    """Translate a G4Beamline BLTrackFile (beam.tmp) to the highemittanceBeamGen format.

    G4BL columns: x y z Px Py Pz t PDGid EvNum TrkId Parent weight
      - x, y, z  : [mm]
      - Px,Py,Pz : [MeV/c]
      - t        : [ns]

    Output columns: T X Y Z xp yp zp P E
      - T        : time [ns]
      - X, Y, Z  : position [m]
      - xp,yp,zp : direction cosines
      - P        : total momentum [GeV/c]
      - E        : total energy [GeV]

    Parameters
    ----------
    input_path  : path to G4BL BLTrackFile
    output_path : path to write beamgen file
    mass_MeV    : particle rest mass [MeV/c²] (default: muon)
    z_override  : if set, replace all Z values with this constant [m]

    Returns
    -------
    int - number of particles written
    """
    cols = None
    with open(input_path) as fh:
        for line in fh:
            if line.startswith('#'):
                tokens = line.lstrip('#').split()
                if tokens and tokens[0].lower() in ('x', 't'):
                    cols = [t.lower() for t in tokens]
            else:
                break

    # fall back to G4BL default column order
    if cols is None:
        cols = ['x', 'y', 'z', 'px', 'py', 'pz', 't', 'pdgid', 'evnum', 'trkid', 'parent', 'weight']

    def col(name):
        return data[:, cols.index(name)]

    data = np.loadtxt(input_path, comments='#')
    if data.ndim == 1:
        data = data[np.newaxis, :]

    x_mm, y_mm, z_mm = col('x'), col('y'), col('z')
    Px, Py, Pz       = col('px'), col('py'), col('pz')
    t                = col('t')

    X = x_mm * 1e-3
    Y = y_mm * 1e-3
    Z = np.full(len(x_mm), z_override) if z_override is not None else z_mm * 1e-3

    P_MeV = np.sqrt(Px**2 + Py**2 + Pz**2)
    xp = Px / P_MeV
    yp = Py / P_MeV
    zp = Pz / P_MeV

    P_GeV = P_MeV * 1e-3
    E_GeV = np.sqrt(P_GeV**2 + (mass_MeV * 1e-3)**2)

    out = np.column_stack([t, X, Y, Z, xp, yp, zp, P_GeV, E_GeV])

    header = "T X Y Z xp yp zp P E"
    np.savetxt(output_path, out, header=header, comments='')

    return len(out)
