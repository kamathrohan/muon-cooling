import numpy as np
from jinja2 import Template


def fmt_array(values, fmt='.6g'):
    if len(values) == 0:
        return '{0.0}'
    formatted = [format(v, fmt) for v in values]
    if len(set(formatted)) == 1:
        return '{' + formatted[0] + '}'
    return '{' + ', '.join(formatted) + '}'


def fmt_mm_array(values_m):
    if len(values_m) == 0:
        return '{0.0*mm}'
    formatted = [f'{v*1e3:.4f}*mm' for v in values_m]
    if len(set(formatted)) == 1:
        return '{' + formatted[0] + '}'
    return '{' + ', '.join(formatted) + '}'


def fmt_str_array(values, injection=None):
    if len(values) == 0:
        default = injection if injection is not None else "G4_U"
        return '{' + f'"{default}"' + '}'
    if len(set(values)) == 1:
        return '{"' + values[0] + '"}'
    return '{' + ', '.join(f'"{v}"' for v in values) + '}'


def generate_sampler_lines(n_samplers: int,
                           total_length_m: float = None,
                           aper_m: float = 5.0,
                           reference_element: str = "mc1",
                           reference_element_number: int = 0,
                           start_m: float = None,
                           end_m: float = None) -> str:
    """Generate BDSIM samplerplacement lines evenly spaced across the channel.

    Parameters
    ----------
    n_samplers              : int   - Number of samplers to place
    total_length_m          : float - Full channel length [m]; samplers span ±half (ignored if start_m/end_m given)
    aper_m                  : float - Rectangular aperture half-size [m]
    reference_element       : str   - BDSIM element name to reference
    reference_element_number: int   - referenceElementNumber value
    start_m                 : float - Start position [m]; overrides total_length_m
    end_m                   : float - End position [m]; overrides total_length_m

    Returns
    -------
    str - Multi-line BDSIM sampler placement block
    """
    if start_m is not None and end_m is not None:
        positions = np.linspace(start_m * 1e3, end_m * 1e3, n_samplers)
    else:
        half_mm = total_length_m / 2 * 1e3
        positions = np.linspace(-half_mm, half_mm, n_samplers)
    lines = []
    for i, s in enumerate(positions, 1):
        lines.append(
            f's{i}: samplerplacement, referenceElement="{reference_element}", '
            f'referenceElementNumber={reference_element_number}, '
            f's={s:.1f}*mm, shape="rectangular",'
            f'aper1={aper_m}*m, aper2={aper_m}*m;'
        )
    return '\n'.join(lines)


def generate_sampler_lines_from_positions(positions_m,
                                          aper_m: float = 5.0,
                                          reference_element: str = "mc1",
                                          reference_element_number: int = 0) -> str:
    """Generate BDSIM samplerplacement lines at arbitrary positions.

    Parameters
    ----------
    positions_m              : array-like - Sampler positions [m]
    aper_m                   : float      - Rectangular aperture half-size [m]
    reference_element        : str        - BDSIM element name to reference
    reference_element_number : int        - referenceElementNumber value

    Returns
    -------
    str - Multi-line BDSIM sampler placement block
    """
    lines = []
    for i, s_m in enumerate(positions_m, 1):
        lines.append(
            f's{i}: samplerplacement, referenceElement="{reference_element}", '
            f'referenceElementNumber={reference_element_number}, '
            f's={s_m * 1e3:.1f}*mm, shape="rectangular",'
            f'aper1={aper_m}*m, aper2={aper_m}*m;'
        )
    return '\n'.join(lines)


def generate_beam_block(mode: str,
                        particle: str = "mu+",
                        momentum_MeV: float = 200.0,
                        x0_mm: float = 0.0,
                        y0_mm: float = 0.0,
                        z0_m: float = 0.0,
                        distr_type: str = "userfile",
                        distr_file: str = "",
                        distr_file_format: str = "t[ns]:x[m]:y[m]:z[m]:xp[rad]:yp[rad]:zp[rad]:-:E[GeV]",
                        nlines_ignore: int = 1) -> str:
    """Generate a BDSIM beam block.

    Parameters
    ----------
    mode             : str   - One of 'reference', 'offset', or 'beam'
                               'reference' : particle + momentum only
                               'offset'    : particle + momentum + X0/Y0/Z0
                               'beam'      : particle + momentum + distribution file
    particle         : str   - Particle type (default "mu+")
    momentum_MeV     : float - Reference momentum [MeV/c]
    x0_mm, y0_mm     : float - Transverse offsets [mm] (offset mode)
    z0_m             : float - Longitudinal offset [m] (offset mode)
    distr_type       : str   - BDSIM distribution type (beam mode)
    distr_file       : str   - Path to distribution file (beam mode)
    distr_file_format: str   - Column format string (beam mode)
    nlines_ignore    : int   - Header lines to skip in distribution file (beam mode)

    Returns
    -------
    str - BDSIM beam command block
    """
    if mode == "reference":
        return (
            f'beam, particle="{particle}",\n'
            f'      momentum={momentum_MeV}*MeV;'
        )
    elif mode == "offset":
        return (
            f'beam, particle="{particle}",\n'
            f'      momentum={momentum_MeV}*MeV,\n'
            f'      X0={x0_mm}*mm,\n'
            f'      Y0={y0_mm}*mm,\n'
            f'      Z0={z0_m}*m;'
        )
    elif mode == "beam":
        return (
            f'beam, particle="{particle}",\n'
            f'      momentum={momentum_MeV}*MeV,\n'
            f'      distrType="{distr_type}",\n'
            f'      distrFile="{distr_file}",\n'
            f'      distrFileFormat = "{distr_file_format}",\n'
            f'      nlinesIgnore={nlines_ignore};'
        )
    else:
        raise ValueError(f"Unknown beam mode '{mode}'. Choose 'reference', 'offset', or 'beam'.")


def render_template(tpl_path: str, context: dict, out_path: str) -> None:
    with open(tpl_path) as f:
        rendered = Template(f.read()).render(**context)
    with open(out_path, 'w') as f:
        f.write(rendered)
