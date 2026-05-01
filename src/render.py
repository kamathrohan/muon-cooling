import numpy as np
from jinja2 import Template

from .elements import Dipole, Absorber, RFCavity
from .beamline import MuonCoolingChannel

MUON_MASS_MEV_C2 = 105.66


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


def render_gmad(beamline: MuonCoolingChannel,
                tpl_path: str,
                out_gmad: str,
                n_samplers: int = 129,
                sampler_aper_m: float = 5.0,
                sampler_start_m: float = None,
                sampler_end_m: float = None,
                beam_mode: str = "beam",
                beam_kwargs: dict = None) -> None:
    """Render a full channel gmad from templates/channel.tpl and a MuonCoolingChannel object.

    Parameters
    ----------
    beamline : MuonCoolingChannel - channel with Solenoid/Dipole/Absorber/RFCavity elements
    tpl_path : str          - path to templates/channel.tpl
    out_gmad : str          - output gmad path
    """
    beamline.compute_rf_time_offsets(beamline.reference_momentum, MUON_MASS_MEV_C2)

    pancake_df = beamline.build_pancake_dataframe()
    df = pancake_df.sort_values(['coil_index', 'pancake_index']).reset_index(drop=True)

    dipoles   = [e for e in beamline.elements if isinstance(e, Dipole)]
    absorbers = [e for e in beamline.elements if isinstance(e, Absorber)]
    rf_cavs   = [e for e in beamline.elements if isinstance(e, RFCavity)]

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

    def fmt_str_array(values):
        if len(values) == 0:
            return '{"G4_U"}'
        if len(set(values)) == 1:
            return '{"' + values[0] + '"}'
        return '{' + ', '.join(f'"{v}"' for v in values) + '}'

    n_coils = len(df)
    context = {
        # channel properties
        'total_length':      beamline.total_length,
        'total_width':       beamline.total_width,
        'on_axis_tolerance': beamline.on_axis_tolerance,
        'coil_material':     fmt_str_array(list(df['material'])),

        # channel field models
        'magnetic_field_model':  beamline.magnetic_field_model,
        'magnetic_field_method': beamline.magnetic_field_method,
        'dipole_field_model':    beamline.dipole_field_model,
        'interpolator':          beamline.interpolator,
        'electric_field_model':  beamline.electric_field_model,

        # coils
        'n_coils':               n_coils,
        'coil_inner_radius':     fmt_mm_array(df['r_in']),
        'coil_radial_thickness': fmt_mm_array(df['r_thick']),
        'coil_length_z':         fmt_array(df['L_pancake']),
        'coil_current':          fmt_array(df['currDensity']),
        'coil_offset_x':         fmt_array(df['coil_offset_x']),
        'coil_offset_y':         fmt_array(df['coil_offset_y']),
        'coil_offset_z':         fmt_array(df['z_center'], '.9g'),
        'coil_tilt_x':           fmt_array(df['coil_tilt_x']),
        'coil_tilt_y':           fmt_array(df['coil_tilt_y']),
        'coil_tilt_z':           fmt_array(df['coil_tilt_z']),

        # dipoles
        'n_dipoles':              len(dipoles),
        'dipole_aperture':        fmt_mm_array([d.aperture      for d in dipoles]),
        'dipole_length_z':        fmt_mm_array([d.length_z      for d in dipoles]),
        'dipole_field_strength':  fmt_array([d.field_strength   for d in dipoles]),
        'dipole_enge_coefficient':fmt_array([d.enge_coefficient for d in dipoles]),
        'dipole_offset_z':        fmt_array([d.z_center         for d in dipoles], '.9g'),

        # absorbers
        'n_absorbers':                len(absorbers),
        'absorber_type':              fmt_str_array([a.absorber_type       for a in absorbers]),
        'absorber_material':          fmt_str_array([a.material            for a in absorbers]),
        'absorber_offset_z':          fmt_array([a.z_center               for a in absorbers], '.9g'),
        'absorber_cylinder_length':   fmt_mm_array([a.cylinder_length      for a in absorbers]),
        'absorber_cylinder_radius':   fmt_mm_array([a.cylinder_radius      for a in absorbers]),
        'absorber_wedge_opening_angle': fmt_array([a.wedge_opening_angle   for a in absorbers]),
        'absorber_wedge_height':      fmt_mm_array([a.wedge_height         for a in absorbers]),
        'absorber_wedge_rotation_angle': fmt_array([a.wedge_rotation_angle for a in absorbers]),
        'absorber_wedge_offset_x':    fmt_array([a.wedge_offset_x          for a in absorbers]),
        'absorber_wedge_offset_y':    fmt_array([a.wedge_offset_y          for a in absorbers]),
        'absorber_wedge_apex_to_base':fmt_mm_array([a.wedge_apex_to_base   for a in absorbers]),

        # rf cavities
        'n_rf_cavities':          len(rf_cavs),
        'rf_offset_z':            fmt_array([r.z_center          for r in rf_cavs], '.9g'),
        'rf_time_offset':         fmt_array([r.time_offset        for r in rf_cavs]),
        'rf_length':              fmt_mm_array([r.length          for r in rf_cavs]),
        'rf_voltage':             fmt_array([r.voltage            for r in rf_cavs]),
        'rf_phase':               fmt_array([r.phase              for r in rf_cavs]),
        'rf_frequency':           fmt_array([r.frequency          for r in rf_cavs]),
        'rf_window_thickness':    fmt_mm_array([r.window_thickness for r in rf_cavs]),
        'rf_window_material':     fmt_str_array([r.window_material for r in rf_cavs]),
        'rf_window_radius':       fmt_mm_array([r.window_radius   for r in rf_cavs]),
        'rf_cavity_material':     fmt_str_array([r.cavity_material for r in rf_cavs]),
        'rf_cavity_vacuum_material': fmt_str_array([r.cavity_vacuum_material for r in rf_cavs]),
        'rf_cavity_radius':       fmt_mm_array([r.cavity_radius   for r in rf_cavs]),
        'rf_cavity_thickness':    fmt_mm_array([r.cavity_thickness for r in rf_cavs]),
    }

    context['sampler_lines'] = generate_sampler_lines(
        n_samplers, beamline.total_length, aper_m=sampler_aper_m,
        start_m=sampler_start_m, end_m=sampler_end_m,
    )
    context['beam_block'] = generate_beam_block(
        beam_mode, momentum_MeV=beamline.reference_momentum, **(beam_kwargs or {})
    )

    with open(tpl_path) as f:
        rendered = Template(f.read()).render(**context)

    with open(out_gmad, 'w') as f:
        f.write(rendered)

    print(f"Written {out_gmad}: {n_coils} coils, {len(dipoles)} dipoles, "
          f"{len(absorbers)} absorbers, {len(rf_cavs)} RF cavities, "
          f"{n_samplers} samplers")
