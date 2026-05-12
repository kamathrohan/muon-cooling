from ..physics.elements import Dipole, Absorber, RFCavity
from ..physics.beamline import MuonCoolingChannel
from ..physics.translate import MUON_MASS_MEV as MUON_MASS_MEV_C2
from .template import (
    fmt_array, fmt_mm_array, fmt_str_array,
    generate_sampler_lines, generate_sampler_lines_from_positions,
    generate_beam_block, render_template,
)


def render_gmad(beamline: MuonCoolingChannel,
                tpl_path: str,
                out_gmad: str,
                sampler_mode: str = "linspace",
                sampler_kwargs: dict = None,
                beam_mode: str = "beam",
                beam_kwargs: dict = None) -> None:
    """Render a full channel gmad from a template and a MuonCoolingChannel object.

    Parameters
    ----------
    beamline      : MuonCoolingChannel
    tpl_path      : str  - path to config/channel.tpl
    out_gmad      : str  - output gmad path
    sampler_mode  : str  - "linspace" or "positions"
        "linspace"  : evenly spaced samplers; sampler_kwargs keys:
                        n        (int,   default 129)
                        aper_m   (float, default 5.0)
                        start_m  (float, optional — defaults to -total_length/2)
                        end_m    (float, optional — defaults to +total_length/2)
        "positions" : arbitrary positions; sampler_kwargs keys:
                        positions_m (array-like, required)
                        aper_m      (float, default 5.0)
    sampler_kwargs: dict - keyword args forwarded to the sampler generator
    beam_mode     : str  - forwarded to generate_beam_block
    beam_kwargs   : dict - keyword args forwarded to generate_beam_block
    """
    half_length = beamline.total_length / 2
    global_offset = half_length
    for elem in beamline.elements:
        if not (-half_length <= elem.z_center <= half_length):
            z_global = elem.z_center + global_offset
            raise ValueError(
                f"Element '{elem.name}' at local z={elem.z_center}m (global z={z_global}m) "
                f"exceeds channel bounds ±{half_length}m (total_length={beamline.total_length}m)"
            )

    rf_cavs = [e for e in beamline.elements if isinstance(e, RFCavity)]
    unphased = [c.name for c in rf_cavs if c.time_offset is None]
    if unphased:
        raise ValueError(
            f"RF cavities are not phased (time_offset is None): {unphased}. "
            "Call compute_rf_time_offsets() or set_rf_time_offsets() before rendering."
        )

    pancake_df = beamline.build_pancake_dataframe()
    df = pancake_df.sort_values(['coil_index', 'pancake_index']).reset_index(drop=True)

    dipoles   = [e for e in beamline.elements if isinstance(e, Dipole)]
    absorbers = [e for e in beamline.elements if isinstance(e, Absorber)]

    n_coils = len(df)
    context = {
        # channel properties
        'total_length':      beamline.total_length,
        'total_width':       beamline.total_width,
        'on_axis_tolerance': beamline.on_axis_tolerance,
        'coil_material':     fmt_str_array(list(df['material']), injection="G4_Cu"),

        # channel field models
        'magnetic_field_model':  beamline.magnetic_field_model,
        'magnetic_field_method': beamline.magnetic_field_method,
        'dipole_field_model':    beamline.dipole_field_model,
        'interpolator':          beamline.interpolator,
        'electric_field_model':  beamline.electric_field_model,
        'grid_points_per_mm':    beamline.grid_points_per_mm,
        'z_period_start':        beamline.z_period_start,
        'z_period_end':          beamline.z_period_end,
        'period_length':         beamline.period_length,

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
        'absorber_type':              fmt_str_array([a.absorber_type       for a in absorbers], injection="wedge"),
        'absorber_material':          fmt_str_array([a.material            for a in absorbers], injection="G4_Cu"),
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
        'rf_window_material':     fmt_str_array([r.window_material for r in rf_cavs], injection="G4_Cu"),
        'rf_window_radius':       fmt_mm_array([r.window_radius   for r in rf_cavs]),
        'rf_cavity_material':     fmt_str_array([r.cavity_material for r in rf_cavs], injection="G4_Cu"),
        'rf_cavity_vacuum_material': fmt_str_array([r.cavity_vacuum_material for r in rf_cavs], injection="G4_Cu"),
        'rf_cavity_radius':       fmt_mm_array([r.cavity_radius   for r in rf_cavs]),
        'rf_cavity_thickness':    fmt_mm_array([r.cavity_thickness for r in rf_cavs]),
    }

    skw = sampler_kwargs or {}
    _ref_elem   = skw.get('reference_element', 'mc1')
    _ref_num    = skw.get('reference_element_number', 0)
    _aper_m     = skw.get('aper_m', 5.0)
    if sampler_mode == "linspace":
        n_samplers = skw.get('n', 129)
        context['sampler_lines'] = generate_sampler_lines(
            n_samplers, beamline.total_length,
            aper_m=_aper_m,
            reference_element=_ref_elem,
            reference_element_number=_ref_num,
            start_m=skw.get('start_m'),
            end_m=skw.get('end_m'),
        )
    elif sampler_mode == "positions":
        positions_m = skw['positions_m']
        n_samplers = len(positions_m)
        context['sampler_lines'] = generate_sampler_lines_from_positions(
            positions_m,
            aper_m=_aper_m,
            reference_element=_ref_elem,
            reference_element_number=_ref_num,
        )
    else:
        raise ValueError(f"Unknown sampler_mode '{sampler_mode}'. Choose 'linspace' or 'positions'.")

    context['beam_block'] = generate_beam_block(
        beam_mode, momentum_MeV=beamline.reference_momentum, **(beam_kwargs or {})
    )

    render_template(tpl_path, context, out_gmad)

    print(f"Written {out_gmad}: {n_coils} coils, {len(dipoles)} dipoles, "
          f"{len(absorbers)} absorbers, {len(rf_cavs)} RF cavities, "
          f"{n_samplers} samplers")
