from src import (
    MuonCoolingChannel,
    build_coil_beamline,
    build_dipole_beamline,
    build_absorber_beamline,
    build_rf_beamline,
    render_gmad,
)


channel = MuonCoolingChannel(
    n_cells            = 10,
    cell_length        = 0.8,
    total_length       = 124.0,
    total_width        = 0.8,
    reference_momentum = 200.0,
    on_axis_tolerance  = 2e-2,
)

coil_1_template = {
    'name': 'Coil1',
    'r_in': 0.285,
    'r_thick': 0.070,
    'N_pancakes': 17,
    'L_pancake': 0.012,
    'L_spacing': 0.004,
    'currDensity': 328.43,
    'N_sheets': 5,
}

coil_2_template = {
    'name': 'Coil2',
    'r_in': 0.185,
    'r_thick': 0.060,
    'N_pancakes': 5,
    'L_pancake': 0.012,
    'L_spacing': 0.003,
    'currDensity': 300.0,
    'N_sheets': 5,
}

dipole_template = {
    'name': 'Dipole',
    'field_strength': 0.2,
    'aperture': 0.2,
    'length_z': 0.1,
    'enge_coefficient': 5.5,
}

absorber_template = {
    'name': 'Absorber',
    'absorber_type': 'wedge',
    'material': 'G4_LITHIUM_HYDRIDE',
    'wedge_opening_angle': 0.17453292519943295,
    'wedge_height': 0.24,
    'wedge_apex_to_base': 0.1134256,
}

rf_template = {
    'name': 'RF',
    'length': 0.18856,
    'voltage': 30.0,
    'phase': -1.5707963267948966 + 0.3490658503988659,
    'frequency': 704e6,
    'cavity_radius': 0.163,
    'cavity_thickness': 0.003,
    'cavity_material': 'G4_Cu',
}


channel = build_coil_beamline(
    z_coils=[0.081, 0.211],
    coil_templates=[coil_2_template, coil_1_template],
    polarities=[1, 1],
    fixed=True,
    channel=channel,
)

channel = build_dipole_beamline(
    coil_cell_z=0.7,
    dipole_template=dipole_template,
    channel=channel,
)

channel = build_absorber_beamline(
    absorber_template=absorber_template,
    channel=channel,
)

channel = build_rf_beamline(
    n_rf_cells=3,
    rf_spacing=0.1946,
    rf_template=rf_template,
    channel=channel,
)

channel.summary()

render_gmad(channel, "templates/channel.tpl", "channel.gmad", n_samplers=10)
