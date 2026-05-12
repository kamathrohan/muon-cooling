import json
import numpy as np
from src import (
    MuonCoolingChannel,
    build_coil_beamline,
    build_dipole_beamline,
    build_absorber_beamline,
    build_rf_beamline,
    render_gmad,
)


channel = MuonCoolingChannel(
    n_cells            = 151,
    cell_length        = 1.0,
    total_length       = 170.0,
    total_width        = 0.8,
    reference_momentum = 200.0,
    on_axis_tolerance  = 2e-5,
    magnetic_field_model="solenoidblock",
    magnetic_field_method="grid",
    grid_points_per_mm = 1.0,
    z_period_start     = -70000,
    z_period_end       = 175000,
    period_length      = 2000,
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
    'wedge_height': 0.28571,
    'wedge_apex_to_base': 0.28571
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
    coil_cell_z=0.025,
    dipole_template=dipole_template,
    channel=channel,
)

channel = build_absorber_beamline(
    n_cells = 129,
    absorber_template=absorber_template,
    channel=channel, wedge_alignment_angle1=0.0, wedge_alignment_angle2=np.pi
)

channel = build_rf_beamline(
    n_cells= 129,
    n_rf_cells=3,
    rf_spacing=0.1946,
    rf_template=rf_template,
    channel=channel,
)





n_coils = sum(1 for e in channel.elements if hasattr(e, 'tilt_x'))
channel.set_tilts({"tiltX": [10.0] * n_coils, "tiltY": [20.0] * n_coils, "tiltZ": [30.0] * n_coils})

with open("config/rfPhases.json") as f:
    rf_phases = json.load(f)
channel.set_rf_time_offsets(rf_phases["closedOrbit"])
render_gmad(channel, "config/channel.tpl", "output/channeltest.gmad",
            sampler_mode="linspace",
            sampler_kwargs={"n": 200, "start_m": -65, "end_m": 65},
            beam_mode="beam", beam_kwargs={"distr_file": "beam_bdsim4.dat"})
