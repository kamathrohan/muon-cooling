Li6: matdef, Z=3, A=6.0151214, density=0.46, state="solid";
MICE_LiH: matdef, density=0.69, state="solid", components=["Li6","G4_Li","G4_H"], componentsFractions={0.814,0.043,0.143};
MICE_LIQUID_HYDROGEN: matdef, Z=1, A=1.008, density=0.07053, state="liquid";

cooldef1: coolingchannel,
	nCoils                    = {{ n_coils }},
	coilInnerRadius           = {{ coil_inner_radius }},
	coilRadialThickness       = {{ coil_radial_thickness }},
	coilLengthZ               = {{ coil_length_z }},
	coilCurrent               = {{ coil_current }},
	coilOffsetX               = {{ coil_offset_x }},
	coilOffsetY               = {{ coil_offset_y }},
	coilOffsetZ               = {{ coil_offset_z }},
	coilTiltX                 = {{ coil_tilt_x }},
	coilTiltY                 = {{ coil_tilt_y }},
	coilTiltZ                 = {{ coil_tilt_z }},
	coilMaterial              = {{ coil_material }},
	onAxisTolerance           = {{ on_axis_tolerance }},
	nDipoles                  = {{ n_dipoles }},
	dipoleAperture            = {{ dipole_aperture }},
	dipoleLengthZ             = {{ dipole_length_z }},
	dipoleFieldStrength        = {{ dipole_field_strength }},
	dipoleEngeCoefficient     = {{ dipole_enge_coefficient }},
	dipoleOffsetZ             = {{ dipole_offset_z }},
	nAbsorbers                = {{ n_absorbers }},
	absorberType              = {{ absorber_type }},
	absorberMaterial          = {{ absorber_material }},
	absorberOffsetZ           = {{ absorber_offset_z }},
	absorberCylinderLength    = {{ absorber_cylinder_length }},
	absorberCylinderRadius    = {{ absorber_cylinder_radius }},
	absorberWedgeOpeningAngle = {{ absorber_wedge_opening_angle }},
	absorberWedgeHeight       = {{ absorber_wedge_height }},
	absorberWedgeRotationAngle = {{ absorber_wedge_rotation_angle }},
	absorberWedgeOffsetX      = {{ absorber_wedge_offset_x }},
	absorberWedgeOffsetY      = {{ absorber_wedge_offset_y }},
	absorberWedgeApexToBase   = {{ absorber_wedge_apex_to_base }},
	nRFCavities               = {{ n_rf_cavities }},
	rfOffsetZ                 = {{ rf_offset_z }},
	rfTimeOffset              = {{ rf_time_offset }},
	rfLength                  = {{ rf_length }},
	rfVoltage                 = {{ rf_voltage }},
	rfPhase                   = {{ rf_phase }},
	rfFrequency               = {{ rf_frequency }},
	rfWindowThickness         = {{ rf_window_thickness }},
	rfWindowMaterial          = {{ rf_window_material }},
	rfWindowRadius            = {{ rf_window_radius }},
	rfCavityMaterial          = {{ rf_cavity_material }},
	rfCavityVacuumMaterial    = {{ rf_cavity_vacuum_material }},
	rfCavityRadius            = {{ rf_cavity_radius }},
	rfCavityThickness         = {{ rf_cavity_thickness }},
	magneticFieldModel        = "{{ magnetic_field_model }}",
	magneticFieldMethod       = "{{ magnetic_field_method }}",
	dipoleFieldModel          = "{{ dipole_field_model }}",
	interpolator              = "{{ interpolator }}",
	electricFieldModel        = "{{ electric_field_model }}";


mc1: muoncooler, l={{ total_length }}*m, horizontalWidth={{ total_width }}*m, coolingDefinition="cooldef1";
d1: drift, l=1*mm;
c1: dump, l=1*mm;
lat: line=(mc1,d1,c1);
use, period=lat;


{{ sampler_lines }}


option, checkOverlaps=1,
        stopSecondaries=1,
        storeTrajectories=1,
        collimatorsAreInfiniteAbsorbers=1,
        maximumStepLength=20*mm,
        integratorSet="geant4",
        physicsList="qgsp_bic em muon";

{{ beam_block }}
      



