Li6: matdef, Z=3, A=6.0151214, density=0.46, state="solid";
MICE_LiH: matdef, density=0.69, state="solid", components=["Li6","G4_Li","G4_H"], componentsFractions={0.814,0.043,0.143};
MICE_LIQUID_HYDROGEN: matdef, Z=1, A=1.008, density=0.07053, state="liquid";

cooldef1:  coolingchannel,
               nCoils={{ n_coils }},
    	        	    coilInnerRadius = {{ coil_inner_radius }},
          	       	    coilRadialThickness = {{ coil_radial_thickness }},
		                coilLengthZ = {{ coil_length_z }},
                        coilCurrent = {{ coil_current }},
                        coilOffsetX= {{ coil_offset_x }},
                        coilOffsetY = {{ coil_offset_y }},
                        coilOffsetZ = {{ coil_offset_z }},
                        coilTiltX = {{ coil_tilt_x }},
                        coilTiltY = {{ coil_tilt_y }},
                        coilTiltZ = {{ coil_tilt_z }},
                        coilMaterial = {"G4_Cu"},
                        onAxisTolerance = 2e-2,
                        nDipoles={{ n_dipoles }},
                        dipoleAperture= {{ dipole_aperture }},
                        dipoleLengthZ = {{ dipole_length_z }},
                        dipoleFieldStrength = {{ dipole_field_strength }},
                        dipoleEngeCoefficient = {{ dipole_enge_coefficient }},
                        dipoleOffsetZ = {{ dipole_offset_z }},
			   		nAbsorbers = {{ n_absorbers }},
			   		absorberType = {{ absorber_type }},
                        absorberMaterial = {{ absorber_material }},
                        absorberOffsetZ = {{ absorber_offset_z }},
                        absorberCylinderLength = {{ absorber_cylinder_length }},
			   		absorberCylinderRadius = {{ absorber_cylinder_radius }},
			   		absorberWedgeOpeningAngle = {{ absorber_wedge_opening_angle }},
			   		absorberWedgeHeight = {{ absorber_wedge_height }},
			   		absorberWedgeRotationAngle = {{ absorber_wedge_rotation_angle }},
                        absorberWedgeOffsetX = {{ absorber_wedge_offset_x }},
                    absorberWedgeOffsetY = {{ absorber_wedge_offset_y }},
			   		absorberWedgeApexToBase = {{ absorber_wedge_apex_to_base }},
			  			nRFCavities = {{ n_rf_cavities }},
                        rfOffsetZ = {{ rf_offset_z }},
                        rfTimeOffset = {{ rf_time_offset }},
                        rfLength={{ rf_length }},
                        rfVoltage={{ rf_voltage }},
                        rfPhase={{ rf_phase }},
                        rfFrequency={{ rf_frequency }},
			   		rfWindowThickness = {{ rf_window_thickness }},
			   		rfWindowMaterial = {{ rf_window_material }},
			   		rfWindowRadius = {{ rf_window_radius }},
			   		rfCavityMaterial = {{ rf_cavity_material }},
			   		rfCavityVacuumMaterial = {{ rf_cavity_vacuum_material }},
			   		rfCavityRadius = {{ rf_cavity_radius }},
			   		rfCavityThickness = {{ rf_cavity_thickness }},

    					magneticFieldModel="{{ magnetic_field_model }}",
                            magneticFieldMethod="{{ magnetic_field_method }}",
                            dipoleFieldModel="{{ dipole_field_model }}",
                            interpolator="{{ interpolator }}",
    					electricFieldModel="{{ electric_field_model }}";


mc1: muoncooler, l={{ total_length }}*m, horizontalWidth={{ total_width }}*m, coolingDefinition="cooldef1";
d1: drift, l=1*mm;
c1: dump, l=1*mm;
lat: line=(mc1,d1,c1);
use, period=lat;


s1: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-51600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s2: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-50800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s3: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-50000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s4: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-49200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s5: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-48400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s6: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-47600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s7: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-46800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s8: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-46000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s9: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-45200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s10: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-44400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s11: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-43600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s12: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-42800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s13: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-42000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s14: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-41200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s15: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-40400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s16: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-39600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s17: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-38800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s18: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-38000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s19: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-37200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s20: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-36400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s21: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-35600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s22: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-34800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s23: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-34000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s24: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-33200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s25: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-32400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s26: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-31600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s27: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-30800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s28: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-30000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s29: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-29200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s30: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-28400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s31: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-27600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s32: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-26800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s33: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-26000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s34: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-25200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s35: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-24400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s36: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-23600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s37: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-22800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s38: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-22000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s39: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-21200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s40: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-20400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s41: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-19600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s42: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-18800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s43: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-18000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s44: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-17200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s45: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-16400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s46: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-15600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s47: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-14800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s48: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-14000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s49: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-13200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s50: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-12400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s51: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-11600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s52: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-10800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s53: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-10000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s54: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-9200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s55: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-8400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s56: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-7600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s57: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-6800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s58: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-6000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s59: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-5200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s60: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-4400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s61: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-3600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s62: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-2800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s63: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-2000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s64: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-1200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s65: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=-400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s66: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s67: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=1200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s68: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=2000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s69: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=2800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s70: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=3600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s71: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=4400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s72: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=5200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s73: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=6000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s74: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=6800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s75: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=7600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s76: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=8400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s77: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=9200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s78: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=10000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s79: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=10800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s80: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=11600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s81: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=12400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s82: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=13200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s83: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=14000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s84: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=14800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s85: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=15600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s86: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=16400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s87: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=17200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s88: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=18000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s89: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=18800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s90: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=19600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s91: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=20400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s92: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=21200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s93: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=22000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s94: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=22800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s95: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=23600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s96: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=24400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s97: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=25200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s98: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=26000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s99: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=26800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s100: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=27600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s101: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=28400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s102: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=29200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s103: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=30000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s104: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=30800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s105: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=31600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s106: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=32400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s107: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=33200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s108: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=34000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s109: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=34800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s110: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=35600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s111: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=36400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s112: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=37200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s113: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=38000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s114: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=38800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s115: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=39600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s116: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=40400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s117: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=41200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s118: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=42000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s119: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=42800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s120: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=43600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s121: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=44400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s122: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=45200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s123: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=46000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s124: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=46800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s125: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=47600.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s126: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=48400.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s127: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=49200.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s128: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=50000.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;
s129: samplerplacement, referenceElement="mc1", referenceElementNumber=0, s=50800.0*mm, shape="rectangular",aper1=5*m, aper2=5*m;


option, checkOverlaps=1,
        stopSecondaries=1,
        storeTrajectories=1,
        collimatorsAreInfiniteAbsorbers=1,
        maximumStepLength=20*mm,
        integratorSet="geant4",
        physicsList="qgsp_bic em muon";

!beam, particle="mu+", momentum=200.0*MeV;

beam, particle="mu+",
      momentum=200*MeV,
      !X0=8.62*mm, !8.66727827153607 8.62
      !Y0=7.02*mm, !7.090012632655423 7.02
      !Z0=10.4*m; !9.6
      distrType="userfile",
      distrFile="highemittanceBeamGen_1.5mm_z10.4_offsetXY.dat",
      !distrFile="highemittanceBeamGen_2.5mm_z10.4_offsetXY.dat",
      !distrFile="lowemittanceBeamGen_z10.4_offsetXY.dat",
      distrFileFormat = "t[ns]:x[m]:y[m]:z[m]:xp[rad]:yp[rad]:zp[rad]:-:E[GeV]",
      nlinesIgnore=1;
      



