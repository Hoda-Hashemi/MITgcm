MITgcm diagnostics_calc_phivel extraction summary
run_dir=/Users/hodahashemi/Documents/Github/Playground/MITgcm/verification/vortexSphere_mitgcm_referenceCase/run
iteration=22

What fields were extracted:
  dynDiag record 1: UVELMASS
  dynDiag record 2: VVELMASS
  dyn_Aux  record 2: PhiVEL
  dyn_Aux  record 3: PsiVEL

Velocity/transport meaning:
  UVELMASS and VVELMASS are the horizontal face-normal input velocities
  used by diagnostics_calc_phivel. They are in m/s, not transports.
  In diagnostics_fill_state they are formed as:
    UVELMASS = uVel * hFacW
    VVELMASS = vVel * hFacS

Transport computed internally by diagnostics_calc_phivel:
  surface DRF(1) = 1.0000000000000000e+03 m
  uTrans = DYG * DRF(1) * UVELMASS
  vTrans = DXG * DRF(1) * VVELMASS
  Units of uTrans/vTrans are m^3/s.

Files written:
  phi_csv=/Users/hodahashemi/Documents/Github/Playground/MITgcm/verification/vortexSphere_mitgcm_referenceCase/Scripts/results_referenceCase/extract_it0000000022/PhiOcean.00000000.csv
  psi_csv=/Users/hodahashemi/Documents/Github/Playground/MITgcm/verification/vortexSphere_mitgcm_referenceCase/Scripts/results_referenceCase/extract_it0000000022/PsiOcean.00000000.csv
  input_velocity_csv=/Users/hodahashemi/Documents/Github/Playground/MITgcm/verification/vortexSphere_mitgcm_referenceCase/Scripts/results_referenceCase/extract_it0000000022/VelMassOcean.00000000.csv
  input_transport_csv=/Users/hodahashemi/Documents/Github/Playground/MITgcm/verification/vortexSphere_mitgcm_referenceCase/Scripts/results_referenceCase/extract_it0000000022/TransportOcean.00000000.csv
  bundle_npz=/Users/hodahashemi/Documents/Github/Playground/MITgcm/verification/vortexSphere_mitgcm_referenceCase/Scripts/results_referenceCase/extract_it0000000022/dyn_fields_it0000000022.npz
  summary_txt=/Users/hodahashemi/Documents/Github/Playground/MITgcm/verification/vortexSphere_mitgcm_referenceCase/Scripts/results_referenceCase/extract_it0000000022/README.txt
