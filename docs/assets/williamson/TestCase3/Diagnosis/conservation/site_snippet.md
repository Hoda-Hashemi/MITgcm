# TC3 Conservation Diagnostics

**Experiment:** Steady compact-support zonal flow

## Suggested Section Copy

Mass is preserved to roundoff (max |rel| 1.172e-15). The mechanical-energy proxy is not preserved (max |rel| 0.0830516). Derived potential enstrophy is not preserved (max |rel| 0.0365314).

These diagnostics use the current MITgcm MDS outputs only. Mass and free-surface drift come from `Eta/ETAN`; energy uses `momKE` when present, otherwise centered `U/V`; PV and potential enstrophy are derived from `U/V/Eta` because no run currently writes a native vorticity or PV diagnostic.

## Availability

- Configured diagnostics: UVEL, VVEL, ETAN, ETANSQ, momKE, dynDiag, WVELMASS, PhiVEL, PsiVEL, dyn_Aux, D, dynStDiag
- Runs with conservation tables: 2

## Assets To Link

- alpha 0: Mass/free surface drift -> `assets/williamson/TestCase3/Diagnosis/conservation/alpha_0/conservation_mass.png`
- alpha 0: Energy proxy drift -> `assets/williamson/TestCase3/Diagnosis/conservation/alpha_0/conservation_energy.png`
- alpha 0: PV/enstrophy drift -> `assets/williamson/TestCase3/Diagnosis/conservation/alpha_0/conservation_pv_enstrophy.png`
- alpha 0: Conservation time series table -> `assets/williamson/TestCase3/Diagnosis/conservation/alpha_0/conservation_timeseries.csv`
- alpha 0: Conservation summary data -> `assets/williamson/TestCase3/Diagnosis/conservation/alpha_0/conservation_summary.json`
- alpha 1.0472: Mass/free surface drift -> `assets/williamson/TestCase3/Diagnosis/conservation/alpha_1.0472/conservation_mass.png`
- alpha 1.0472: Energy proxy drift -> `assets/williamson/TestCase3/Diagnosis/conservation/alpha_1.0472/conservation_energy.png`
- alpha 1.0472: PV/enstrophy drift -> `assets/williamson/TestCase3/Diagnosis/conservation/alpha_1.0472/conservation_pv_enstrophy.png`
- alpha 1.0472: Conservation time series table -> `assets/williamson/TestCase3/Diagnosis/conservation/alpha_1.0472/conservation_timeseries.csv`
- alpha 1.0472: Conservation summary data -> `assets/williamson/TestCase3/Diagnosis/conservation/alpha_1.0472/conservation_summary.json`

## Per-Run Verdicts

- alpha 0: mass preserved to roundoff; quantity small drift; energy small drift; enstrophy small drift; health finite state fields.
- alpha 1.0472: mass preserved to roundoff; quantity not preserved; energy not preserved; enstrophy not preserved; health finite state fields.
