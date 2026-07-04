# TC7 Conservation and Output-Health Checks

**Experiment:** Analyzed 500 mb initial state

## Suggested Section Copy

Mass is preserved to roundoff (max |rel| 1.782e-16). The mechanical-energy proxy is not preserved (max |rel| 0.0114949). Derived potential enstrophy is not preserved (max |rel| 0.293603).

These diagnostics use the current MITgcm MDS outputs only. Mass and free-surface drift come from `Eta/ETAN`; energy uses `momKE` when present, otherwise centered `U/V`; PV and potential enstrophy are derived from `U/V/Eta` because no run currently writes a native vorticity or PV diagnostic.

## Note

TC7 has staged analyzed input and needs completed run output before conservation diagnostics can be evaluated.

## Availability

- Configured diagnostics: UVEL, VVEL, ETAN, ETANSQ, momKE, dynDiag, WVELMASS, PhiVEL, PsiVEL, dyn_Aux, D, dynStDiag
- Runs with conservation tables: 3

## Assets To Link

- alpha c1_19781221_0000: Mass/free surface drift -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c1_19781221_0000/conservation_mass.png`
- alpha c1_19781221_0000: Energy proxy drift -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c1_19781221_0000/conservation_energy.png`
- alpha c1_19781221_0000: PV/enstrophy drift -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c1_19781221_0000/conservation_pv_enstrophy.png`
- alpha c1_19781221_0000: Conservation time series table -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c1_19781221_0000/conservation_timeseries.csv`
- alpha c1_19781221_0000: Conservation summary data -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c1_19781221_0000/conservation_summary.json`
- alpha c2_19790116_0000: Mass/free surface drift -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c2_19790116_0000/conservation_mass.png`
- alpha c2_19790116_0000: Energy proxy drift -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c2_19790116_0000/conservation_energy.png`
- alpha c2_19790116_0000: PV/enstrophy drift -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c2_19790116_0000/conservation_pv_enstrophy.png`
- alpha c2_19790116_0000: Conservation time series table -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c2_19790116_0000/conservation_timeseries.csv`
- alpha c2_19790116_0000: Conservation summary data -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c2_19790116_0000/conservation_summary.json`
- alpha c3_19790109_0000: Mass/free surface drift -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c3_19790109_0000/conservation_mass.png`
- alpha c3_19790109_0000: Energy proxy drift -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c3_19790109_0000/conservation_energy.png`
- alpha c3_19790109_0000: PV/enstrophy drift -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c3_19790109_0000/conservation_pv_enstrophy.png`
- alpha c3_19790109_0000: Conservation time series table -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c3_19790109_0000/conservation_timeseries.csv`
- alpha c3_19790109_0000: Conservation summary data -> `assets/williamson/TestCase7/Diagnosis/conservation/alpha_c3_19790109_0000/conservation_summary.json`

## Per-Run Verdicts

- alpha c1_19781221_0000: mass preserved to roundoff; quantity not preserved; energy not preserved; enstrophy not preserved; health finite state fields.
- alpha c2_19790116_0000: mass preserved to roundoff; quantity not preserved; energy not preserved; enstrophy not preserved; health finite state fields.
- alpha c3_19790109_0000: mass preserved to roundoff; quantity not preserved; energy not preserved; enstrophy not preserved; health finite state fields.
