# TC2 Conservation Diagnostics

**Experiment:** Steady solid-body geostrophic flow

## Suggested Section Copy

Mass is preserved to roundoff (max |rel| 8.499e-16). The mechanical-energy proxy is not preserved (max |rel| 0.147859). Derived potential enstrophy is not preserved (max |rel| 0.0220451).

These diagnostics use the current MITgcm MDS outputs only. Mass and free-surface drift come from `Eta/ETAN`; energy uses `momKE` when present, otherwise centered `U/V`; PV and potential enstrophy are derived from `U/V/Eta` because no run currently writes a native vorticity or PV diagnostic.

## Availability

- Configured diagnostics: UVEL, VVEL, ETAN, ETANSQ, momKE, dynDiag, WVELMASS, PhiVEL, PsiVEL, dyn_Aux, D, dynStDiag
- Runs with conservation tables: 4

## Assets To Link

- alpha 0: Conservation time series table -> `assets/williamson/TestCase2/Diagnosis/conservation/alpha_0/conservation_timeseries.csv`
- alpha 0: Conservation summary data -> `assets/williamson/TestCase2/Diagnosis/conservation/alpha_0/conservation_summary.json`
- alpha 0.05: Conservation time series table -> `assets/williamson/TestCase2/Diagnosis/conservation/alpha_0.05/conservation_timeseries.csv`
- alpha 0.05: Conservation summary data -> `assets/williamson/TestCase2/Diagnosis/conservation/alpha_0.05/conservation_summary.json`
- alpha 1.52: Conservation time series table -> `assets/williamson/TestCase2/Diagnosis/conservation/alpha_1.52/conservation_timeseries.csv`
- alpha 1.52: Conservation summary data -> `assets/williamson/TestCase2/Diagnosis/conservation/alpha_1.52/conservation_summary.json`
- alpha 1.57: Conservation time series table -> `assets/williamson/TestCase2/Diagnosis/conservation/alpha_1.57/conservation_timeseries.csv`
- alpha 1.57: Conservation summary data -> `assets/williamson/TestCase2/Diagnosis/conservation/alpha_1.57/conservation_summary.json`

## Per-Run Verdicts

- alpha 0: mass preserved to roundoff; quantity small drift; energy small drift; enstrophy small drift; health finite state fields.
- alpha 0.05: mass preserved to roundoff; quantity not preserved; energy noticeable drift; enstrophy not preserved; health finite state fields.
- alpha 1.52: mass preserved to roundoff; quantity not preserved; energy not preserved; enstrophy not preserved; health finite state fields.
- alpha 1.57: mass preserved to roundoff; quantity not preserved; energy not preserved; enstrophy not preserved; health finite state fields.
