# TC5 Conservation Diagnostics

**Experiment:** Zonal flow over an isolated mountain

## Suggested Section Copy

Output health is invalid because at least one Eta/U/V/scalar state contains non-finite values (minimum finite fraction 0). The first bad saved record is day 1. Do not use the conservation verdict as final until the run is regenerated with finite state fields.

These diagnostics use the current MITgcm MDS outputs only. Mass and free-surface drift come from `Eta/ETAN`; energy uses `momKE` when present, otherwise centered `U/V`; PV and potential enstrophy are derived from `U/V/Eta` because no run currently writes a native vorticity or PV diagnostic.

## Availability

- Configured diagnostics: UVEL, VVEL, ETAN, ETANSQ, momKE, dynDiag, WVELMASS, PhiVEL, PsiVEL, dyn_Aux, D, dynStDiag
- Runs with conservation tables: 1

## Assets To Link

- alpha 0: Mass/free surface drift -> `assets/williamson/TestCase5/Diagnosis/conservation/alpha_0/conservation_mass.png`
- alpha 0: Energy proxy drift -> `assets/williamson/TestCase5/Diagnosis/conservation/alpha_0/conservation_energy.png`
- alpha 0: PV/enstrophy drift -> `assets/williamson/TestCase5/Diagnosis/conservation/alpha_0/conservation_pv_enstrophy.png`
- alpha 0: Conservation time series table -> `assets/williamson/TestCase5/Diagnosis/conservation/alpha_0/conservation_timeseries.csv`
- alpha 0: Conservation summary data -> `assets/williamson/TestCase5/Diagnosis/conservation/alpha_0/conservation_summary.json`

## Per-Run Verdicts

- alpha 0: mass invalid output; quantity invalid output; energy invalid output; enstrophy invalid output; health invalid output: non-finite state fields; first non-finite state at day 1.
