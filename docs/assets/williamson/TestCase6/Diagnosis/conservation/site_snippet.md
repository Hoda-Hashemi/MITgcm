# TC6 Conservation Diagnostics

**Experiment:** Rossby-Haurwitz wave

## Suggested Section Copy

Mass is preserved to roundoff (max |rel| 8.433e-16). The mechanical-energy proxy is not preserved (max |rel| 0.00122662). Derived potential enstrophy is noticeable drift (max |rel| 4.695e-04).

These diagnostics use the current MITgcm MDS outputs only. Mass and free-surface drift come from `Eta/ETAN`; energy uses `momKE` when present, otherwise centered `U/V`; PV and potential enstrophy are derived from `U/V/Eta` because no run currently writes a native vorticity or PV diagnostic.

## Availability

- Configured diagnostics: UVEL, VVEL, ETAN, ETANSQ, momKE, dynDiag, WVELMASS, PhiVEL, PsiVEL, dyn_Aux, D, dynStDiag
- Runs with conservation tables: 1

## Assets To Link

- alpha 0: Mass/free surface drift -> `/scratch/p8190783-hmh85/MITgcm/Sandbox/output/TestCase6/Diagnosis/conservation/alpha_0/conservation_mass.png`
- alpha 0: Energy proxy drift -> `/scratch/p8190783-hmh85/MITgcm/Sandbox/output/TestCase6/Diagnosis/conservation/alpha_0/conservation_energy.png`
- alpha 0: PV/enstrophy drift -> `/scratch/p8190783-hmh85/MITgcm/Sandbox/output/TestCase6/Diagnosis/conservation/alpha_0/conservation_pv_enstrophy.png`
- alpha 0: Conservation time series table -> `/scratch/p8190783-hmh85/MITgcm/Sandbox/output/TestCase6/Diagnosis/conservation/alpha_0/conservation_timeseries.csv`
- alpha 0: Conservation summary data -> `/scratch/p8190783-hmh85/MITgcm/Sandbox/output/TestCase6/Diagnosis/conservation/alpha_0/conservation_summary.json`

## Per-Run Verdicts

- alpha 0: mass preserved to roundoff; quantity noticeable drift; energy not preserved; enstrophy noticeable drift; health finite state fields.
