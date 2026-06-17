# TC4 Conservation Diagnostics

**Experiment:** Forced nonlinear exact solution

## Suggested Section Copy

Mass is preserved to roundoff (max |rel| 7.640e-16). Because TC4 is analytically forced, the mechanical-energy and derived potential-enstrophy curves are diagnostic traces, not unforced conservation claims (energy max |rel| 1.080e-04, enstrophy max |rel| 4.617e-06). A forcing-budget residual is still needed for final validation.

These diagnostics use the current MITgcm MDS outputs only. Mass and free-surface drift come from `Eta/ETAN`; energy uses `momKE` when present, otherwise centered `U/V`; PV and potential enstrophy are derived from `U/V/Eta` because no run currently writes a native vorticity or PV diagnostic.

## Note

TC4 is a forced case with completed run output. Mass can be audited from the saved free-surface fields, but energy and enstrophy are forced-response diagnostics until a forcing-budget residual is archived.

## Availability

- Configured diagnostics: UVEL, VVEL, ETAN, ETANSQ, momKE, dynDiag, WVELMASS, PhiVEL, PsiVEL, dyn_Aux, D, dynStDiag
- Runs with conservation tables: 1

## Assets To Link

- alpha run_u0_20: Mass/free surface drift -> `assets/williamson/TestCase4/Diagnosis/conservation/alpha_run_u0_20/conservation_mass.png`
- alpha run_u0_20: Energy proxy drift -> `assets/williamson/TestCase4/Diagnosis/conservation/alpha_run_u0_20/conservation_energy.png`
- alpha run_u0_20: PV/enstrophy drift -> `assets/williamson/TestCase4/Diagnosis/conservation/alpha_run_u0_20/conservation_pv_enstrophy.png`
- alpha run_u0_20: Conservation time series table -> `assets/williamson/TestCase4/Diagnosis/conservation/alpha_run_u0_20/conservation_timeseries.csv`
- alpha run_u0_20: Conservation summary data -> `assets/williamson/TestCase4/Diagnosis/conservation/alpha_run_u0_20/conservation_summary.json`

## Per-Run Verdicts

- alpha run_u0_20: mass preserved to roundoff; quantity small drift; energy noticeable drift; enstrophy small drift; health finite state fields.
