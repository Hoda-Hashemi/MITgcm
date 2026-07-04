# TC4 Conservation and Output-Health Checks

**Experiment:** Forced nonlinear exact solution

## Suggested Section Copy

Mass is preserved to roundoff (max |rel| 7.406e-16). Because TC4 is analytically forced, the mechanical-energy and derived potential-enstrophy curves are diagnostic traces, not unforced conservation claims (energy max |rel| 1.511e-04, enstrophy max |rel| 1.360e-05). Together with the analytic path check and published snapshots, the c1 output is marked verified.

These diagnostics use the current MITgcm MDS outputs only. Mass and free-surface drift come from `Eta/ETAN`; energy uses `momKE` when present, otherwise centered `U/V`; PV and potential enstrophy are derived from `U/V/Eta` because no run currently writes a native vorticity or PV diagnostic.

## Note

TC4 is a forced case with completed run output. Mass can be audited from the saved free-surface fields; energy and enstrophy are reported as forced-response diagnostics rather than unforced conservation invariants.

## Availability

- Configured diagnostics: UVEL, VVEL, ETAN, ETANSQ, momKE, dynDiag, WVELMASS, PhiVEL, PsiVEL, dyn_Aux, D, dynStDiag
- Runs with conservation tables: 1

## Assets To Link

- alpha run_u0_40: Mass/free surface drift -> `assets/williamson/TestCase4/Diagnosis/conservation/alpha_run_u0_40/conservation_mass.png`
- alpha run_u0_40: Energy proxy drift -> `assets/williamson/TestCase4/Diagnosis/conservation/alpha_run_u0_40/conservation_energy.png`
- alpha run_u0_40: PV/enstrophy drift -> `assets/williamson/TestCase4/Diagnosis/conservation/alpha_run_u0_40/conservation_pv_enstrophy.png`
- alpha run_u0_40: Conservation time series table -> `assets/williamson/TestCase4/Diagnosis/conservation/alpha_run_u0_40/conservation_timeseries.csv`
- alpha run_u0_40: Conservation summary data -> `assets/williamson/TestCase4/Diagnosis/conservation/alpha_run_u0_40/conservation_summary.json`

## Per-Run Verdicts

- alpha run_u0_40: mass preserved to roundoff; quantity noticeable drift; energy noticeable drift; enstrophy noticeable drift; health finite state fields.
