# TC4 Conservation Diagnostics

**Experiment:** Forced nonlinear exact solution

## Suggested Section Copy

Unavailable: scaffold / pending validation: no validated TC4 MDS run is available; the analytic forcing hook must be implemented before this can be verified

These diagnostics use the current MITgcm MDS outputs only. Mass and free-surface drift come from `Eta/ETAN`; energy uses `momKE` when present, otherwise centered `U/V`; PV and potential enstrophy are derived from `U/V/Eta` because no run currently writes a native vorticity or PV diagnostic.

## Note

TC4 is a forced case, so an unforced energy/enstrophy conservation verdict needs run output plus forcing-budget diagnostics.

## Availability

- Configured diagnostics: UVEL, VVEL, ETAN, ETANSQ, momKE, dynDiag, WVELMASS, PhiVEL, PsiVEL, dyn_Aux, D, dynStDiag
- Runs with conservation tables: 0

## Assets To Link

- No conservation assets are available for this case yet.

## Per-Run Verdicts

- alpha unavailable: unavailable (scaffold / pending validation: no validated TC4 MDS run is available; the analytic forcing hook must be implemented before this can be verified).
