# TC4 Conservation Diagnostics

**Experiment:** Forced nonlinear exact solution

## Suggested Section Copy

Unavailable: pending validation: the TC4 analytic forcing hook is prepared, but no completed MDS output is archived yet

These diagnostics use the current MITgcm MDS outputs only. Mass and free-surface drift come from `Eta/ETAN`; energy uses `momKE` when present, otherwise centered `U/V`; PV and potential enstrophy are derived from `U/V/Eta` because no run currently writes a native vorticity or PV diagnostic.

## Note

TC4 is a forced case, so an unforced energy/enstrophy conservation verdict needs run output plus forcing-budget diagnostics.

## Availability

- Configured diagnostics: UVEL, VVEL, ETAN, ETANSQ, momKE, dynDiag, WVELMASS, PhiVEL, PsiVEL, dyn_Aux, D, dynStDiag
- Runs with conservation tables: 0

## Assets To Link

- No conservation assets are available for this case yet.

## Per-Run Verdicts

- alpha unavailable: unavailable (pending validation: the TC4 analytic forcing hook is prepared, but no completed MDS output is archived yet).
