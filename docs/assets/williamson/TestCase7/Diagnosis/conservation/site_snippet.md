# TC7 Conservation Diagnostics

**Experiment:** Analyzed 500 mb initial state

## Suggested Section Copy

Unavailable: needs data: TC7 requires the analyzed 500 mb initial-condition NPZ and a completed MITgcm run before conservation can be checked

These diagnostics use the current MITgcm MDS outputs only. Mass and free-surface drift come from `Eta/ETAN`; energy uses `momKE` when present, otherwise centered `U/V`; PV and potential enstrophy are derived from `U/V/Eta` because no run currently writes a native vorticity or PV diagnostic.

## Note

TC7 needs a completed run before conservation diagnostics can be evaluated.

## Availability

- Configured diagnostics: UVEL, VVEL, ETAN, ETANSQ, momKE, dynDiag, WVELMASS, PhiVEL, PsiVEL, dyn_Aux, D, dynStDiag
- Runs with conservation tables: 0

## Assets To Link

- No conservation assets are available for this case yet.

## Per-Run Verdicts

- alpha unavailable: unavailable (needs data: TC7 requires the analyzed 500 mb initial-condition NPZ and a completed MITgcm run before conservation can be checked).
