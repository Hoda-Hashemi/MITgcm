# TC7 Conservation Diagnostics

**Experiment:** Analyzed 500 mb initial state

## Suggested Section Copy

Unavailable: pending validation: TC7 analyzed input is staged, but no completed MDS output is archived yet

These diagnostics use the current MITgcm MDS outputs only. Mass and free-surface drift come from `Eta/ETAN`; energy uses `momKE` when present, otherwise centered `U/V`; PV and potential enstrophy are derived from `U/V/Eta` because no run currently writes a native vorticity or PV diagnostic.

## Note

TC7 has staged analyzed input and needs completed run output before conservation diagnostics can be evaluated.

## Availability

- Configured diagnostics: UVEL, VVEL, ETAN, ETANSQ, momKE, dynDiag, WVELMASS, PhiVEL, PsiVEL, dyn_Aux, D, dynStDiag
- Runs with conservation tables: 0

## Assets To Link

- No conservation assets are available for this case yet.

## Per-Run Verdicts

- alpha unavailable: unavailable (pending validation: TC7 analyzed input is staged, but no completed MDS output is archived yet).
