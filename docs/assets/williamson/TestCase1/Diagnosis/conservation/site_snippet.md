# TC1 Conservation Diagnostics

**Experiment:** Cosine bell passive-tracer advection

## Suggested Section Copy

Mass is preserved to roundoff (max |rel| 7.709e-16). The TC1 tracer amount is well preserved (max |rel| 6.213e-08). The energy curve is diagnostic only for this prescribed-flow tracer test.

These diagnostics use the current MITgcm MDS outputs only. Mass and free-surface drift come from `Eta/ETAN`; energy uses `momKE` when present, otherwise centered `U/V`; PV and potential enstrophy are derived from `U/V/Eta` because no run currently writes a native vorticity or PV diagnostic.

## Note

TC1 uses prescribed nondivergent velocity, so the tracer amount is the primary conserved quantity. Energy is reported only as a diagnostic proxy from the available fields.

## Availability

- Configured diagnostics: UVEL, VVEL, ETAN, ETANSQ, SALT, SALTSQ, dynDiag, WVELMASS, PhiVEL, PsiVEL, dyn_Aux, D, UVELSQ, VVELSQ, dynStDiag
- Runs with conservation tables: 4

## Assets To Link

- alpha 0: Conservation time series table -> `assets/williamson/TestCase1/Diagnosis/conservation/alpha_0/conservation_timeseries.csv`
- alpha 0: Conservation summary data -> `assets/williamson/TestCase1/Diagnosis/conservation/alpha_0/conservation_summary.json`
- alpha 0.05: Conservation time series table -> `assets/williamson/TestCase1/Diagnosis/conservation/alpha_0.05/conservation_timeseries.csv`
- alpha 0.05: Conservation summary data -> `assets/williamson/TestCase1/Diagnosis/conservation/alpha_0.05/conservation_summary.json`
- alpha 1.52: Conservation time series table -> `assets/williamson/TestCase1/Diagnosis/conservation/alpha_1.52/conservation_timeseries.csv`
- alpha 1.52: Conservation summary data -> `assets/williamson/TestCase1/Diagnosis/conservation/alpha_1.52/conservation_summary.json`
- alpha 1.57: Conservation time series table -> `assets/williamson/TestCase1/Diagnosis/conservation/alpha_1.57/conservation_timeseries.csv`
- alpha 1.57: Conservation summary data -> `assets/williamson/TestCase1/Diagnosis/conservation/alpha_1.57/conservation_summary.json`

## Per-Run Verdicts

- alpha 0: mass preserved to roundoff; quantity preserved to roundoff; energy preserved to roundoff; enstrophy preserved to roundoff; health finite state fields.
- alpha 0.05: mass preserved to roundoff; quantity well preserved; energy small drift; enstrophy well preserved; health finite state fields.
- alpha 1.52: mass preserved to roundoff; quantity well preserved; energy small drift; enstrophy small drift; health finite state fields.
- alpha 1.57: mass preserved to roundoff; quantity well preserved; energy small drift; enstrophy small drift; health finite state fields.
