## Conservation diagnostics

These diagnostics are computed from MITgcm MDS output. `Eta/ETAN` supplies mass and free-surface volume, `S/SALT` supplies the TC1 transported tracer when present, `momKE` or centered `U/V` supplies the mechanical-energy proxy, and `U/V/Eta` supply derived vorticity, potential vorticity, and potential enstrophy.

| Case | Status | Mass | Target quantity | Energy proxy | PV/enstrophy | Output health |
| --- | --- | --- | --- | --- | --- | --- |
| TC1 | Verified | preserved to roundoff | well preserved | small drift | small drift | finite saved state fields |
| TC2 | Pending validation | preserved to roundoff | not preserved | not preserved | not preserved | finite saved state fields |
| TC3 | Pending validation | preserved to roundoff | not preserved | not preserved | not preserved | finite saved state fields |
| TC4 | Pending output | unavailable | unavailable | unavailable | unavailable | pending validation: the TC4 analytic forcing hook is prepared, but no completed MDS output is archived yet |
| TC5 | Issues | invalid output | invalid output | invalid output | invalid output | invalid; first bad day 1 |
| TC6 | Pending validation | preserved to roundoff | noticeable drift | not preserved | noticeable drift | finite saved state fields |
| TC7 | Pending output | unavailable | unavailable | unavailable | unavailable | pending validation: TC7 analyzed input is staged, but no completed MDS output is archived yet |

### Required output fields

- `Eta` or `ETAN`: mass, free-surface volume, surface potential energy.
- `U`/`V` or `UVEL`/`VVEL`: velocity health, derived KE if `momKE` is absent, relative vorticity, PV, and potential enstrophy.
- `momKE`: preferred kinetic-energy diagnostic when written.
- `S` or `SALT`: TC1 tracer amount.

### Per-run notes

#### TC1

Mass is preserved to roundoff (max |rel| 7.709e-16). The TC1 tracer amount is well preserved (max |rel| 6.213e-08). The energy curve is diagnostic only for this prescribed-flow tracer test.

- alpha `0`: mass preserved to roundoff; quantity preserved to roundoff; energy preserved to roundoff; PV/enstrophy preserved to roundoff
- alpha `0.05`: mass preserved to roundoff; quantity well preserved; energy small drift; PV/enstrophy well preserved
- alpha `1.52`: mass preserved to roundoff; quantity well preserved; energy small drift; PV/enstrophy small drift
- alpha `1.57`: mass preserved to roundoff; quantity well preserved; energy small drift; PV/enstrophy small drift

#### TC2

Mass is preserved to roundoff (max |rel| 8.499e-16). The mechanical-energy proxy is not preserved (max |rel| 0.147859). Derived potential enstrophy is not preserved (max |rel| 0.0220451).

- alpha `0`: mass preserved to roundoff; quantity small drift; energy small drift; PV/enstrophy small drift
- alpha `0.05`: mass preserved to roundoff; quantity not preserved; energy noticeable drift; PV/enstrophy not preserved
- alpha `1.52`: mass preserved to roundoff; quantity not preserved; energy not preserved; PV/enstrophy not preserved
- alpha `1.57`: mass preserved to roundoff; quantity not preserved; energy not preserved; PV/enstrophy not preserved

#### TC3

Mass is preserved to roundoff (max |rel| 1.172e-15). The mechanical-energy proxy is not preserved (max |rel| 0.0830516). Derived potential enstrophy is not preserved (max |rel| 0.0365314).

- alpha `0`: mass preserved to roundoff; quantity small drift; energy small drift; PV/enstrophy small drift
- alpha `1.0472`: mass preserved to roundoff; quantity not preserved; energy not preserved; PV/enstrophy not preserved

#### TC4

Unavailable: pending validation: the TC4 analytic forcing hook is prepared, but no completed MDS output is archived yet

- alpha `unavailable`: unavailable: pending validation: the TC4 analytic forcing hook is prepared, but no completed MDS output is archived yet

#### TC5

Output health is invalid because at least one Eta/U/V/scalar state contains non-finite values (minimum finite fraction 0). The first bad saved record is day 1. Do not use the conservation verdict as final until the run is regenerated with finite state fields.

- alpha `0`: invalid output; 15 non-finite saved records, first bad day 1

#### TC6

Mass is preserved to roundoff (max |rel| 8.433e-16). The mechanical-energy proxy is not preserved (max |rel| 0.00122662). Derived potential enstrophy is noticeable drift (max |rel| 4.695e-04).

- alpha `0`: mass preserved to roundoff; quantity noticeable drift; energy not preserved; PV/enstrophy noticeable drift

#### TC7

Unavailable: pending validation: TC7 analyzed input is staged, but no completed MDS output is archived yet

- alpha `unavailable`: unavailable: pending validation: TC7 analyzed input is staged, but no completed MDS output is archived yet
