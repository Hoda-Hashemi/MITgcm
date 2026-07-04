## Conservation diagnostics

These diagnostics are computed from MITgcm MDS output. `Eta/ETAN` supplies mass and free-surface volume, `S/SALT` supplies the TC1 transported tracer when present, `momKE` or centered `U/V` supplies the mechanical-energy proxy, and `U/V/Eta` supply derived vorticity, potential vorticity, and potential enstrophy.

Derived potential enstrophy is reconstructed in postprocessing; it is not a native MITgcm output field in these runs.

$$h = H + \eta,\qquad M = \rho_0\int_\Omega h\,dA$$

$$E_m = \rho_0\int_\Omega hK\,dA + \frac{1}{2}\rho_0 g\int_\Omega \eta^2\,dA$$

$$\zeta = \frac{1}{a\cos\theta}\left(\frac{\partial v}{\partial\lambda} - \frac{\partial(u\cos\theta)}{\partial\theta}\right)$$

$$q = \frac{\zeta + f}{h},\qquad Z = \frac{1}{2}\int_\Omega hq^2\,dA$$

$$\Delta_r X(t) = \frac{X(t)-X(0)}{|X(0)|}$$

`K` is native `momKE` when available, otherwise `(u^2+v^2)/2`. `f` comes from `fCoriC.bin` when a rotated run writes it; otherwise the diagnostic uses `2*omega*sin(latitude)`.

| Case | Status | Mass | Target quantity | Energy proxy | PV/enstrophy | Output health |
| --- | --- | --- | --- | --- | --- | --- |
| TC1 | Verified | preserved to roundoff | well preserved | small drift | small drift | finite saved state fields |
| TC2 | Pending validation | preserved to roundoff | not preserved | not preserved | not preserved | finite saved state fields |
| TC3 | Pending validation | preserved to roundoff | small drift | small drift | small drift | finite saved state fields |
| TC4 | Verified | preserved to roundoff | noticeable drift | noticeable drift | noticeable drift | finite saved state fields |
| TC5 | Pending validation | preserved to roundoff | not preserved | not preserved | not preserved | finite saved state fields |
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

Mass is preserved to roundoff (max |rel| 6.699e-16). The mechanical-energy proxy is small drift (max |rel| 6.415e-07). Derived potential enstrophy is small drift (max |rel| 5.460e-07).

- alpha `0`: mass preserved to roundoff; quantity small drift; energy small drift; PV/enstrophy small drift

#### TC4

Mass is preserved to roundoff (max |rel| 7.406e-16). Because TC4 is analytically forced, the mechanical-energy and derived potential-enstrophy curves are diagnostic traces, not unforced conservation claims (energy max |rel| 1.511e-04, enstrophy max |rel| 1.360e-05). Together with the analytic path check and published snapshots, the c1 output is marked verified.

- alpha `run_u0_40`: mass preserved to roundoff; quantity noticeable drift; energy noticeable drift; PV/enstrophy noticeable drift

#### TC5

Mass is preserved to roundoff (max |rel| 1.786e-16). The mechanical-energy proxy is not preserved (max |rel| 0.444233). Derived potential enstrophy is not preserved (max |rel| 0.347357).

- alpha `0`: mass preserved to roundoff; quantity not preserved; energy not preserved; PV/enstrophy not preserved

#### TC6

Mass is preserved to roundoff (max |rel| 8.433e-16). The mechanical-energy proxy is not preserved (max |rel| 0.00122662). Derived potential enstrophy is noticeable drift (max |rel| 4.695e-04).

- alpha `0`: mass preserved to roundoff; quantity noticeable drift; energy not preserved; PV/enstrophy noticeable drift

#### TC7

Unavailable: pending validation: TC7 analyzed input is staged, but no completed MDS output is archived yet

- alpha `unavailable`: unavailable: pending validation: TC7 analyzed input is staged, but no completed MDS output is archived yet
