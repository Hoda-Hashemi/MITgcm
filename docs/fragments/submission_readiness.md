# Submission Readiness

Audit date: 2026-06-17.

| Case | Tonight action | Status |
| --- | --- | --- |
| TC1 | Submit/rerun after optional cleanup | Ready; current run ends normally, but TC1 alpha `1.57` has stale extra iterations/assets |
| TC2 | Archive stale run dirs, then submit | Templates fixed; current run dirs lack `selectCoriMap=3` and `fCori*.bin` |
| TC3 | Archive stale run dirs, then submit | Templates fixed; current run dirs lack `selectCoriMap=3` and `fCori*.bin` |
| TC4 | Completed on normal/onode08 | `run_u0_20` completed through day 5; fields are finite and the translating low follows the analytic path. Formal forcing-budget diagnostics remain pending |
| TC5 | Debug-only | Current fields are finite at iter `0000000000`, then non-finite from iter `0000001440` onward |
| TC6 | Submit/rerun | Ready; final fields finite through iter `0000040320` |
| TC7 | Do not submit | Missing `input/raw/tc7_initial_conditions.npz` |

TC4 is no longer waiting on output, but still should not be called finally validated until formal forcing-budget diagnostics are archived. No final submission for TC5 or TC7.
