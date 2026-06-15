# Submission Readiness

Audit date: 2026-06-16.

| Case | Tonight action | Status |
| --- | --- | --- |
| TC1 | Submit/rerun after optional cleanup | Ready; current run ends normally, but TC1 alpha `1.57` has stale extra iterations/assets |
| TC2 | Archive stale run dirs, then submit | Templates fixed; current run dirs lack `selectCoriMap=3` and `fCori*.bin` |
| TC3 | Archive stale run dirs, then submit | Templates fixed; current run dirs lack `selectCoriMap=3` and `fCori*.bin` |
| TC4 | Do not submit | Forced case is still scaffold/pending validation |
| TC5 | Debug-only | Current fields are finite at iter `0000000000`, then non-finite from iter `0000001440` onward |
| TC6 | Submit/rerun | Ready; final fields finite through iter `0000040320` |
| TC7 | Do not submit | Missing `input/raw/tc7_initial_conditions.npz` |

Final-submit commands tonight: TC1, TC2 after archive/restage, TC3 after archive/restage, and TC6. No final submission for TC4, TC5, or TC7.

