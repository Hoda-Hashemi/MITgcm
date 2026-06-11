# 🐍 Python/ 
README for the `Python/` folder.

### Folder structure 📂

```text
Python/
│
├── data/                           #  ParaView exports
│   └── ParaviewValues.csv          # CSV exported from ParaView
│   └── interpolatedDepthParaview.csv 
│
├── output/                         # Organized simulation outputs
│   └── run_YYYYMMDD_HHMMSS/        # Timestamped run directories
│       ├── plots.png               # Example plot(s)
│       └── log.txt                 # Terminal or kernel outputs
│
├── src/                            # Source code and reports
│   ├── __pycache__/                
│   ├── CheckInterpolationPlots.py 
│   ├── gridOperations.py           
│   ├── config.py                   
│   ├── plotting.py                 
│   ├── SphereSolver_1D.py          
│   ├── SphereSolver_2D.py         
│   ├── utils.py                    
│   └── OceanmainVortexSphere.py    # Main simulation entry point
│
├── README.md                      
└── requirements.txt                # Python dependencies
```

---
### Quick start ⚡

1. **Create & activate** a virtual environment:

```bash
python -m venv venv
# macOS / Linux
source venv/bin/activate
```

2. **Install dependencies** (includes `ipykernel` for running cells on vscode, shift+enter):

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

---
### Data input for comparison (ParaView) 🌏
* To export data from ParaView: open your dataset, then `File -> Save Data...` and choose CSV. Place the exported CSV inside `Python/data/` — e.g. `data/ParaviewValues.csv`.
* The code expects CSVs to contain named columns used by `src/utils.py` / `src/plotting.py`. If you change column names, update callers accordingly.

---

### Output
Each simulation run creates a timestamped folder under `output/` (e.g. `run_20251031_131001/`). Typical run contents:

* `compare_side_1D.png`, `Plots_<src.SphereSolver_1D.SphereSolver1D object at 0x11514fb60>>.png`, ... — visualization PNGs
* `log.txt` — run log and printed diagnostics

### Latest run — visual summary

<p align="center">
  <img src="output/run_20251031_151639/Plots_<src.SphereSolver_2D.SphereSolver2D object at 0x114d1b8c0>.png" alt="Main plot" width="500" />
</p>

<!-- Subfigures table -->
<table align="center">
  <tr>
    <td align="center">
      <img src="output/run_20251031_151639/compare_overlay_1D.png" width="500" alt="1D overlay" /><br>
      <em>1D overlay</em>
    </td>
    <td align="center">
      <img src="output/run_20251031_151639/compare_side_1D.png" width="500" alt="1D side" /><br>
      <em>1D side-by-side</em>
    </td>
    </tr>
</table>

<table align="center">
  <tr>
    <td align="center">
      <img src="output/run_20251031_151639/compare_overlay_2D.png" width="500" alt="2D overlay" /><br>
      <em>2D overlay</em>
    </td>
    <td align="center">
      <img src="output/run_20251031_151639/compare_side_2D.png" width="500" alt="2D side" /><br>
      <em>2D side-by-side</em>
    </td>
  </tr>
</table>

---

### Configuration 👾
* To change parameter values persistently, edit `src/config.py`. 

---
### Example workflow

1. Prepare or export input data to `data/`.
2. Change parameters in `src/config.py` as needed.
3. Run `python src/OceanmainVortexSphere.py` or open the Jupyter notebook.
4. Inspect `output/run_<timestamp>/` for `plots.png` and `log.txt`.

---

### Recommended ⁕

* Use the VS Code ipykernel for an optimal workflow: it allows you to run, validate, and inspect outputs directly in the same window.
* Keep `config.py` minimal: avoid hardcoding absolute paths. Use relative paths and create run directories using a helper function.

---

### Troubleshooting

* If dependencies fail to install: upgrade `pip` and retry: `python -m pip install --upgrade pip` then `pip install -r requirements.txt`.
* If plots are empty or constant: verify that input CSV column names match expected names and that `config.py` grid sizes are set correctly.

---

