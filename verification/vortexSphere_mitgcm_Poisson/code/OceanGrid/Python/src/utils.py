import pandas as pd

#!COMMENT:  Load a CSV file (saved data from paraview), clean backspace characters, and return a DataFrame. 
import os
import pandas as pd
import tempfile

def ParaviewSavedData(file_path: str) -> pd.DataFrame:
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File '{file_path}' not found.")
    
    # Read file, remove backspace chars
    with open(file_path, "r", errors="ignore") as f:
        content = f.read().replace("\x08", "")
    
    # Write to secure temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False, encoding='utf-8') as tmp:
        tmp.write(content)
        temp_path = tmp.name
    
    try:
        df = pd.read_csv(temp_path)

        df.columns = [
        c.replace('ψ₁', 'psi')
        .replace('ψ', 'psi')
        .replace('d₀', 'depth')
        .replace('θ', 'theta')
        .replace('lambda', 'phi')
        .replace('ƛ', 'phi')
        .replace('λ', 'phi')
        .replace('f(theta)', 'Coriolis')
        .replace('f(θ)', 'Coriolis')
        .replace('u₁(theta,phi):0', 'u0')
        .replace('u₁(θ,ƛ):0', 'u0')
        .replace('u₁(theta,phi):1', 'u1')
        .replace('u₁(θ,ƛ):1', 'u1')
        .replace('u₁(theta,phi):2', 'u2')
        .replace('u₁(θ,ƛ):2', 'u2')
        .replace('Cell Type', 'cell_type')
        .strip()
        for c in df.columns
    ]

        # Optional: strip whitespace
        df.columns = df.columns.str.strip()
        
    finally:
        os.unlink(temp_path)  # Always clean up
    
    return df

 #!ROCKET: Clean Logging 
import os
import sys
from datetime import datetime

def setup_logging(base_dir="output"):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(base_dir, f"run_{timestamp}")
    os.makedirs(run_dir, exist_ok=True)

    log_path = os.path.join(run_dir, "log.txt")

    class DualWriter:
        def __init__(self, filename):
            self.terminal = sys.stdout
            self.log = open(filename, "w", encoding="utf-8")

        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)

        def flush(self):
            self.terminal.flush()
            self.log.flush()

    sys.stdout = DualWriter(log_path)
    print(f"\nLogging to {log_path}\n")

    return run_dir

