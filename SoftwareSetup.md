# Softwares

---

## MATLAB

**Tested on:** R2023b
**Installation:** [https://www.mathworks.com/downloads/](https://www.mathworks.com/downloads/)

### macOS

Add to PATH:

```bash
echo 'export PATH="/Applications/MATLAB_R2023b.app/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

```bash
matlab -batch "disp(version)"
# 9.14.0... (R2023b)
```

### Linux

Add to PATH:

```bash
echo 'export PATH="/usr/local/MATLAB/R2023b/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

```bash
matlab -batch "disp(version)"
# 9.14.0... (R2023b)
```

### Windows

Added to PATH automatically by the installer.

```powershell
matlab -batch "disp(version)"
# 9.14.0... (R2023b)
```

## Gurobi

**Tested on:** 13.0.2
**Installation:** [https://www.gurobi.com/downloads/](https://www.gurobi.com/downloads/)
**License:** free academic license — register on the Gurobi portal, then run the `grbgetkey` command it gives you.

### macOS

Add to PATH:

```bash
echo 'export GUROBI_HOME="/Library/gurobi1302/macos_universal2"' >> ~/.zshrc
echo 'export PATH="$GUROBI_HOME/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

Activate license:

```bash
grbgetkey YOUR-KEY-HERE
```

Link to MATLAB:

```matlab
addpath('/Library/gurobi1302/macos_universal2/matlab')
savepath
```

```bash
gurobi_cl --version
# Gurobi 13.0.2 (mac64[arm], academic)
```

### Linux

```bash
sudo tar -xzf gurobi13.0.2_linux64.tar.gz -C /opt/
echo 'export GUROBI_HOME="/opt/gurobi1302/linux64"' >> ~/.bashrc
echo 'export PATH="$GUROBI_HOME/bin:$PATH"' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH="$GUROBI_HOME/lib:$LD_LIBRARY_PATH"' >> ~/.bashrc
source ~/.bashrc
```

Activate license:

```bash
grbgetkey YOUR-KEY-HERE
```

Link to MATLAB:

```matlab
addpath('/opt/gurobi1302/linux64/matlab')
savepath
```

```bash
gurobi_cl --version
# Gurobi 13.0.2 (linux64, academic)
```

### Windows

Added to PATH automatically by the installer — open a new PowerShell window afterward.

Activate license:

```powershell
grbgetkey YOUR-KEY-HERE
```

Link to MATLAB:

```matlab
addpath('C:\gurobi1302\win64\matlab')
savepath(fullfile(userpath, 'pathdef.m'))
```

```powershell
gurobi_cl --version
# Gurobi 13.0.2 (win64, academic)
```

WSL (Windows Subsystem for Linux) is required to run the pipeline's shell scripts. See [Changes](#windows-pipeline-scripts-run-under-wsl). For installing WSL, see [Python (Miniforge)](#python-miniforge) → Windows.

Do not activate a separate license inside WSL. Point it at the license file this native install already produced:

```bash
echo 'export GRB_LICENSE_FILE="/mnt/c/Users/<name>/gurobi.lic"' >> ~/.bashrc
source ~/.bashrc
```

## Python (Miniforge)

**Installation:** [https://github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge)

### macOS

```bash
bash Miniforge3-MacOSX-arm64.sh
```

```bash
conda --version
python3 --version
```

### Linux

```bash
bash Miniforge3-Linux-x86_64.sh
```

```bash
conda --version
python3 --version
```

### Windows

Install WSL:

```powershell
wsl --install -d Ubuntu
```

Then, inside WSL, follow the Linux steps above.

---

# Packages and Environments

---

## COBRA Toolbox

**Installation:** [https://github.com/opencobra/cobratoolbox](https://github.com/opencobra/cobratoolbox)

### macOS / Linux

```bash
git clone --depth=1 https://github.com/opencobra/cobratoolbox.git ~/cobratoolbox
```

```matlab
addpath('~/cobratoolbox')
savepath
initCobraToolbox(false)
```

### Windows

```powershell
git clone --depth=1 https://github.com/opencobra/cobratoolbox.git $env:USERPROFILE\cobratoolbox
```

```matlab
addpath(fullfile(getenv('USERPROFILE'), 'cobratoolbox'))
savepath(fullfile(userpath, 'pathdef.m'))
initCobraToolbox(false)
```

## Python Environment

Installed automatically by `setup_environment.sh`, along with the rest of the pipeline's Python dependencies.

### macOS / Linux

From the Pipeline's root directory:

```bash
bash setup_environment.sh
```

### Windows

From the Pipeline's root directory, inside WSL:

```bash
bash setup_environment.sh
```

---

# Troubleshooting and known issues

---

## MATLAB (Apple Silicon)

| Error | Cause | Fix |
|---|---|---|
| `TranslateSBML` / FBC extension fails | COBRA Toolbox's SBML MEX binaries have no native ARM64 build | Load the `.mat` files produced by `genome_to_draftmodels.sh`'s CobraPy step directly, instead of `readCbModel` on the SBML. |

## Gurobi (Windows)

| Error | Cause | Fix |
|---|---|---|
| `gurobi_cl` not recognized | `bin` not on PATH | Check: `$env:PATH -split ';' \| Select-String "gurobi"`. If missing, add `C:\gurobi1302\win64\bin` via System Properties → Environment Variables, then open a new PowerShell. |
| `No Gurobi license found` | `gurobi.lic` missing | Re-run `grbgetkey YOUR-KEY-HERE`. |
| `Invalid MEX-file... module could not be found` | MATLAB launched before PATH was updated | Confirm `C:\gurobi1302\win64\bin` is on PATH, then restart MATLAB. |
| `Solver not found` in MATLAB | Gurobi MATLAB folder not on path | Run `addpath('C:\gurobi1302\win64\matlab')` then `savepath(fullfile(userpath, 'pathdef.m'))`. |

## CarveMe

| Error | Cause | Fix |
|---|---|---|
| `carve`: Unable to run diamond | `diamond`'s directory isn't on PATH for CarveMe's internal subprocess call | `config.sh` exports the env's `bin` directory onto PATH — confirm it was regenerated by `setup_environment.sh`. |

---

# Changes

---

## Gurobi replaces CPLEX

The original pipeline used IBM CPLEX. CPLEX 12.10 ships an x86-only MEX file, which does not load on Apple Silicon, and its PyPI Python API is the Community Edition, size-limited and unable to solve CarveMe's problems. Gurobi runs natively on both architectures with no size limit under the academic license, and is used for both CarveMe and every MATLAB COBRA step.

## Windows: pipeline scripts run under WSL

The pipeline's scripts (`setup_environment.sh` and the stages built on it) are bash, and depend on conda and Gurobi's Python bindings. Native Windows has no bash, so this part of the setup runs inside WSL instead, following the Linux instructions. MATLAB, its Gurobi link, and COBRA Toolbox stay on native Windows since they don't depend on bash.

Gurobi itself is only installed once, natively on Windows. `gurobipy` bundles its own solver library, so WSL does not need a separate Gurobi Optimizer install — it only needs the license file that native install's `grbgetkey` produced.
