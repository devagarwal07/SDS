# SDS Streamlit App

This repo contains a Streamlit app that generates Safety Data Sheets (SDS) from a SMILES string.

## Deploy on Streamlit Cloud (recommended)

Streamlit Cloud build files:
- `requirements.txt` pins compatible versions and installs RDKit on Linux/macOS via `rdkit-pypi`.
- `runtime.txt` pins Python to 3.10.13 so RDKit wheels resolve.

Steps:
1. Push this repo to GitHub.
2. In Streamlit Cloud, create an app:
   - Repo: this repo
   - Branch: main
   - App file: `SDS/app.py`
3. Deploy.

## Run locally on Windows with RDKit

On Windows, install RDKit via conda (no official pip wheel). Use the provided environment file:

```powershell
# Install Miniconda or Mambaforge first
conda env create -f SDS/environment.yml
conda activate sds
python -m streamlit run SDS/app.py
```

## Run locally without conda (no RDKit)

If you prefer plain pip (no RDKit locally):

```powershell
python -m pip install -r SDS/requirements.txt
python -m streamlit run SDS/app.py
```

Note: Some features (like RDKit-derived properties) require RDKit; they’ll work fully on Streamlit Cloud.
