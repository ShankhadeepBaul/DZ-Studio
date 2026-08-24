#!/usr/bin/env bash
# ===================================================================
#  DZ Studio - one-time build script (macOS)
#  Run once on a Mac with Python 3.10+:   bash build_macos.sh
#  Produces: dist/DZ Studio.app  -- give it to any Mac user, no Python needed.
#  (A Mac app can only be built on a Mac.)
# ===================================================================
set -e
echo "[1/3] Installing build tools (once)..."
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt pyinstaller

echo "[2/3] Building the app (a few minutes)..."
pyinstaller --noconfirm --windowed --onedir \
    --name "DZ Studio" \
    --add-data "example_data.xlsx:." \
    --hidden-import openpyxl \
    --hidden-import openpyxl.cell._writer \
    --hidden-import matplotlib.backends.backend_pdf \
    --hidden-import matplotlib.backends.backend_svg \
    --hidden-import matplotlib.backends.backend_agg \
    --exclude-module PyQt5 \
    --exclude-module PyQt6 \
    --exclude-module tkinter \
    dz_app.py

echo "[3/3] Done. Your app: dist/DZ Studio.app"
echo "Zip it (right-click -> Compress) and share. Recipients may need to"
echo "right-click -> Open the first time (unsigned app warning)."
