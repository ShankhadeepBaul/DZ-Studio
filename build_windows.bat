@echo off
REM ===================================================================
REM  DZ Studio - one-time build script (WINDOWS)
REM
REM  Run this ONCE on a Windows PC that has Python 3.10+ installed
REM  (tick "Add Python to PATH" when you install Python).
REM
REM  It produces:  dist\DZ Studio\  -- a self-contained folder.
REM  Zip that folder and give it to anyone. They double-click
REM  "DZ Studio.exe" inside it. They do NOT need Python, or anything.
REM
REM  You can only build a WINDOWS app on Windows. For a Mac app,
REM  run build_macos.sh on a Mac.
REM ===================================================================
echo.
echo [1/3] Installing build tools (once)...
python -m pip install --upgrade pip
python -m pip install -r requirements.txt pyinstaller
if errorlevel 1 goto :fail

echo.
echo [2/3] Building the app (this takes a few minutes)...
pyinstaller --noconfirm --windowed --onedir ^
    --name "DZ Studio" ^
    --add-data "example_data.xlsx;." ^
    --hidden-import openpyxl ^
    --hidden-import openpyxl.cell._writer ^
    --hidden-import matplotlib.backends.backend_pdf ^
    --hidden-import matplotlib.backends.backend_svg ^
    --hidden-import matplotlib.backends.backend_agg ^
    --exclude-module PyQt5 ^
    --exclude-module PyQt6 ^
    --exclude-module tkinter ^
    dz_app.py
if errorlevel 1 goto :fail

echo.
echo [3/3] Done.
echo.
echo   Your app is here:   dist\DZ Studio\
echo   Run it by opening:  dist\DZ Studio\DZ Studio.exe
echo.
echo   To share: right-click the "DZ Studio" folder, Send to -^> Compressed
echo   (zipped) folder, and send that .zip. The other person unzips and
echo   double-clicks DZ Studio.exe.
echo.
pause
exit /b 0

:fail
echo.
echo BUILD FAILED. Copy the red error text above and send it for help.
echo Most common cause: Python not installed or not on PATH.
pause
exit /b 1
