@echo off
REM Compile einverted with G-U wobble patch for Windows
REM Requires MSYS2 or Cygwin with gcc, make, patch, curl/wget

echo Compiling einverted with G-U wobble patch for Windows...
echo.

REM Check if we're in MSYS2 or Cygwin environment
where bash >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: bash not found. Please install MSYS2 or Cygwin and add to PATH.
    echo.
    echo To install MSYS2:
    echo   1. Download from https://www.msys2.org/
    echo   2. Install and run MSYS2 MinGW 64-bit terminal
    echo   3. Install build tools: pacman -S mingw-w64-x86_64-gcc make patch curl tar
    echo   4. Run this script from MSYS2 terminal
    exit /b 1
)

REM Run the actual compilation in bash
bash compile_einverted_windows.sh %*