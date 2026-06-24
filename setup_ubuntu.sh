#!/bin/bash
# ============================================================================
# CO2 + SW68 Non-Equilibrium Dissolution Model — Ubuntu 22.04 Setup Script
#
# Usage:
#   chmod +x setup_ubuntu.sh
#   ./setup_ubuntu.sh
#
# This installs dependencies, compiles both steps, and runs the full pipeline.
# ============================================================================

set -e

echo "=== CO2 Oil Non-Equilibrium Model — Ubuntu Setup ==="
echo ""

# ── 1. Install system packages ──────────────────────────────────────────────
echo "[1/4] Installing system dependencies..."
sudo apt update -qq
sudo apt install -y -qq build-essential g++ cmake python3-pip python3-tk

# ── 2. Install Python packages ──────────────────────────────────────────────
echo "[2/4] Installing Python packages..."
pip3 install --quiet numpy pandas matplotlib scipy

# ── 3. Compile C++ programs ─────────────────────────────────────────────────
echo "[3/4] Compiling C++ programs..."

g++ -std=c++17 -O2 -I lib src/calc_weq_lut.cpp src/saft_mie.cpp -o calc_weq_lut
echo "  -> calc_weq_lut  compiled"

g++ -std=c++17 -O2 -I lib src/calc_taums.cpp src/saft_mie.cpp -o calc_taums
echo "  -> calc_taums   compiled"

g++ -std=c++17 -O2 -I lib src/test.cpp src/saft_mie.cpp -o test_saft
echo "  -> test_saft    compiled"

# ── 4. Run the pipeline ─────────────────────────────────────────────────────
echo "[4/4] Running computation pipeline..."
echo ""

echo "  Step 1: Equilibrium solubility LUT..."
./calc_weq_lut

echo ""
echo "  Step 2: Mesoscopic transport tau_ms LUT..."
./calc_taums

echo ""
echo "  Generating figures..."
python3 scripts/plot_weq.py
python3 scripts/plot_taums.py

echo ""
echo "============================================================"
echo "  Setup complete!"
echo ""
echo "  Output files:"
echo "    weq_lut.csv              — solubility LUT"
echo "    taums_lut.csv            — transport LUT"
echo "    fig_*.png                — visualisation figures"
echo "    docs/FIGURES_README.md   — figure-by-figure analysis"
echo "    docs/TAUMS_ANALYSIS.md   — tau_ms result analysis"
echo ""
echo "  To re-run after changing parameters:"
echo "    ./calc_weq_lut && ./calc_taums"
echo "============================================================"
