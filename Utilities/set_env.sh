#!/bin/bash

echo ""
echo ">>> Setting necessary environment variables"
echo "==========================================="

export DATA_DIR=${HOME}/analysis-omega/out/GetSimpleTuple/data
export SIM_DIR=${HOME}/analysis-omega/out/GetSimpleTuple/sim
export HAPRAD_CPP=${HOME}/software/HAPRAD_CPP
# CERNLIB is NOT required. Set CERN_LIB and build with `make USE_CERNLIB=1`
# only if you restore the PDF x FF semi-inclusive model (which needs PDFLIB's
# PDFSET/STRUCTM via init_pdf.f). Example for Debian/Ubuntu:
# export CERN_LIB=/lib/x86_64-linux-gnu/

echo ">>> DATA_DIR   = ${DATA_DIR}"
echo ">>> SIM_DIR    = ${SIM_DIR}"
echo ">>> HAPRAD_CPP = ${HAPRAD_CPP}"
echo ">>> CERN_LIB   = ${CERN_LIB:-<unset, not required>}"
echo ""
