#!/bin/bash

#####################################################
#                                                   #
#  Script made to execute the chain of programs     #
#        to calculate hadron radiative corrections  #
#                                                   #
#####################################################

function print_help() {
    echo "SCRIPT: exec_rad-corr_chain.sh"
    echo "======================="
    echo "./exec_rad-corr_chain.sh --pid <pid>"
    echo "Where:"
    echo "  <pid>  = selects particle according to its PDG number (211, -211, 2212)"
    echo "Example:"
    echo "  ./exec_rad-corr_chain.sh --pid 211"
    exit
}

function process_args() {
    arr=("$@")
    ic=0
    while [[ $ic -le $((${#arr[@]}-1)) ]]; do
        if [[ "${arr[$ic]}" == "--pid" ]]; then
            pid=${arr[$((ic+1))]}
        else
            echo "ERROR: unrecognized argument: ${arr[$((ic))]}.";
            print_help;
        fi
        ((ic+=2))
    done
}

function print_args() {
    echo "SCRIPT: exec_rad-corr_chain.sh"
    echo "=============================="
    echo "pid = ${pid}"
}

################
###   Main   ###
################

if [[ -z "${DATA_DIR}" ]]; then
    echo "ERROR: variable DATA_DIR is unset."
    exit 1
fi

if [[ -z "${HAPRAD_CPP}" ]]; then
    echo "ERROR: variable HAPRAD_CPP is unset."
    exit 1
fi

# CERN_LIB is optional: it is only consulted when building with USE_CERNLIB=1.
if [[ -n "${CERN_LIB}" ]]; then
    echo "INFO: CERN_LIB is set (${CERN_LIB}); build with USE_CERNLIB=1 to link it."
fi

if [[ ${#} -ne 2 ]]; then
    echo "ERROR: ${#} arguments were provided, they should be 2."
    print_help
fi

argArray=("$@")
process_args "${argArray[@]}"
print_args

targets=("D_C" "D_Fe" "D_Pb" "C" "Fe" "Pb")

# move to main dir
cd ${HAPRAD_CPP}/Utilities

# (related to centroids)
# create `binning.csv` that contains the multi-dimensional edges for each bin
./bin/GetBinning

# loop over targets
for tar in "${targets[@]}"; do

    # (related to centroids)
    # calculates the center-of-mass (or centroids) of each bin
    ./bin/GetCentroids -t${tar} -p${pid}

    # (related to structure functions)
    # fit the `PhiPQ` distributions
    ./bin/FitPhiPQ -t${tar} -p${pid}

    # (final step)
    # get radiative correction factors.
    #
    # The structure-function amplitudes now come from a table file, loaded at
    # run time (see TSemiInclusiveModel.h). This used to source the integrated
    # fit from FitPhiPQ, sed the three constants into TStructFunctionArray.cxx,
    # recompile, run, and sed them back out. Those constants made sigma_Born
    # independent of z and p_t (hurdle H1) and caused the phi = 180 pole (H4).
    # The per-cell producer is PLAN.md Phase 1; until it exists, supply tables
    # by hand as ${SI_TABLE_DIR}/newphihist_<target>.root.
    table="${SI_TABLE_DIR:-${HAPRAD_CPP}/Utilities}/newphihist_${tar}.root"
    if [[ -f "${table}" ]]; then
        ./bin/GetRC -t${tar} -e${BEAM_ENERGY:-5.015} -m"${table}"
    else
        echo "WARNING: no structure-function table ${table}; skipping GetRC for ${tar}."
    fi

    # return to Utilities
    cd ${HAPRAD_CPP}/Utilities
done
