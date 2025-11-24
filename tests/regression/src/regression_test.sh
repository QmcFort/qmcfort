#!/usr/bin/env bash
set -Eeuo pipefail

# ==============================================================================
# QmcFort regression tests helper functions:
#   - helper functions used in reg_tests/tests/<test>/test.sh
# ==============================================================================

# Detect which MPI implementation mpirun belongs to
function detect_mpi_lib() {
  local v
  v="$(mpirun --version 2>&1 | tr '[:upper:]' '[:lower:]')"
  if [[ "$v" == *"open mpi"* ]]; then
    echo "openmpi"
    return
  elif [[ "$v" == *"intel(r) mpi"* || "$v" == *"intel mpi"* ]]; then
    echo "intelmpi"
    return
  fi
  # Fallbacks if version string is odd
  if command -v ompi_info >/dev/null 2>&1; then
    echo "openmpi"; return
  fi
  if command -v impi_info >/dev/null 2>&1 || [[ -n "$I_MPI_ROOT" ]]; then
    echo "intelmpi"; return
  fi
  echo "unknown"
}

function qmcfort_run() {
    local mpi_ranks=$1
    local omp_ranks=$2

    local mpi_lib="$(detect_mpi_lib)"

    if [ $mpi_lib == "intelmpi" ] 
    then
        mpirun  -np $mpi_ranks \
                -genv OMP_NUM_THREADS $omp_ranks \
                -genv OMP_PLACES cores \
                -genv OMP_PROC_BIND close \
                -genv I_MPI_PIN 1 \
                -genv I_MPI_PIN_RESPECT_CPUSET 0 \
                $BIN
    elif [ $mpi_lib == "openmpi" ]
    then
        mpirun  -np $mpi_ranks \
                --bind-to none \
                -x OMP_NUM_THREADS=$omp_ranks \
                $BIN
    else
        echo "unrecognized mpi library - try with a usual mpirun call"
        export OMP_NUM_THREADS=$omp_ranks
        mpirun -np $mpi_ranks $BIN
    fi
}

function qmcfort_clean() {
    rm -f ham* overlap orbitals eigenval cas_space
    rm -f afqmc_chkpt weights reconf phaser
    rm -f ci_coeff_noci orbitals_trial_noci
    rm -f qmcfort_out* qmcfort_log* timing* qmcfort_proj qmcfort_in values.txt
    rm -f hist_*
}

function qmcfort_clean_noham() {
    rm -f afqmc_chkpt weights reconf phaser
    rm -f qmcfort_out* qmcfort_log* timing* qmcfort_in values.txt
    rm -f hist_*
}