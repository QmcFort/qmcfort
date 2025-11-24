#!/usr/bin/env bash
set -Eeuo pipefail

# read qmcfort_run bash function
source utils.sh

# define qmcfort binary
BIN=$qmcfort_root/bin/qmcfort

# run Hartree-Fock calculation
function run_pyscf() {
  source $qmcfort_root/.venv/bin/activate
  python pyscf_itf.py
  deactivate
}

# run Cholesky decomposition
function run_cholesky() {
  cp qmcfort_in_cholesky qmcfort_in

  mpi_ranks=1
  omp_ranks=8
  blas_ranks=1
  qmcfort_run $mpi_ranks $omp_ranks $blas_ranks

  cp qmcfort_log qmcfort_log_cholesky
  cp qmcfort_out qmcfort_out_cholesky
  cp timing timing_cholesky
}

# run AFQMC calculation
function run_afqmc() {
  cp qmcfort_in_afqmc qmcfort_in
  
  mpi_ranks=16
  omp_ranks=1
  blas_ranks=1
  rm -f afqmc_chkpt
  qmcfort_run $mpi_ranks $omp_ranks $blas_ranks

  cp qmcfort_log qmcfort_log_afqmc
  cp qmcfort_out qmcfort_out_afqmc
  cp timing timing_afqmc
}

# edit here
