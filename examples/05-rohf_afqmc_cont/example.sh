#!/usr/bin/env bash
set -Eeuo pipefail

# read qmcfort_run bash function
source utils.sh

# define qmcfort binary
BIN=$qmcfort_root/bin/qmcfort

# run Hartree-Fock calculation
function run_hf() {
  cp qmcfort_in_hf qmcfort_in
  
  mpi_ranks=1
  omp_ranks=8
  blas_ranks=1
  qmcfort_run $mpi_ranks $omp_ranks $blas_ranks

  cp qmcfort_log qmcfort_log_hf
  cp qmcfort_out qmcfort_out_hf
  cp timing timing_hf
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

# run AFQMC calculation
function run_afqmc_cont() {
  cp qmcfort_in_afqmc_cont qmcfort_in
  
  mpi_ranks=16
  omp_ranks=1
  blas_ranks=1
  qmcfort_run $mpi_ranks $omp_ranks $blas_ranks

  cp qmcfort_log qmcfort_log_afqmc_cont
  cp qmcfort_out qmcfort_out_afqmc_cont
  cp timing timing_afqmc_cont
}

# edit here
