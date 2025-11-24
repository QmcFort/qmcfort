#!/usr/bin/env bash
set -Eeuo pipefail

# read qmcfort_run bash function
source utils.sh

# define qmcfort binary
BIN=$qmcfort_root/bin/qmcfort

# run Cholesky decomposition
function run_cholesky() {
  cp qmcfort_in_cholesky qmcfort_in

  mpi_ranks=1
  omp_ranks=16
  blas_ranks=1
  qmcfort_run $mpi_ranks $omp_ranks $blas_ranks

  cp qmcfort_log qmcfort_log_cholesky
  cp qmcfort_out qmcfort_out_cholesky
  cp timing timing_cholesky
}

# run Hartree-Fock calculation
function run_hf() {
  cp qmcfort_in_hf qmcfort_in

  mpi_ranks=1
  omp_ranks=1
  blas_ranks=1
  qmcfort_run $mpi_ranks $omp_ranks $blas_ranks

  cp qmcfort_log qmcfort_log_hf
  cp qmcfort_out qmcfort_out_hf
  cp timing timing_hf
}

# run Frozen-core approximation
function run_fs() {
  cp qmcfort_in_fs qmcfort_in

  mpi_ranks=1
  omp_ranks=1
  blas_ranks=1
  qmcfort_run $mpi_ranks $omp_ranks $blas_ranks

  cp qmcfort_log qmcfort_log_fs
  cp qmcfort_out qmcfort_out_fs
  cp timing timing_fs
}

# run MP2 calculation
function run_mp2() {
  cp qmcfort_in_mp2 qmcfort_in
  
  mpi_ranks=1
  omp_ranks=1
  blas_ranks=1
  qmcfort_run $mpi_ranks $omp_ranks $blas_ranks

  cp qmcfort_log qmcfort_log_mp2
  cp qmcfort_out qmcfort_out_mp2
  cp timing timing_mp2
}

# edit here