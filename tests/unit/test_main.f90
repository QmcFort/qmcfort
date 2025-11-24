! SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
! SPDX-License-Identifier: Apache-2.0

program tests

  use fruit

  use mpi, only: comm_world
  use test_afqmc_proj, only: test_afqmc_proj_set
  use test_array_bin_file_io, only: test_array_bin_file_io_set
  use test_array_file_io_factories, only: test_array_file_io_factories_set
  use test_array_numpy_file_io, only: test_array_numpy_file_io_set
  use test_array_txt_file_io, only: test_array_txt_file_io_set
  use test_boys_functions, only: test_boys_functions_set
  use test_diis_mixer, only: test_diis_mixer_set
  use test_energy_types, only: test_energy_types_set
  use test_expm_cheb_mod, only: test_expm_cheb_mod_set
  use test_expm_krylov_mod, only: test_expm_krylov_mod_set
  use test_expm_mod, only: test_expm_mod_set
  use test_expm_taylor_mod, only: test_expm_taylor_mod_set
  use test_file_handle, only: test_file_handle_set
  use test_gauss_base, only: test_gauss_base_set
  use test_hamilton_operations, only: test_hamilton_operations_set
  use test_hfproc, only: test_hfproc_set
  use test_krylov, only: test_krylov_set
  use test_lapack, only: test_lapack_set
  use test_lcg32_basic_rng, only: test_lcg32_basic_rng_set
  use test_lcg48_basic_rng, only: test_lcg48_basic_rng_set
  use test_lcg64_basic_rng, only: test_lcg64_basic_rng_set
  use test_mc_descriptor, only: test_mc_descriptor_set
  use test_mc_murchie_davidson, only: test_mc_murchie_davidson_set
  use test_method_base, only: test_method_base_set
  use test_qmcfort_in, only: test_qmcfort_in_set
  use test_qmcfort_pos, only: test_qmcfort_pos_set
  use test_regression, only: test_regression_set
  use test_slater, only: test_slater_set
  use test_slater_condon_rules, only: test_slater_condon_rules_set
  use test_statistics, only: test_statistics_set
  use test_sparse, only: test_sparse_set
  use test_standalone, only: test_standalone_set
  use test_string, only: test_string_set
  use test_vsl_basic_rng, only: test_vsl_basic_rng_set

  implicit none

  integer :: failed_count

  call init_fruit()
  call comm_world%init()

  !src/afqmc/afqmc_proj.f
  call test_afqmc_proj_set()

  !src/base_io/array_bin_file_io.f
  call test_array_bin_file_io_set()

  !src/base_io/array_file_io_factories.f
  call test_array_file_io_factories_set()

  !src/base_io/array_numpy_file_io.f
  call test_array_numpy_file_io_set()

  !src/base_io/array_txt_file_io.f
  call test_array_txt_file_io_set()

  !src/gto/boys_functions.f
  call test_boys_functions_set()

  !src/base_io/file_handle
  call test_file_handle_set()

  !src/diis/diis_mixer
  call test_diis_mixer_set()

  !src/energy_types.f
  call test_energy_types_set()

  !src/expm/expm_cheb_mod.f
  call test_expm_cheb_mod_set()

  !src/expm/expm_krylov_mod.f
  call test_expm_krylov_mod_set()

  !src/expm/expm_mod.f9
  call test_expm_mod_set()

  !src/expm/expm_taylor_mod.f
  call test_expm_taylor_mod_set()

  !src/gauss_base.f
  call test_gauss_base_set()

  !src/hamilton_operations.f
  call test_hamilton_operations_set()

  !src/hfproc.f
  call test_hfproc_set()

  !src/krylov.f 
  call test_krylov_set()

  !src/lapack.f
  call test_lapack_set()

  !src/random/basic_rng/lcg32_basic_rng.f
  call test_lcg32_basic_rng_set()

  !src/random/basic_rng/lcg48_basic_rng.f
  call test_lcg48_basic_rng_set()

  !src/random/basic_rng/lcg48_basic_rng.f
  call test_lcg64_basic_rng_set()

  !src/mc_descriptor.f
  call test_mc_descriptor_set()

  !src/gto/mc_murchie_davidson.f
  call test_mc_murchie_davidson_set()

  !src/method_base.f
  call test_method_base_set()

  !src/base_io/qmcfort_in
  call test_qmcfort_in_set() 

  !src/qmcfort_pos.f
  call test_qmcfort_pos_set()

  !src/regression.f
  call test_regression_set()

  !src/slater.f
  call test_slater_set()

  !src/slater_condon_rules.f
  call test_slater_condon_rules_set()

  !src/sparse.f
  call test_sparse_set()

  !src/standalone.f
  call test_standalone_set()

  !src/statistics.f
  call test_statistics_set()

  !src/string.f
  call test_string_set()

  !src/random/basic_rng/vsl_basic_rng.f
  call test_vsl_basic_rng_set()


  call get_failed_count(failed_count)
  call fruit_summary
  call fruit_finalize
  if (failed_count > 0) stop 1
end program tests
