import os
import sys

test_modules = ["test_string.f90",
                "test_sparse.f90",
                "test_standalone.f90", 
                "test_statistics.f90",
                "test_file_handle.f90",
                "test_array_file_io.f90",
                "test_method_base.f90",
                "test_qmcfort_pos.f90",
                "test_qmcfort_in.f90",
                "test_mc_descriptor.f90",
                "test_growth_estimator.f90",
                "test_energy_types.f90",
                "test_brng.f90",
                "test_lapack.f90",
                "test_regression.f90",
                "test_expm.f90",
                "test_krylov.f90",
                "test_diis_mixer.f90",
                "test_gauss_base.f90",
                "test_boys_functions.f90",
                "test_mc_murchie_davidson.f90",
                "test_hfproc.f90",
                "test_slater.f90",
                "test_slater_condon_rules.f90",
                "test_wave_trial_md.f90",
                "test_afqmc_proj.f90"]

driver_exe = "fruit_driver"
ext = ".f90"
driver = driver_exe + ext
build_command = "make " + driver_exe

arglen = len(sys.argv)
if (arglen >= 2):
    vers = sys.argv[1]
else:
    vers = "std"
build_command += " VERSION=" + vers

if (arglen >=3):
    nprocs = " -j4"
else:
    nprocs = " -j1"
build_command += nprocs

version = int(sys.version[0])
subversion = int(sys.version[2])

if (version == 2):
    import imp
    try:
        imp.find_module("FRUIT")
        found = True
    except ImportError:
        found = False
elif (version == 3):
    if (subversion >= 9):
        import importlib
        spam_spec = importlib.util.find_spec("FRUIT")
        found = spam_spec is not None
    else:
        import importlib
        spam_spec = importlib.find_loader("FRUIT")
        found = spam_spec is not None
else:
    print("python version =", version, " is not supported")

found = False
if (found):
    from FRUIT import *
    suite = test_suite(test_modules)
    suite.build_run(driver, build_command)
    suite.summary()
else:
    print("FRUITpy module is not found - defualt pure fortran driver will be used")
    os.system(build_command)
    #os.system("mpirun -np 1 ./" + driver_exe)
    os.system("./" + driver_exe)