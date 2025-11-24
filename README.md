[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17702833.svg)](https://doi.org/10.5281/zenodo.17702833)

<img src="logo.png" alt="QmcFort Logo" width="300"/>

# Quantum Monte Carlo with Fortran

**QmcFort** is an open-source software package for performing **Auxiliary-Field Quantum Monte Carlo (AFQMC)** calculations.
Written in modern **Fortran**, it is designed for high-accuracy electronic structure simulations
and achieves near-peak performance on modern CPU architectures.

**QmcFort** provides robust implementations of **phaseless** and **free-projection** AFQMC,
as well as the alternative constraints. 
It supports both the **local energy** and **hybrid energy** formalisms, 
and  offers several **propagators** for controlling **time-step errors**.
In addition, **non-orthogonal configuration ineraction (NOCI)** selection algorithm based on the **AFQMC random walk** 
enables systematic improvement of the **AFQMC trial wave function**.

While primarily focused on AFQMC, **QmcFort** also implements:
- Evaluation of Gaussian-type orbital (GTO) integrals using the **McMurchie-Davidson** scheme
- **Hartree-Fock** mean-field theory (**RHF**, **UHF**, **ROHF**)
- Second-order **Møller–Plesset** perturbation theory (**MP2**)

## Download

To download the code, run:
```bash
git clone https://github.com/QmcFort/qmcfort.git
```
This creates a new directory named `qmcfort/` in you current working directory. 

Navigate into it with:
```bash
cd qmcfort
```

## Prerequisites

To build and run QmcFort, the following software must be available on your system:

- **Fortran compiler** (e.g. `ifort`, `gfortran`)
- **BLAS/LAPACK** (e.g. Intel MKL, OpenBLAS) 
- **MPI library** (e.g. Intel MPI, OpenMPI)
- **Python 3** (required for compilation with dependencies and for running regression tests)
- **python3-venv** (optional, used to create an isolated Python environment for QmcFort)

## Building QmcFort

QmcFort uses a top-level `Makefile` located in the project's root directory, together with
several sub-makfiles stored in the `makefiles/` directory.

The entire build and installation process is controlled by this top-level `Makefile`.
Users can build, test, and install QmcFort directly from the root directory using standard `make` commands.

To view all available build targets, run:

```bash
make help
```

### Toolchain 

To build QmcFort, a **toolchain** must be available on your system. 

The toolchain refers to the combination of:

- Fortran compiler
- BLAS/LAPACK library
- MPI libraray

The toolchain is specified as `TOOLCHAIN=<compiler>-<blas_lapack_lib>-<mpi_lib>`.

The following toolchains have been tested so far:
- `TOOLCHAIN=ifort-mkl-intelmpi`   *(default)*
- `TOOLCHAIN=ifort-mkl-openmpi`
- `TOOLCHAIN=gfortran-openblas-openmpi`

You can select a toolchain by editing its value at the beginning of the `Makefile`, or
by passing it explicitly when invoking `make`:

```bash
make TOOLCHAIN=<toolchain> ...
```

Alternatively, set it via an environment variable:  

```bash
export TOOLCHAIN=<toolchain>
```

If nothing is explicitly selected, the default `TOOLCHAIN=ifort-mkl-intelmpi` will be used.

> **__NOTE:__** When using **OpenBLAS**, set the environment variable `OPENBLAS_ROOT`, for example:
```bash
export OPENBLAS_ROOT=/opt/OpenBLAS
```

### Compilation

With a chosen toolchain, you can use parallel compilation with dependencies:

```bash
make DEPS=1 TOOLCHAIN=<toolchain> qmcfort -j4
```

or serial compilation with:

```bash
make TOOLCHAIN=<toolchain> qmcfort
```

> **__NOTE__** `python3` is required because dependency generation relies on Python scripts.
Without these dependencies, parallel compilation may fail.

If completed successfully, the `qmcfort` binary will be created in both `build/<toolchain>/bin/` and `bin/`.

## Testing

QmcFort includes both **unit tests** and **regression tests**, located in the `tests/` directory.

### Unit Tests

To compile and run unit tests, execute:

```bash
make DEPS=1 TOOLCHAIN=<toolchain> ut -j4
```

### Regression Tests

After the `qmcfort` binary has been created, you can run regression tests (`python3` required).
List all available regression tests:

```bash
make rt-list
```

and run regression tests as follows:

```bash
make rt                                          # run all tests
make rt CASES="h2_hf_vdz h2o_afqmc_vdz"          # run two specified tests 
make rt CASES="*afqmc*"                          # run only AFQMC tests
make rt CASES="*h2o*"                            # run only tests on water molecule
```

You may select specific tests by setting the `CASES` variable, which supports simple pattern matching as illustrated in examples above.

The regression tests assume a maximum of 4 processes in total, including MPI ranks and OpenMP threads.
For exact reproducibility, this setup should not be altered.
With a correct setup, the full regression-test suite completes in under five minutes.

## Getting Started & Documentation

A complete user and developer documentation will be provided in one of the upcoming releases.

For now, QmcFort includes an [`examples/`](examples/) directory containing several ready-to-run input files and workflows.
Each example is accompannied by a dedicated `README.rst` file that provides a detailed explantion of the setup,
reuqired input, and expected output.

Since all examples rely on a common set of input and output files, we provide a brief overview of the files used by QmcFort below.

### Overview of the QmcFort Files 

- **`qmcfort_in`**  
    Main input file containing key-value pairs that specify methods, basis sets and algorithms.

- **`qmcfort.xyz`**  
    Molecular geometry in standard `.xyz` file.

- **`basis_set`**  
    List of Gaussian basis functions for each element appearing in `qmcfort.xyz`.

- **`qmcfort_out`**  
    Full, human-readable output file containing energies, diagnostics, and run information.

- **`qmcfort_log`**  
    Compact log file with essential output messages.

- **`timing`**  
    Time-profiling information;

- **`afqmc_chkpt`**  
    Binary checkpoint file used to continue or restart AFQMC simulations.

- **Hamiltonian** Files:
    - **`ham_info`**  
        Text header summarizing the Hamiltonian structure and metadata.

    - **`ham0`**  
        Scalar energy shift (e.g., core energy) - NumPy `.npy` file storing an array of shape $(1)$.

    - **`ham1`**  
        One-electron matrix elements - Numpy `.npy` file storing an array of shape $(N,N,N_s^{(1)})$.

    - **`ham2`**  
        Two-electron integrals:  
            NumPy `.npy` file storing an array of shape $(N,N,N,N,N_s^{(2)})$ if `integral_mode = eri`;  
            NumPy `.npy` file storing an array of shape $(N,N,N_g,N_s^{(1)})$ if `integral_mode = cholesky`.

    - **`overlap`**  
        Overlap integrals - NumPy `.npy` file storing an array of shape $(N,N,N_s^{(1)})$.

    - **`orbitals`**  
        Orbital coefficients - NumPy `.npy` file storing an array of shape $(N,N,N_s)$.

    - **`eigenval`**  
        Orbital energies - NumPy `.npy` file storing an array of shape $(N,N_s)$.

- **NOCI** Files:
    - **`ci_coeff_noci`**  
        CI coefficients of the NOCI wave function - NumPy `.npy` file storing an array of shape $(N_d)$.

    - **`orbitals_trial_noci`**  
        Orbital coefficients for all NOCI determinants - NumPy `.npy` file storing an array of shape $(N,N_{\alpha},N_s,N_d)$.

The variables used to describe array shapes are summarized in the table below.

| Symbol        | Meaning                                              |
|---------------|------------------------------------------------------|
| $N$           | Number of basis functions                            |
| $N_g$         | Number of Cholesky vectors (auxiliary fields)        |
| $N_d$         | Number of Slater determinants                        |
| $N_{\alpha}$  | Number of alpha (spin up) electrons                  |
| $N_s$         | 1 for restricted; 2 for unrestricted orbitals        |
| $N_s^{(1)}$   | Number of spin channels in one-electron integrals    |
| $N_s^{(2)}$   | Number of spin channels in two-electron integrals    |

## QmcFort Tools and Python Virtual Environment

QmcFort uses Python for several tasks, such as **generating build dependencies** and **running regression tests**.
These tasks rely only on the Python standard library and do **not** require any external packages. 

However, additional tools, such as **downloading Gaussian basis sets** or **interfacing QmcFort with PySCF**, depend on external Python libraries.
To avoid interfering with the system Python installation, these tools are meant to run inside a dedicated Python virtual environment.

QmcFort provides `make` rules that automatically create and configure a virtual environment in the `.venv` directory.
Running any of the following commands will create the virtual environment and install the required packages:

```bash
make venv-all                                   # equivalent to: make venv-bse venv-qmcfortpy
make venv-bse                                   # install basis_set_exchange library
make venv-qmcfortpy                             # install Qmcfort-PySCF interface
```
### Working with a Virtual Environment

Activate the virtual environment with:

```bash
source .venv/bin/activate
```

After activation, the commands

```bash
which python
which python3
```
should return paths inside `.venv/bin/`, confirming that the virtual environment is configured properly.

To list installed packages, run:

```bash
pip list
```

To deactivate the virtual environment, run:

```bash
deactivate
```

To remove the entire virtual environment, run:
```bash
make venv-clean
```

## Interface to other Electronic Structure Codes

Although QmcFort includes its own implementation of GTO integral evaluation and a Hartree–Fock mean-field solver---both avoiding the explicit storage of two-electron four-orbital integrals---these components are still experimental.
For production-level calculations, we recommend using a well-established electronic structure package to verify mean-field results or to generate the input files required for AFQMC simulations.


QmcFort currently provides a Python package called **QmcFortPy**, which interfaces QmcFort with [PySCF](https://pyscf.org/index.html).
To install this interface into your virtual environment, run:

```bash
make venv-qmcfortpy
```

as described in the section [QmcFort Tools and and Python Virtual Environment](#QmcFort-Tools-and-Python-Virtual-Environment).

Once installed, the following minimal Python script: 

```python
from pyscf import gto, scf
from qmcfortpy.pyscf import qmcfort_pyscf_itf

mol = gto.M(atom = "H 0 0 0.869; F 0 0 -0.046",
            basis = "cc-pvdz",
            spin = 0)

mf = scf.RHF(mol)
mf.kernel()

qmcfort_pyscf_itf(mol, mf)
```

will generate all necesary input files for running an AFQMC simulation with QmcFort.

Concrete usage instructions and workflow details can be found in [`examples/02-pyscf_qmcfort`](examples/02-pyscf_qmcfort).

## Download / Update Gaussian Basis Sets

The Gaussian basis sets required for the evaluation of molecular integrals are stored under `data/basis_sets`.
To update the database, edit the script `tools/update_basis_sets.py` and run:

```bash
make basis_sets
``` 

This command executes `tools/update_basis_sets.py` script, which downloads the requested basis sets from the 
[Basis Set Exchange](https://www.basissetexchange.org/) using the Python package `basis_set_exchange`.

## How to Contribute

Please read the contribution guidelines in [CONTRIBUTING.md](CONTRIBUTING.md).

## License

QmcFort is distributed under the **Apache License 2.0**.

See the [LICENSE](LICENSE) file for details.

## Citations 

If you use any part of QmcFort software, please cite:

- [Z. Sukurma, et al. *JCTC* **19** (15), 4921-4934 (2023)](https://pubs.acs.org/doi/full/10.1021/acs.jctc.3c00322)
- [Z. Sukurma, et al. *JCTC* **20** (10), 4205-4217 (2024)](https://pubs.acs.org/doi/full/10.1021/acs.jctc.4c00304)

If you use QmcFort for NOCI selection and AFQMC/NOCI calculations, please cite:

- [Z. Sukurma, et al. *JCTC* **21** (9), 4481-4493 (2025)](https://pubs.acs.org/doi/full/10.1021/acs.jctc.5c00127)

## Contact 

QmcFort is primarily developed and maintained by Zoran Sukurma.

If you encounter any issues with installation, compilation, or running the code, you can:

- Open an issue on [GitHub](https://github.com/QmcFort/qmcfort):
- Or send an email to [qmcfort@gmail.com](mailto:qmcfort@gmail.com)

I am happy to help with bug reports, technical questions, and contributions.
