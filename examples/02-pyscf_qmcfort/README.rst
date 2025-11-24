Example 02: PySCF-QmcFort Interface 
===================================

Task
----

Use PySCF to generate input files for an AFQMC calculation on the HF molecule using a restricted Hartree-Fock
(RHF) trial wave function 

Prerequisites
-------------

* **QmcFortPy**

  Installable using ``make venv-qmcfortpy`` or ``make venv-all``.

* **Example 01**

  For a detailed explanations of the AFQMC procedure, see **Example 01** first.

Overview
--------

1. PySCF Mean-Field Calculation with QmcFortPy Export
2. Cholesky Decomposition of Electron Repulsion Integrals (ERIs)
3. ph-AFQMC Calculation with RHF Trial Wave Function 

To run the example, either follow step-by-step instructions provided below, or add the following lines 
to the bottom of ``example.sh`` script:

.. code-block:: bash

    clean_dir
    run_pyscf
    run_cholesky
    run_afqmc

and execute it (or submit it as a batch job).

Step 1: PySCF Mean-Field Calculation with QmcFortPy Export
----------------------------------------------------------

Input
~~~~~

* ``qmcfort.xyz``

  .. include:: qmcfort.xyz

* ``pyscf_itf.py``

  .. include:: pyscf_itf.py

.. code-block:: python

    # pyscf_itf.py

    from pyscf import gto, scf
    from qmcfortpy.pyscf import qmcfort_pyscf_itf

    mol = gto.M(atom = "qmcfort.xyz",
                basis = "cc-pvdz",
                spin = 0)

    mf = scf.RHF(mol)
    mf.kernel()

    qmcfort_pyscf_itf(mol, mf, basis="ao_orth")

Run
~~~

Run PySCF calculation:

.. code-block:: bash

    source example.sh
    clean_dir
    run_pyscf

Output
~~~~~~

A successful run should print something like:

.. code-block:: text

    converged SCF energy = -100.019481490578
    System Info
    ===========
    n                     =  19
    nbtot                 =  19
    ispin                 =  1
    ispin1                =  1
    ispin2                =  1
    nel                   =  5 5
    basis                 =  ao_orth
    ==========================================

QmcFortPy exports molecular integrals in the orthogonalized AO basis (``basis = ao_orth``), or in MO basis (``basis = mo``). 
These are written into the Hamiltonian files: ``ham_info``, ``ham0``, ``ham1``, ``ham2``, ``orbitals``, ``overlap``, ``eigenval``.

Note that PySCF dumps out the full list of electron-repulsion integrals, whcih becomes
memory intensive for large systems or basis sets.
Since AFQMC requires the low-rank factorized ERIs, we perform Cholesky decomposition in the next step.

Step 2: Cholesky Decomposition of Electron Repulsion Integrals (ERIs)
---------------------------------------------------------------------

Input
~~~~~

* ``qmcfort_in_cholesky``

.. code-block:: text

    # qmcfort_in_cholesky

    write_files = 2

    integral_mode = cholesky
    compress_h2 = .true.
    chol_tol = 1.0E-05

Run
~~~

Perform Cholesky decomposition:

.. code-block:: bash

    source example.sh
    run_cholesky

Output
~~~~~~

After the QmcFort run, ``qmcfort_out``, ``qmcfort_log`` and ``timing`` files are created.
The ``ham2`` file is overwritten to store the Choelsky vectors instead of ERIs.

Searching for ``number of cholesky vectors`` in ``qmcfort_out`` or ``qmcfort_log`` should yield
a number around 115 (the exact value may vary slightly due to the batching).

Step 3: ph-AFQMC Calculation with RHF Trial Wave Function 
---------------------------------------------------------

Input
~~~~~

* ``qmcfort_in_afqmc``

.. code-block:: text

    # qmcfort_in_afqmc

    write_files = 0

    brng = lcg48
    brng_seed = 274302913829064

    ispin = 1
    spin = 0

    afqmc = .true.
      projection = phaseless
      hybrid = 1.0
      afqmc_mix = 1.0

      subtract_mean_field = .true.
      importance_sampling_shift = mix

      spin_projected_walkers = .true.
  
      nwalkers = 6400
      tau = 0.01
      nblocks = 80
      eqblocks = 20
      steps_per_block = 10
      sample = 1
      samplex = 1
      reorth = 1
      pop = 10

      prop = s2
      expm = taylor
      expm_order = 6

      exchange_mode = eri 

      p_capping = .true.
      cut_shift = remove

For a detailed explanation of the AFQMC keywords, see **Example 01**.

Run
~~~

To start the AFQMC calculation, run:

.. code-block:: bash

    source example.sh
    run_afqmc

Output
~~~~~~

After the AFQMC calculation, ``qmcfort_out``, ``qmcfort_log`` and ``timing`` files are upated.

Searching for ``afqmc_etot`` in ``qmcfort_out`` or ``qmcfort_log``, should show:

.. code-block:: text

    afqmc_energy of the electronic-ionic system:
    ---------------------------------------------------------
    contribution                     Re E            Im E          stddev      corr. length   stddev final 
    ------------                     ----            ----          ------      ------------   ------------ 
    afqmc_enuc =                    5.20502212
    afqmc_e0   =                    0.00000000
    afqmc_e1   =                 -150.64234215    0.00000000    2.229496E-04       77.50      1.962701E-03
    afqmc_eh   =                   55.54113678    0.00000000    2.875758E-04       92.90      2.771824E-03
    afqmc_ex   =                  -10.33558290    0.00000000    4.938958E-05       91.05      4.712632E-04
    afqmc_e    =                 -105.43678827    0.00000000    1.059410E-04       51.06      7.570297E-04
    afqmc_etot =                 -100.23176615    0.00000000    1.059410E-04       51.06      7.570297E-04

Performance
~~~~~~~~~~~

Searching for ``average FLOPS`` in the ``timing`` file should yield a number around 20 GFLOPS.