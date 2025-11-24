Example 01: Unrestricted Hartree-Fock & Auxiliary-Field Quantum Monte Carlo
===========================================================================

Task
----

AFQMC calculation using the unrestricted Hartree-Fock (UHF) trial wave function

Overview
--------

.. 1. :ref:`UHF Calculation with Evaluation of Molecular Integrals <uhf_step>`
.. 2. :ref:`ph-AFQMC Calculation with UHF Trial Wave Function <afqmc_step>`

1. UHF Calculation with Evaluation of Molecular Integrals
2. ph-AFQMC Calculation with UHF Trial Wave Function

To run the example, either follow step-by-step instructions provided below, or add the following lines 
to the bottom of ``example.sh`` script:

.. code-block:: bash

    clean_dir
    run_hf 
    run_afqmc

and execute it (or submit it as a batch job).

.. _uhf_step:

Step 1: UHF calculation with Evaluation of Molecular Integrals
--------------------------------------------------------------

Input
~~~~~

* ``qmcfort.xyz``

  .. include:: qmcfort.xyz

* ``basis_set``

  .. include:: basis_set

* ``qmcfort_in_hf``

..  .. include:: qmcfort_in_hf

.. code-block:: text

    # qmcfort_in_hf

    integral_mode = cholesky
    chol_tol = 1.0E-05
    basis = ao_orth

    write_files = 2

    ispin = 2
    spin = 3

    hfock = .true.
    hf_maxiter = 100
    ldiis = .true.
    rdm_tol = 1.0E-06
    hf_tol = 1.0E-10

Here, ``integral_mode = cholesky`` and ``chol_tol = 1.0E-05`` instruct QmcFort to use a low-rank factorization of the ERIs.
Since AFQMC theory assumes an orthogonal basis functions, the keyword ``basis = ao_orth`` is used to orthogonalize AO integrals. 
Alternatively, one may use ``basis = mo`` to transform integrals into the MO basis, which is orthogonal per construction. 
The keyword ``ispin = 2`` requires unrestricted spin model, with ``spin`` being :math:`N_{\alpha} - N_{\beta}`. 

Run
~~~

Run the Hartree-Fock calculation:

.. code-block:: bash 

    source example.sh
    clean_dir
    run_hf

Output
~~~~~~

This will generate Hamiltonian files ``ham0``, ``ham1``, ``ham2``, as well as ``orbitals``, ``overlap``, and ``eigenval`` files. 
It also generates the output files ``qmcfort_out``, ``qmcfort_log``, and ``timing``.

Searching for ``hf_etot`` in ``qmcfort_out`` or ``qmcfort_log`` file, should reveal:

.. code-block:: text

    hf_energy of the electronic-ionic system:
    ---------------------------------------------------------

    contribution                     Re E
    ------------                     ----
    hf_enuc =                       0.00000000
    hf_e0   =                       0.00000000
    hf_e1   =                     -73.96117561
    hf_eh   =                      26.17999484
    hf_ex   =                      -6.60994988
    hf_e    =                     -54.39113065
    hf_etot =                     -54.39113065

Here, the energy components have following meaning:

* ``enuc`` - the classical nuclear repulsion energy between point charges,
* ``e0``   - the core energry, 
* ``e1``   - the one-body electronic energy (kinetic energy + electron-ion potential),
* ``eh``   - the classical Hartree energy of the electron density,
* ``ex``   - the exchange energy.
* ``e``    - the total electronic energy, defined as ``e = e0 + e1 + eh + ex``, and
* ``etot`` - the total molecular energy, ``etot = enuc + e``.

.. _afqmc_step:

Step 2: ph-AFQMC Calculation with UHF trial Wave Function
---------------------------------------------------------

Input
~~~~~

* ``qmcfort_in_afqmc``

..  .. include:: qmcfort_in_afqmc

.. code-block:: text

    # qmcfort_in_afqmc
    
    write_files = 0
    
    brng = lcg48
    brng_seed = 128702611622883
    
    ispin = 2
    spin = 3
    
    afqmc = .true.
      projection = phaseless
      hybrid = 1.0
      afqmc_mix = 1.0
    
      subtract_mean_field = .true.
      importance_sampling_shift = mix
    
      spin_projected_walkers = .true.
    
      nwalkers = 800
      tau = 0.01
      nblocks = 80
      eqblocks = 20
      steps_per_block = 100
      sample = 1
      samplex = 1
      reorth = 1
      pop = 10
    
      prop = s2
      expm_mode = taylor
      expm_order = 6
    
      exchange_mode = eri 
    
      p_capping = .true.
      cut_shift = remove

In this step we perform a phaseless AFQMC calculation (``projection = phaseless``) in a hybrid energy mode (``hybrid = 1.0``).
We use mean-field subtraction (``mean_field_subtraction = .true.``) and a walker-dependent force bias (``importance_sampling_shift = mix``).
The total number of walkers is set via ``nwalkers = 800`` and the time step is ``tau = 0.01``.

The total simulation length is organized into blocks:
``steps_per_block = 100`` specifies the block size,. 
``eqblocks = 20`` defines the number of equilibration blocks (not used for measurements), and
``nblocks = 80`` sets the number of sampling blocks.

Additional keywords control the frequency of specific operations:
``reorth`` (walker reorthogonalization),
``pop`` (walker population control),
``sample`` (energy measurement), and
``samplex`` (exchange energy subsampling relative to ``sample``).

Finally, a pseudorandom number generator (PRNG) is an important component of the AFQMC algorithm.
Here we use ``brng = lcg48`` together with a fixed random seed (``brng_seed``) for reproducibility.
For production calculations, it is recommended to use either
``brng = vsl_brng_mt19937`` (if compiled with Intel MKL),  or
``brng = intrinsic`` (Fortran compiler intrinsic PRNG).
To randomize the simulation, set ``brng_seed = random``. 
To repeat the calculation with exactly the same random number sequence, search for ``brng_seed`` in the 
``qmcfort_out`` file and copy its value into ``qmcfort_in``.

Run
~~~

To perform the AFQMC calculation, run:

.. code-block:: bash

    source example.sh
    run_afqmc

Output
~~~~~~

After the AFQMC calculation, the files ``qmcfort_out``, ``qmcfort_log`` and ``timing`` are upated.

Searching for ``afqmc_etot`` in ``qmcfort_out`` or ``qmcfort_log``, should show:

.. code-block:: text

    afqmc_energy of the electronic-ionic system:
    ---------------------------------------------------------
    contribution                     Re E            Im E          stddev      corr. length   stddev final
    ------------                     ----            ----          ------      ------------   ------------
    afqmc_enuc =                    0.00000000
    afqmc_e0   =                    0.00000000
    afqmc_e1   =                  -73.94889939    0.00000000    9.965658E-05       79.28      8.873189E-04
    afqmc_eh   =                   25.99869290    0.00000000    1.530477E-04       50.11      1.083347E-03
    afqmc_ex   =                   -6.52916896    0.00000000    8.955283E-05        1.36      1.046050E-04
    afqmc_e    =                  -54.47937544    0.00000000    3.528039E-05       46.63      2.409159E-04
    afqmc_etot =                  -54.47937544    0.00000000    3.528039E-05       46.63      2.409159E-04

Here, ``Re E`` and ``Im E`` denote the real and imaginary parts of the energy vlaue.
The ``stddev`` column is the naive estimator of the standard deviation of the mean (SDM),
which ignores the autocorrelation in the data.
The ``corr. length`` reports the correlation length estimated using block averaging,
and ``stddev_final``  is the corresponding corrected SDM that accounts for this correlation.
Thus, the AFQMC energy estimate is intpreted as ``Re E`` + i * ``Im E``  +- ``stddev_final``.

Performance
~~~~~~~~~~~

The ``timing`` file gives profiling information.
For compute-intensive tasks, QmcFort also reports FLOPS performance.
Searching for ``average FLOPS`` shows the average number of floating-point operations per second (FLOPS) per MPI process.
For such a small system size, the expected value is around 10 GFLOPS.