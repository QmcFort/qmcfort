Example 04: Free-Projection AFQMC Calculation with a NOCI Trial Wave Function
=============================================================================

Task
----

Perform an AFQMC calculation using a non-orthognal configuration interaction (NOCI) trial wave function

Prerequisites
-------------

* **Example 01**

  For details about Hartree-Fock and AFQMC calculations

* **Example 03** 

  For details about the frozen-core approximation

Overview
--------

1. HF Calculation with Evaluation of Molecular Integrals and Frozen-Core Approximation
2. NOCI Wave Function Selection
3. Free Projection (fp)-AFQMC Calculation using the NOCI Trial Wave Function

To run the example, either follow step-by-step instructions provided below, or add the following lines 
to the bottom of ``example.sh``:

.. code-block:: bash

    clean_dir
    run_hf
    run_noci
    run_afqmc

and execute the script (or submit it as a batch job).

Step 1: HF Calculation with Evaluation of Molecular Integrals and Frozen-Core Approximation
-------------------------------------------------------------------------------------------

Input
~~~~~

* ``qmcfort.xyz``

* ``basis_set``

* ``qmcfort_in_hf``

.. code-block:: text

    # qmcfort_in_hf

    integral_mode = cholesky
    chol_tol = 1.0E-05
    
    write_files = 2
    
    ispin = 1
    spin = 0
    
    hfock = .true.
    hf_maxiter = 100
    ldiis = .true.
    rdm_tol = 1.0E-06
    hf_tol = 1.0E-10
    
    basis = mo
    frozen_core = .true.

This input performs multiple tasks in a sequence:

1. Evaluation of molecular integrals using Cholesky decomposition,  
2. Hartree–Fock calculation,  
3. Transformation of molecular integrals to the MO basis, and
4. Application of the frozen-core approximation.

For more details on Hartree-Fock calculation and frozen-core approximation, see **Examples 01 and 03**, respectively.

Run
~~~

Execute:

.. code-block:: bash

    source example.sh
    clean_dir
    run_hf

Output
~~~~~~

Hamiltonian and output files are generated.

Searching for ``hf_etot`` in ``qmcfort_out`` or ``qmcfort_log`` should yield:

.. code-block:: text

    hf_energy of the electronic-ionic system:
    ---------------------------------------------------------

    contribution                     Re E
    ------------                     ----
    hf_enuc =                       9.20490701
    hf_e0   =                       0.00000000
    hf_e1   =                    -122.38648498
    hf_eh   =                      47.32427797
    hf_ex   =                      -9.10556571
    hf_e    =                     -84.16777273
    hf_etot =                     -74.96286572

while searching for ``frozen_core_etot`` yields:

.. code-block:: text

    frozen_core_energy of the electronic-ionic system:
    ---------------------------------------------------------

    contribution                     Re E
    ------------                     ----
    frozen_core_enuc =              0.00000000
    frozen_core_e0   =              0.00000000
    frozen_core_e1   =            -65.40875294
    frozen_core_eh   =              9.48896611
    frozen_core_ex   =             -4.74448305
    frozen_core_e    =            -60.66426989
    frozen_core_etot =            -60.66426989

Step 2: NOCI Wave Function Selection
------------------------------------

Input
~~~~~

* ``qmcfort_in_noci``

.. code-block:: text

    # qmcfort_in_noci

    write_files = 0

    brng = lcg48
    brng_seed = 177602913603057

    ispin = 1
    spin = 0

    noci_afqmc = .true.
      ndet_max = 2000
      nepochs = 10
      steps_per_epoch = 100
      ovlp_thresh = 0.6
      energy_thresh_max = 0.0005
      energy_thresh_min = 0.00001
      sigma_eloc = 4.0
      update_afqmc_trial = .true.
    
      #afqmc
        projection = phaseless
        hybrid = 1.0
        afqmc_mix = 1.0

        subtract_mean_field = .true.
        importance_sampling_shift = mix

        nwalkers = 3200
        tau = 0.01
        nblocks = 1
        eqblocks = 1
        steps_per_block = 1
        sample = 1
        samplex = 1
        reorth = 1
        pop = 10

        prop = s2

        exchange_mode = cholesky

        p_capping = .true.
        shift_cut = none


This step performs a short AFQMC simulation while promoting the most important walkers
into a NOCI wave function (``noci_afqmc = .true.``).  
For more details, see the `original publication <https://pubs.acs.org/doi/full/10.1021/acs.jctc.5c00127>`_. 

Run
~~~

Perform NOCI selection:

.. code-block:: bash

    source example.sh
    run_noci

Output
~~~~~~

The otuput files ``qmcfort_out``, ``qmcfort_log``, and ``timing`` are updated.
The NOCI coefficients, and NOCI orbitals are written to: 

* ``ci_coeff_noci``
* ``orbitals_trial_noci``.

Searching for ``Final number of NOCI determinants`` in ``qmcfort_out``, yields otuput similar to:

.. code-block:: text

    Final number of NOCI determinants           11
                det no.         |c_i|^2        weight          H_ii           E_i
                =======         =======        ======          ====           ===
     noci:f:     1              0.9083         0.2652        -74.9629       -75.0077
     noci:f:     2              0.0487         0.1410        -74.3610       -74.0164
     noci:f:     3              0.0055         0.1250        -74.2208       -73.8498
     noci:f:     4              0.0124         0.0534        -73.9160       -73.7706
     noci:f:     5              0.0125         0.0941        -74.2014       -73.6417
     noci:f:     6              0.0020         0.0390        -73.6421       -73.3980
     noci:f:     7              0.0032         0.0038        -73.2826       -73.3255
     noci:f:     8              0.0020         0.0435        -73.8507       -73.2282
     noci:f:     9              0.0024         0.1239        -74.1717       -73.0462
     noci:f:    10              0.0010         0.0717        -74.0220       -72.8544
     noci:f:    11              0.0020         0.0395        -73.9232       -72.6458

The number of determinants :math:`N_d` are listed, together with the NOCI energies in the :math:`E_i` column.
The **first row of the E_i column** corresponds to the **variational NOCI energy**.

Step 3: Free Projection (fp)-AFQMC Calculation with a NOCI Trial Wave Function
------------------------------------------------------------------------------

Input
~~~~~

* ``ci_coeff_noci``

* ``orbitals_trial_noci``

* ``qmcfort_in_afqmc``

.. code-block:: text

    # qmcfort_in_afqmc

    write_files = 0

    brng = lcg48
    brng_seed = 274302913829064

    ispin = 1
    spin = 0

    afqmc = .true.
      projection = free
      hybrid = 1.0
      afqmc_mix = 1.0

      subtract_mean_field = .true.
      importance_sampling_shift = none

      spin_projected_walkers = .true.
    
      nwalkers = 3200
      tau = 0.01
      nblocks = 30
      eqblocks = 10
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

Here, ``projection = free`` activates a free-projection AFQMC simulation.
It is recommended to disable the force bias (``importance_sampling_shift = none``) in fp-AFQMC simulations.

As long as ``ci_coeff_noci`` and ``orbitals_trial_noci`` are present in the working directory, the 
NOCI wave function will be used as a trial wave function.

Run
~~~

Perform fp-AFQMC simulation:

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
    afqmc_enuc =                    9.20490701
    afqmc_e0   =                  -60.66426989
    afqmc_e1   =                  -41.47007958   -0.00097123    7.433053E-04       32.17      4.215715E-03
    afqmc_eh   =                   21.94270648    0.00172250    9.078570E-04       33.48      5.252932E-03
    afqmc_ex   =                   -4.02524267   -0.00075500    1.964268E-04       47.46      1.353269E-03
    afqmc_e    =                  -84.21688566   -0.00000373    1.153650E-05       50.54      8.201316E-05
    afqmc_etot =                  -75.01197865   -0.00000373    1.153650E-05       50.54      8.201316E-05

To double check that the NOCI trial wave function is used, search for ``Trial wave function descriptor``
in the ``qmcfor_out`` file, which should reveal:

.. code-block:: text

    Trial wave function descriptor
    -------------------------------
       trial wave function type                     = WaveTrialNOCI
       size of the object                           =    1.0MB
       number of determinants                       =       11
       Exchange energy mode                         = eri
       block size for Hartree terms                 =      200
       block size for exchange terms                =      200
       inplce energy evaluation                     =      T

or search for ``trial_etot`` (trial energy :math:`<T|H|T>`):

.. code-block:: text

    trial_energy of the electronic-ionic system:
    ---------------------------------------------------------
    contribution                     Re E            Im E
    ------------                     ----            ----
    trial_enuc =                    9.20490701
    trial_e0   =                  -60.66426989
    trial_e1   =                  -41.46488930   -0.00000000
    trial_eh   =                   21.94551067   -0.00000000
    trial_ex   =                   -4.02898779    0.00000000
    trial_e    =                  -84.21263631   -0.00000000
    trial_etot =                  -75.00772930   -0.00000000