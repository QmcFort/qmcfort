Restricted Open-Shell Hartree-Fock & Restarting AFQMC simulation from a Checkpoint
==================================================================================

Task
----

AFQMC calculation using a restricted open-shell Hartree-Fock (ROHF) trial wave function,
and then continue the simulation from a checkpoint to reduce the statistical error

Prerequisites
-------------

* **Example 01**

  For details about Hartree-Fock and AFQMC calculations

* **Example 03** 

  For details about the frozen-core approximation

Overview
--------

1. Restricted Open-Shell Hartree-Fock (ROHF) Calculation
2. AFQMC Calculation with ROHF Trial Wave Function
3. Continuing an AFQMC Simulation from a Checkpoint

To run the example, either follow step-by-step instructions provided below, or add the following lines 
to the bottom of ``example.sh`` script:

.. code-block:: bash

    clean_dir
    run_hf 
    run_afqmc
    run_afqmc_cont

and execute it (or submit it as a batch job).

Step 1: Restricted Open-Shell Hartree-Fock (ROHF) Calculation
-------------------------------------------------------------

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

    basis = mo

    frozen_core = .true. 

    ispin = 1
    spin = 2

    hfock = .true.
    hf_maxiter = 100
    ldiis = .false.
    rdm_mix = 1.0
    rdm_tol = 1.0E-06
    hf_tol = 1.0E-10

To activate an ROHF calculation, set ``ispin = 1`` together with a non-zero ``spin`` value.
For ROHF it is advisable to disable DIIS (``ldiis = .false.``);
instead, simple density matrix mixing (``rdm_mix``) can be used to stabilize convergence.

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
    hf_enuc =                      28.08768485
    hf_e0   =                       0.00000000
    hf_e1   =                    -261.65459495
    hf_eh   =                     100.25302030
    hf_ex   =                     -16.29454843
    hf_e    =                    -177.69612308
    hf_etot =                    -149.60843823

Step 2: AFQMC Calculation with ROHF Trial Wave Function
-------------------------------------------------------

Input
~~~~~

* ``qmcfort_in_afqmc``

.. code-block:: text

    # qmcfort_in_afqmc

    write_files = 0

    brng = lcg48
    brng_seed = 128702611622883

    ispin = 2
    spin = 2

    afqmc = .true.
      projection = phaseless
      hybrid = 1.0
      afqmc_mix = 1.0

      subtract_mean_field = .true.
      importance_sampling_shift = mix

      spin_projected_walkers = .true.
    
      nwalkers = 800
      tau = 0.05
      nblocks = 40
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

AFQMC does not support the restricted open-shell spin model directly, so ``ispin = 2`` must be used.
For detailed explanations of AFQMC keywords, see **Example 01**.

Run
~~~

Run the AFQMC simulation:

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
    afqmc_enuc =                   28.08768485
    afqmc_e0   =                 -130.50045550
    afqmc_e1   =                  -83.61902288    0.00000000    2.021878E-04       11.27      6.786342E-04
    afqmc_eh   =                   42.06413161    0.00000000    5.592842E-04        5.33      1.290709E-03
    afqmc_ex   =                   -6.00974486    0.00000000    3.084420E-04        1.28      3.494567E-04
    afqmc_e    =                 -178.06509163    0.00000000    2.306260E-04        6.64      5.942112E-04
    afqmc_etot =                 -149.97740678    0.00000000    2.306260E-04        6.64      5.942112E-04

Step 3: Continuing AFQMC Simulation
-----------------------------------

In many situations, one may want to continue the AFQMC simulation to reduce the statistical uncertainty.
This is done by restarting from the checkpoint file ``afqmc_chkpt``.
As long as the ``afqmc_chkpt`` file is present, QmcFort will try to restart/continue AFQMC simulation.

Input
~~~~~

* ``qmcfort_in_afqmc_cont``

.. code-block:: text

    # qmcfort_in_afqmc_cont

    write_files = 0

    brng = lcg48
    brng_seed = 255534416878492

    ispin = 2
    spin = 2

    afqmc = .true.
      projection = phaseless
      hybrid = 1.0
      afqmc_mix = 1.0

      subtract_mean_field = .true.
      importance_sampling_shift = mix

      spin_projected_walkers = .true.
    
      nwalkers = 800
      tau = 0.05
      nblocks = 90
      eqblocks = 10
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

In this particular case, we had 2000 equilibration and 4000 sampling steps in the initial run. 
We now decide to rearrange this into 1000 equilibration steps and 5000 sampling steps, and to 
run 4000 sampling steps more, i.e. 1000 equilibration and 9000 sampling steps in total.
We accomplish this by setting ``eqblocks * steps_per_block = 1000`` and ``nblocks * steps_per_block = 9000``.

Run 
~~~

Continue the AFQMC simulation: 

.. code-block:: bash

    source example.sh
    run_afqmc_cont

Output
~~~~~~

Searching again for ``afqmc_etot`` in ``qmcfort_out`` or ``qmcfort_log``, yields:

.. code-block:: text

    afqmc_energy of the electronic-ionic system:
    ---------------------------------------------------------
    contribution                     Re E            Im E          stddev      corr. length   stddev final 
    ------------                     ----            ----          ------      ------------   ------------ 
    afqmc_enuc =                   28.08768485
    afqmc_e0   =                 -130.50045550
    afqmc_e1   =                  -83.61884490    0.00000000    1.353893E-04       10.91      4.471161E-04
    afqmc_eh   =                   42.06306575    0.00000000    3.872492E-04        5.63      9.185041E-04
    afqmc_ex   =                   -6.00931152    0.00000000    2.179011E-04        1.28      2.466189E-04
    afqmc_e    =                 -178.06554618    0.00000000    1.598493E-04        7.15      4.275356E-04
    afqmc_etot =                 -149.97786133    0.00000000    1.598493E-04        7.15      4.275356E-04

The standard deviation is reduced by a factor of ``2.306260E-04 / 1.598493E-04 = 1.44``, which agrees well  
with the expected value of ``sqrt(9000 / 4000) = 1.5`` given by the **central limit theorem**.