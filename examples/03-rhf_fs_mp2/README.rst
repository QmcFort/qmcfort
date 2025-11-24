Example 03: MP2 with the Frozen-Core Approximation
==================================================

Task
----

Perform MP2 calculation on the benzene molecule within the frozen-core approximation 

Overview
--------

1. Direct Calculation of the Cholesky Vectors
2. Hartree-Fock Calculation
3. Frozen-Core Approximation
4. MP2 Calculation

To run the example, either follow step-by-step instructions provided below, or add the following lines 
to the bottom of ``example.sh``:

.. code-block:: bash

    clean_dir
    run_cholesky
    run_hf
    run_fs
    run_mp2

and execute the script (or submit it as a batch job).

Step 1: Direct Calculation of the Cholesky Vectors
--------------------------------------------------

Input
~~~~~

* ``qmcfort.xyz``
* ``basis_set``
* ``qmcfort_in_cholesky``

.. code-block:: text

    # qmcfort_in_cholesky

    write_files = 2
    integral_mode = cholesky
    chol_tol = 1.0E-06

It is important to set ``integral_mode = cholesky`` to avoid full (and expensive) evaluation of the electron-repulsion integrals (ERIs).

Run
~~~

Calculate molecular integrals:

.. code-block:: bash 

    source example.sh
    clean_dir
    run_cholesky

Output
~~~~~~

After the QmcFort run, ``qmcfort_out``, ``qmcfort_log`` and ``timing`` files are created. 
Further, Hamiltonian files ``ham0``, ``ham1``, ``ham2``, ``overlap`` are also generated.

Searching for ``number of cholesky vectors`` in ``qmcfort_out`` or ``qmcfort_log`` should yield
a value around 900 (the exact number may vary slightly due to the batching).

Step 2: Hartree-Fock Calculation
--------------------------------

Input
~~~~~

* ``qmcfort_in_hf``

.. code-block:: text

    # qmcfort_in_hf

    hfock = .true.
    hf_maxiter = 100
    ldiis = .true.
    rdm_tol = 1.0E-06
    hf_tol = 1.0E-10

    write_files = 1

The keyword ``write_files = 1`` ensures that only the ``orbitals`` file is overwritten during
the HF calculation, since the molecular integrals remain unchanged.

Run
~~~

Start the HF calculation:

.. code-block:: bash

    source example.sh
    run_hf

Output
~~~~~~

After the run, ``qmcfort_out``, ``qmcfort_log``, ``timing``, and ``obitals`` files are updated.

Searching for ``hf_etot`` in ``qmcfort_out`` or ``qmcfort_log`` file, should reveal:

.. code-block:: text

    hf_energy of the electronic-ionic system:
    ---------------------------------------------------------

    contribution                     Re E      
    ------------                     ----      
    hf_enuc =                     203.15352437
    hf_e0   =                       0.00000000
    hf_e1   =                    -712.71474646
    hf_eh   =                     312.07774535
    hf_ex   =                     -33.23834107
    hf_e    =                    -433.87534218
    hf_etot =                    -230.72181781

Step 3: Frozen-Core Approximation
---------------------------------

Input
~~~~~

* ``qmcfort_in_fs`` 

.. code-block:: text

    # qmcfort_in_fs

    basis = mo
    frozen_core = .true.
    nfrozen = 6

    write_files = 2

Here, ``frozen_core = .true.`` activates the frozen-core approximation, which requires the molecular integrals 
to be in the MO basis (``basis = mo``).
The number of frozen states can be specified explicitly via ``nfrozen`` keyword,
or determined automatically by the code.
The number of core electrons per element follows the same convention as the
`ORCA <https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/frozencore.html>`_ default frozen-core scheme.

.. note::

    The number of frozen electrons equals **two times** the number of frozen states ``nfrozen``.

Run
~~~

Perform the frozen-core approximation:

.. code-block:: bash

    source example.sh
    run_fs

Output
~~~~~~

After the run, ``qmcfort_out``, ``qmcfort_log``, ``timing``, and Hamiltonian files are updated.
The following items in the ``ham_info`` file should change from

.. code-block:: text

    n             =    114
    nbtot         =    114
    nel           =     21    21

to

.. code-block:: text

    n             =    108
    nbtot         =    108
    nel           =     15    15

Searching for ``frozen_core_etot`` shows the energy of the frozen core :math:`<F|H|F>`:

.. code-block:: text

    frozen_core_energy of the electronic-ionic system:
    ---------------------------------------------------------

    contribution                     Re E      
    ------------                     ----      
    frozen_core_enuc =              0.00000000
    frozen_core_e0   =              0.00000000
    frozen_core_e1   =           -331.96539605
    frozen_core_eh   =             58.67505741
    frozen_core_ex   =            -21.02996246
    frozen_core_e    =           -294.32030110
    frozen_core_etot =           -294.32030110

Step 4: MP2 Calculation
-----------------------

Input
~~~~~

* ``qmcfort_in_mp2``

.. code-block:: text

    # qmcfort_in_mp2

    write_files = 0
    mp2 = .true.


Run
~~~

Perform the MP2 calculation:

.. code-block:: bash

    source example.sh
    run_mp2

Output
~~~~~~

Output files ``qmcfort_out``, ``qmcfort_log`` and ``timing`` are updated again.

Searching for ``mp2_etot`` in ``qmcfort_out`` or ``qmcfort_log`` should reveal:

.. code-block:: text

    mp2_energy of the electronic-ionic system:
    ---------------------------------------------------------

    contribution                     Re E      
    ------------                     ----      
    mp2_enuc =                    203.15352437
    mp2_e0   =                   -294.32030110
    mp2_e1   =                   -258.94331781
    mp2_eh   =                    129.59602008
    mp2_ex   =                    -10.99167597
    mp2_e    =                   -434.65927480
    mp2_etot =                   -231.50575043