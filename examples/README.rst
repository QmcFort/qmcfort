Examples
========

The ``examples/`` directory provides a collection of ready-to-run input files,
reference outputs, and bash scripts illustrating how to use **QmcFort** for
different types of electronic-structure calculations.  
All examples follow a similar directory layout and include a dedicated ``README.rst``
describing the files and the concrete workflow.

Prerequisites
-------------

* The environment variable ``qmcfort_root`` must point to the QmcFort repository;
* The ``qmcfort`` binary must exist in ``$qmcfort_root/bin/`` directory;
* ``mpirun`` must be available on your system.

Running an Example
------------------

To keep your ``examples/`` directory clean, we recommend copying the entire example
folder into your working directory ``my_dir``:

.. code-block:: bash

    cd my_dir
    cp -r $qmcfort_root/examples .

Each example contains an ``example.sh`` script that runs the full workflow.
Alternatively, you may follow the step-by-step instructions provided in the
example-specific ``README.rst`` file.

Parallel Execution
------------------

All examples use at most 16 CPU cores.
For optimal performance, make sure that the product

``mpi_ranks * omp_ranks``

does not exceed the number of available CPU cores on your machine.