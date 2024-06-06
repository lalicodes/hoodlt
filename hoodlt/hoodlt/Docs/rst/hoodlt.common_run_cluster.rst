.. _HOODLTRunCluster:

How to Run a HOODLT Job in a Cluster
====================================

.. _XSEDE: https://portal.xsede.org/

**Logging in to cluster**

This is illustrated for the Expanse system in `XSEDE`_ :

    #. Create your user account on `XSEDE`_ (it takes a few days for XSEDE administrator to set up your account):
    #. If it is your first time logging to the cluster, login to the single sign-on login hub first, and then loggin to the cluster (you'll need to use XSEDE Duo service)
    #. If you've suceeded login on to the single sign-on login hub, from now on you can login to the cluster directly

The following bash script shows how to login to Expanse as an example:

.. code-block:: sh

    ssh <your_username>@login.xsede.org
    gsissh expanse
    exit

    ssh <your_username>@login.expanse.sdsc.edu

Tips: You can generate public/private key pair (e.g. *ssh-keygen*) and append your public key to *~/.ssh/authorized_key* file on the cluster, to avoid having to type your password everytime you login; you can also add login.expanse.sdsc.edu to *~/.ssh/config* file on you local computer to give an alias for the cluster.

**Installing required packages on the cluster**

    #. Install Miniconda (or anaconda)
    #. Create a conda environment (*conda create - -name conda_env_name*)
    #. Install packages: conda install -c conda-forge sphinx git openmpi numpy cmake mpi4py scipy pandas openpyxl gsd
    #. Install hoomd-blue with MPI support (go to `HOOMD Documentation <https://hoomd-blue.readthedocs.io/en/latest/installation.html>`_ for more information.))
    #. Install hoodlt (*git clone https://your_username@bitbucket.org/trvsst/hoodlt.git*, and *python setup.py develop*)

**Submitting jobs**

Following is an example jobscript to be submitted to *gpu-shared* partition on Expanse, for jobscripts on other clusters, go to `XSEDE`_ -> Documentation -> User Guides,

.. code-block:: sh
    :linenos:
    :lineno-start: 2

    #!/bin/bash

    #SBATCH --job-name="name_of_job"
    #SBATCH --output="name_of_job.%j.%N.out"
    #SBATCH --partition=gpu-debug
    #SBATCH --nodes=2
    #SBATCH --gpus=8
    #SBATCH --ntasks-per-node==4
    #SBATCH --export=ALL
    #SBATCH --account=ios116
    #SBATCH --no-requeue
    #SBATCH --time=00:10:00

    source /etc/profile.d/modules.sh
    module purge
    module load gpu
    module load slurm
    module load openmpi
    module load cuda10.2/toolkit

    source /home/your_username/(mini,ana)conda3/etc/profile.d/conda.sh
    conda activate conda_env_name
    source /home/your_username/hoomd-venv/bin/activate

    srun --mpi=pmix -n 8 python my_python_script.py

.. hlist::
    :columns: 1

    * Line 3 specify the job name;
    * Line 4 specify the name of the output file;
    * Line 5 specify on which partition you run, for single GPU jobs, use *gpu-shared* partition, for multiple GPU jobs, use *gpu* partition; for test use *gpu-debug*.
    * Line 6 specify number of nodes, always 1 on *gpu-shared* partition, can be 1, 2, et. al. on *gpu* partition;
    * Line 7 specify number of gpus, always 1 on *gpu-shared* partition, should be 4*nodes on *gpu* partition;
    * Line 8 specify the number of task per node, 4 in this case
    * Line 9 specify which environment variables are loaded;
    * Line 10 is the account name on XSEDE, account name is *ios116*;
    * Line 11 specify whether to requeue the current job if there is a node failure. If there is a node failure, the restart gsd file may not be correctly written into, you'll need to first check restart gsd file and then rerun the job;
    * Line 12 specify the limit for the run time, there is a 48 hour limit on Expanse, 30 minuts in gpu-debug;
    * Line 14 change the shell
    * Line 15-19 load Expanse modules;
    * Line 21-23 activate the appropriate conda environment and other pre-requisites.
    * Last line: launches the job

A simple *my_python_script.py* script (which can be easily adapted to run any simulation with hoodlt) may serve to test that
everything works as it is supposed. This script runs instantaneously.

.. code-block:: python

    from mpi4py import MPI
    import os
    import glob
    import numpy as np
    import hoodlt.Data.Modelconfigurations.Saver as Sr
    import hoodlt.HOOMD.SimParameters as Sp
    from hoodlt.HOOMD.SimulationWithBonds import SimulationWithBonds

    # name of the forcefield
    ff = 'dry-ncs'

    # sorted filenames and a_vals
    names = [os.path.basename(x)[:-4] for x in sorted(glob.glob('cAu*uAngAmuEv.gsd'))]
    a_val = np.array([float(name.split('_a')[1].split('_p')[0]) for name in names])
    arg = np.argsort(a_val)
    a_val = a_val[arg]
    names = [names[arg[i]] for i in range(len(arg))]

    # HoomdSimulation parameters
    steps_wind = 500
    steps_log = 10
    quant_log = ['temperature', 'potential_energy', 'pressure']  # for equilibration process
    steps_dump = int(steps_wind/20)
    sim_params = Sp.SimParameters(steps_wind, steps_log, quant_log, steps_dump)

    # start running in parallel here
    # what follows is a very simplified code meant that should run without errors
    # it should be modified with appropriate HOODLT commands
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    txt = 'file_no_%d'%rank

    fid = open(txt, 'w')

    fid.write('I am in process %d'%rank)

    fid.close()

    print('By %d'%rank)


To submit jobs, use the sbatch command as follows:

.. code-block:: sh

    sbatch jobscriptfile

To check job status, use the following sbatch commands:

.. code-block:: sh

    squeue -p gpu-shared
    squeue -u your_username
