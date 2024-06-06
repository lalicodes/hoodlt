"""
:module: CreateExpanseScript.
:platform: Unix, Windows
:synopsis: Generates an Expanse script with the proper parameters.

.. moduleauthor:: Andres Garcia <andresg@iastate.edu> August 2021
.. history:
..                Andres Garcia <andresg@iastate.edu> August 2021
..                  -Created and tested the functions.
..
"""

# Imports: General.
import datetime
import os


class CreateExpanseScript:
    """ Creates an object class that will generate a basic Expanse cluster
        script with the given parameters. It can be later modified.

        :param sefl.account: The name of the account under which the job will
        run.

        :param self.conda_environment_name: The name of the conda environment
        where the Python packages are installed, that include HOODLT.

        :param self.export: Constant that represents what files are to be kept
        when the simulations finish. This is set to ALL.

        :param self.gpus: The number of requested GPUS per node to run the
        simulation.

        :param self.include_date: True if the date must be included in the
        job name and in the output file name. False, otherwise.

        :param self.job_name: The name of the job given to the simulation for
        identification purposes. The date is automatically set and must be kept
        that way.

        :param self.memory: The requested memory that is taken by the job, in
        Gigabytes. The suggested amount of memory is 90 gigabytes, but it can
        be set to a lower number.

        :param self.nodes: The number of nodes requested for the simulation.

        :param self.output: The name of the output file. Must include the
        already added flags.

        :param self.partition: The partition in which the job should be
        processed. The script is currently set for the gpu-shared partition.

        :param self.requeue: If the job should be requeued if the cluster
        fails.

        :param self.script_name_python: The name of the python script to be run
        by the cluster.

        :param self.simulation_time: A 3-tuple of the simulation time given in
        the format (hour, minutes, seconds).

        :param self.tasks_per_node: The number of tasks per node.

        :param self.v_env_name: The virtual environment where HOOMD is
        installed.
    """

    # --------------------------------------------------------------------------
    # Getters, setters and deleters.
    # --------------------------------------------------------------------------

    @property
    def account(self):
        """ Returns the account name parameter.
        """
        return self.__account

    @account.setter
    def account(self, account: str):
        """ Sets the account parameter, that is a string.
        """

        self.__account = f"{account}"

    @account.deleter
    def account(self):
        """ Deletes the account parameter.
        """
        del self.__account

    # --------------------------------------------------------------------------

    @property
    def conda_environment_name(self):
        """ Returns the conda environment name parameter.
        """
        return self.__conda_environment_name

    @conda_environment_name.setter
    def conda_environment_name(self, conda_environment_name: str):
        """ Sets the conda_environment_name parameter, that must be a string.
        """

        self.__conda_environment_name = conda_environment_name

    @conda_environment_name.deleter
    def conda_environment_name(self):
        """ Deletes the conda environment name parameter.
        """
        del self.__conda_environment_name

    # --------------------------------------------------------------------------

    @property
    def export(self):
        """ Returns the export parameter.
        """
        return self.__export

    @export.setter
    def export(self, _):
        """ The export parameter is a constant with the value ALL.
        """

        self.__export = f"ALL"

    @export.deleter
    def export(self):
        """ Deletes the export parameter.
        """
        del self.__export

    # --------------------------------------------------------------------------

    @property
    def gpus(self):
        """ Returns the gpus parameter.
        """
        return self.__gpus

    @gpus.setter
    def gpus(self, gpus: int):
        """ Sets the gpus parameter, that must be an integer.
        """

        # Validate the number of gpus parameters.
        if gpus <= 0:
            raise ValueError(f"The number of gpus must be greater than or equal to 1. Current value: {gpus}")

        self.__gpus = f"{gpus}"

    @gpus.deleter
    def gpus(self):
        """ Deletes the gpus parameter.
        """
        del self.__gpus

    # --------------------------------------------------------------------------

    @property
    def include_date(self):
        """ Returns the include date parameter.
        """
        return self.__include_date

    @include_date.setter
    def include_date(self, include_date: bool):
        """ Sets the include date parameter, that must be an boolean.
        """

        self.__include_date = include_date

    @include_date.deleter
    def include_date(self):
        """ Deletes the include date parameter.
        """
        del self.__include_date

    # --------------------------------------------------------------------------

    @property
    def job_name(self):
        """ Returns the job name parameter.
        """
        return self.__job_name

    @job_name.setter
    def job_name(self, job_name: str):
        """ Sets the job name parameter, that must be a string.
        """

        self.__job_name = f"{job_name}"

    @job_name.deleter
    def job_name(self):
        """ Deletes the job_name parameter.
        """
        del self.__job_name

    # --------------------------------------------------------------------------

    @property
    def memory(self):
        """ Returns the memory parameter.
        """
        return self.__memory

    @memory.setter
    def memory(self, memory: int):
        """ Sets the memory parameter, that must be an integer.
        """

        # Validate the memory parameter.
        if memory <= 0:
            raise ValueError(f"The memory parameter must be greater than or equal to 1. Current value: {memory}")

        self.__memory = f"{memory}G" if memory < 1000 else f"{memory // 1000}T"

    @memory.deleter
    def memory(self):
        """ Deletes the memory parameter.
        """
        del self.__memory

    # --------------------------------------------------------------------------

    @property
    def nodes(self):
        """ Returns the nodes parameter.
        """
        return self.__nodes

    @nodes.setter
    def nodes(self, nodes: int):
        """ Sets the nodes parameter, that must be an integer.
        """

        # Validate the number of nodes parameters.
        if nodes <= 0:
            raise ValueError(f"The number of nodes must be greater than or equal to 1. Current value: {nodes}")

        self.__nodes = f"{nodes}"

    @nodes.deleter
    def nodes(self):
        """ Deletes the nodes parameter.
        """
        del self.__nodes

    # --------------------------------------------------------------------------

    @property
    def output(self):
        """ Returns the output parameter.
        """
        return self.__output

    @output.setter
    def output(self, output: str):
        """ Sets the output parameter, that must be a string.
        """

        self.__output = f'{output}'

    @output.deleter
    def output(self):
        """ Deletes the output parameter.
        """
        del self.__output

    # --------------------------------------------------------------------------

    @property
    def partition(self):
        """ Returns the partition parameter.
        """
        return self.__partition

    @partition.setter
    def partition(self, partition: str):
        """ Sets the partition parameter, that must be a string.
        """

        # The possible values for the partition.
        partition_list = ["compute", "shared", "gpu", "gpu-shared", "large-shared", "debug", "gpu-debug", "preempt", "gpu-preempt"]

        # Validate the partition.
        if partition.lower().strip() not in partition_list:
            raise ValueError(f"The partition must be in the, not case-sensitive, list {partition_list}")

        self.__partition = f'"' + partition.lower().strip() + f'"'

    @partition.deleter
    def partition(self):
        """ Deletes the partition parameter.
        """
        del self.__partition

    # --------------------------------------------------------------------------

    @property
    def requeue(self):
        """ Returns the requeue parameter.
        """
        return self.__requeue

    @requeue.setter
    def requeue(self, requeue: bool):
        """ Sets the requeue parameter, that must be a string.
        """

        self.__requeue = requeue

    @requeue.deleter
    def requeue(self):
        """ Deletes the requeue parameter.
        """
        del self.__requeue

    # --------------------------------------------------------------------------

    @property
    def script_name(self):
        """ Returns the tasks per node parameter.
        """
        return self.__script_name

    @script_name.setter
    def script_name(self, script_name: str):
        """ Sets the script name parameter, that must be a string.
        """
        self.__script_name = script_name

    @script_name.deleter
    def script_name(self):
        """ Deletes the script_name_python parameter.
        """
        del self.__script_name

    # --------------------------------------------------------------------------

    @property
    def script_name_python(self):
        """ Returns the tasks per node parameter.
        """
        return self.__script_name_python

    @script_name_python.setter
    def script_name_python(self, script_name_python: str):
        """ Sets the script name python parameter, that must be a string.
        """
        self.__script_name_python = script_name_python

    @script_name_python.deleter
    def script_name_python(self):
        """ Deletes the script_name_python parameter.
        """
        del self.__script_name_python

    # --------------------------------------------------------------------------

    @property
    def simulation_time(self):
        """ Returns the tasks per node parameter.
        """
        return self.__simulation_time

    @simulation_time.setter
    def simulation_time(self, simulation_time: tuple):
        """ Sets the tasks per node parameter, that must be an integer.
        """

        # Validate the memory parameter.
        if not len(simulation_time) == 3:
            raise ValueError(f"The simulation time must be a 3-tuple of integers. Current value: {simulation_time}")

        elif not isinstance(simulation_time[0], (int,)) or simulation_time[0] < 0 or simulation_time[0] >= 48:
            raise ValueError(f"The simulation_time[0] entry is not valid, is must be an integer greater than zero " +
                             f"and less than 48. Current value: {simulation_time[0]},"
                             f" Type: {type(simulation_time[0])}.")

        elif not isinstance(simulation_time[1], (int,)) or simulation_time[1] < 0 or simulation_time[1] > 59:
            raise ValueError(f"The simulation_time[1] entry is not valid, is must be an integer greater than zero " +
                             f"and less than 59. Current value: {simulation_time[1]},"
                             f" Type: {type(simulation_time[1])}.")

        elif not isinstance(simulation_time[2], (int,)) or simulation_time[2] < 0 or simulation_time[2] > 59:
            raise ValueError(f"The simulation_time[2] entry is not valid, is must be an integer greater than zero " +
                             f"and less than 59. Current value: {simulation_time[2]},"
                             f" Type: {type(simulation_time[2])}.")

        self.__simulation_time = simulation_time

    @simulation_time.deleter
    def simulation_time(self):
        """ Deletes the simulation time parameter.
        """
        del self.__simulation_time

    # --------------------------------------------------------------------------

    @property
    def tasks_per_node(self):
        """ Returns the tasks per node parameter.
        """
        return self.__tasks_per_node

    @tasks_per_node.setter
    def tasks_per_node(self, tasks_per_node: int):
        """ Sets the tasks per node parameter, that must be an integer.
        """

        # Validate the memory parameter.
        if tasks_per_node <= 0:
            raise ValueError(f"The tasks per node parameter must be greater than or equal to 1."
                             f" Current value: {tasks_per_node}")

        self.__tasks_per_node = f"{tasks_per_node}"

    @tasks_per_node.deleter
    def tasks_per_node(self):
        """ Deletes the tasks per node parameter.
        """
        del self.__tasks_per_node

    # --------------------------------------------------------------------------

    @property
    def v_env_name(self):
        """ Returns the virtual environment name parameter.
        """
        return self.__v_env_name

    @v_env_name.setter
    def v_env_name(self, v_env_name: str):
        """ Sets the virtual environment name parameter, that must be a string.
        """

        self.__v_env_name = v_env_name

    @v_env_name.deleter
    def v_env_name(self):
        """ Deletes the virtual environment name parameter.
        """
        del self.__v_env_name

    # --------------------------------------------------------------------------
    # Methods.
    # --------------------------------------------------------------------------

    def generate_script(self, destination_path):
        """ Saves the given script in the provided location with the given name.

            :param destination_path: The absolute path of the folder where the
            script must be saved.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def validate_path():
            # Check that the saving location exists.
            if not os.path.isdir(destination_path):
                raise ValueError(f"The folder where the script is to be saved does not exist."
                                 f" Current name {location}")

            # Check that the destination is not taken.
            if os.path.isdir(full_name) or os.path.isfile(full_name):
                raise ValueError(f"The file already exists, use a different name."
                                 f" Current name {full_name}")

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Get the correct destination path.
        location = destination_path.strip()
        location += "" if location[-1] == os.sep else os.sep

        # Full name.
        full_name = location + self.script_name

        # Validate the path.
        validate_path()

        # Write the file.
        with open(full_name, "w") as fl:
            fl.write(self.__str__())

    # --------------------------------------------------------------------------
    # Constructor(s).
    # --------------------------------------------------------------------------

    def __init__(self):
        """ Constructs an Expanse script object that contains all the needed
            parameters for a basic simulation in the gpu-shared partition. Extra
            parameters can later be added or changed depending on the desired
            type of simulation.
        """

        # ----------------------------------------------------------------------
        # Optional set parameters.
        # ----------------------------------------------------------------------

        # Whether to include the date or not.
        self.include_date = True

        # This is the requeue parameter must be True or False.
        self.requeue = False

        # ----------------------------------------------------------------------
        # General parameters.
        # ----------------------------------------------------------------------

        # This is the name of the account being used.
        self.account = "ios116"

        # This is a constant, for all intents and purposes.
        self.export = "ALL"

        # The number of gpus is initially 1.
        self.gpus = 1

        # Initially the job name is empty.
        self.job_name = ""

        # The requeted memory in gigabytes.
        self.memory = 90

        # The number of nodes is initially 1.
        self.nodes = 1

        # Initially the output name is empty.
        self.output = ""

        # Initially the partition is set to a shared partition.
        self.partition = "gpu-shared"

        # This is the time parameter, that must be a tuple with the time.
        self.simulation_time = (1, 0, 0)

        # The number of tasks per node.
        self.tasks_per_node = 1

        # ----------------------------------------------------------------------
        # Environment names.
        # ----------------------------------------------------------------------

        # Set the name of the conda environment where HOODLT is located.
        self.conda_environment_name = f"/home/<username>/miniconda3/etc/profile.d/conda.sh"

        # Set the name of the virtual environment where HOOMD is located.
        self.v_env_name = f"/home/<username>/hoomd-venv/bin/activate"

        # ----------------------------------------------------------------------
        # Name of the script.
        # ----------------------------------------------------------------------

        # Set the name under which the script will be saved.
        self.script_name = "example"

        # Set the name of the python to be run by the cluster.
        self.script_name_python = "example"

    # --------------------------------------------------------------------------
    # Dunder methods.
    # --------------------------------------------------------------------------

    def __repr__(self):
        """ A quick summary of the value of the settable parameters.
        """

        # Initial message.
        strng = f"Settable parameters for script {self.script_name}:\n"

        # Simulation parameters.
        strng += f"\tself.account: {self.account}\n"
        strng += f"\tself.conda_environment: {self.conda_environment_name}\n"
        strng += f"\tself.export: {self.export}\n"
        strng += f"\tself.gpus: {self.gpus}\n"
        strng += f"\tself.include_date: {self.include_date}\n"
        strng += f"\tself.job_name: {self.job_name}\n"
        strng += f"\tself.memory: {self.memory}\n"
        strng += f"\tself.nodes: {self.nodes}\n"
        strng += f"\tself.output: {self.output}\n"
        strng += f"\tself.partition: {self.partition}\n"
        strng += f"\tself.requeue: {self.requeue}\n"
        strng += f"\tself.script_name: {self.script_name}\n"
        strng += f"\tself.script_name_python: {self.script_name_python}.py\n"
        strng += f"\tself.simulation_time (hours, minutes, seconds): {self.simulation_time}\n"
        strng += f"\tself.source: hoodlt\n"
        strng += f"\tself.tasks_per_node: {self.tasks_per_node}\n"
        strng += f"\tself.v_env_name: {self.v_env_name}\n"

        return strng

    def __str__(self):
        """ Returns a string representation of the current state of the script.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def get_date():
            """ Gets the current date in the format ddmmYYYY.

                :return: The current date in the format ddmmYYYY.
            """
            return datetime.datetime.today().strftime("%d%m") + "20" + datetime.datetime.today().strftime("%y")

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Simulation parameters.
        date = ""
        if self.include_date:
            date = f"-{get_date()}"

        # ----------------------------------------------------------------------
        # Print the script.
        # ----------------------------------------------------------------------

        # Shabang to make it a bash script.
        strng = f"#!bin/bash\n"

        strng += f'#SBATCH --job-name="{self.job_name}{date}"\n'
        strng += f'#SBATCH --output="{self.output}{date}-%j-%N.out"\n'
        strng += f"#SBATCH --partition={self.partition}\n"
        strng += f"#SBATCH --nodes={self.nodes}\n"
        strng += f"#SBATCH --gpus={self.gpus}\n"
        strng += f"#SBATCH --mem={self.memory}\n"
        strng += f"#SBATCH --ntasks-per-node={self.tasks_per_node}\n"
        strng += f"#SBATCH --export={self.export}\n"
        strng += f"#SBATCH --account={self.account}\n"
        if not self.requeue:
            strng += f"#SBATCH --no-requeue\n"

        strnh = f"{self.simulation_time[0]}" if self.simulation_time[0] >= 10 else f"0{str(self.simulation_time[0])}"
        strnm = f"{self.simulation_time[1]}" if self.simulation_time[1] >= 10 else f"0{str(self.simulation_time[1])}"
        strns = f"{self.simulation_time[2]}" if self.simulation_time[2] >= 10 else f"0{str(self.simulation_time[2])}"
        strng += f"#SBATCH --time={strnh}:{strnm}:{strns}\n\n"

        # Location of the cluster modules.
        strng += f"source /etc/profile.d/modules.sh\n\n"

        # Modules to be loaded and imported.
        strng += f"module purge\n"
        strng += f"module load gpu\n"
        strng += f"module load slurm\n"
        strng += f"module load openmpi\n"
        strng += f"module load cuda10.2/toolkit\n\n"

        # Activate the environments.
        strng += f"source {self.conda_environment_name}\n"
        strng += f"source activate hoodlt\n"
        strng += f"source {self.v_env_name}\n\n"

        # The python script to be run.
        strng += f"python {self.script_name}.py"

        return strng
