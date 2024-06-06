"""
:module: Working directory.
:platform: Unix, Windows
:synopsis: Context manager to set a temporary working directory.

.. moduleauthor:: Andres Garcia <andresg@iastate.edu> September 2021
.. history:
..                Andres Garcia <andresg@iastate.edu> September 2021
..                  -Created and tested the functions.
..
"""

# Imports: General.
import copy as cp
import os


class SetupWorkingDirectory:
    """ Context manager to set the temporary working directory to the given
        path.

        It is useful when you need to set the working directory in a different
        place than the default one.

        :param self.directory_path: The path of the directory to be set as the
        working directory.

        :param self.current_working_directory: The current working directory,
        i.e., the directory before the new one is set.
    """

    def __init__(self, directory_path):
        """ Constructs a context manager, i.e., creates the required variables
            to change the working directory.

            :param directory_path: The path of the directory to be set as the
            working directory.
        """
        self.directory_path = directory_path
        self.current_working_directory = cp.deepcopy(os.getcwd())

    def __enter__(self):
        """ Set the working directory to the path and return the string
            that represents the path to the current working directory.

            :return self.directory_path: The path that is set as the current
            working directory.
        """
        os.chdir(self.directory_path)
        return self.directory_path

    def __exit__(self, exc_type, exc_val, exc_tb):
        """ Set the working directory to the original one. Captures and
            returns exceptions, if there are any.

            :param exc_type: Holds the exception type, if any.

            :param exc_val: Holds the value of the exception, if any.

            :param exc_tb: Holds the traceback of the exception, if any.
        """
        os.chdir(self.current_working_directory)