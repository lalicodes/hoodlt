"""
:module: FileCleaning
:platform: Unix, Windows
:synopsis: Utility functions for file manipulation and streamlining processes.

.. moduleauthor:: Andres Garcia <andresg@iastate.edu> August 2021
.. history:
..                Andres Garcia <andresg@iastate.edu> August 2021
..                  -Created and tested the functions.
..
"""

# Imports: General.
import copy
import csv
import glob
import numpy as np
import pandas
import os
import re
import shutil
import zipfile

from matplotlib import pyplot
from matplotlib.gridspec import GridSpec

# Imports: File management.
from hoodlt.FileManagement.CreateExpanseScript import CreateExpanseScript as Es


class FileCleaningParameters:
    """ Set the parameters to use with the clean_files function.

        :param self.expanse_scripts: A list of properly instantiated expanse
        scripts.

        :param self.files_per_folder: The number of files per folder for the
        lattice simulations to be distributed. This determines how many folders
        will be created.

        :param self.folder_origin: The name of the folder where the files to be
        cleaned are located.

        :param self.folder_simulations: The base name of the folders where the
        simulation files will be grouped.

        :param self.folder_zip: The name of the zip folder where the simulation
        files that will NOT be moved will be stored.

        :param self.offset_folder_simulations: The offset of the counter to
        generate the numbered folders of the simulation group folders.

        :param self.files_copy_list: The list of extra files that must be copied
        to ALL the simulation folders.
    """

    # --------------------------------------------------------------------------
    # Getters, Setters and Deleters.
    # --------------------------------------------------------------------------

    @property
    def expanse_scripts(self):
        """ Returns the expanse scripts, i.e., the list of expanse scripts to be
            generated.
        """
        return self.__expanse_scripts

    @expanse_scripts.setter
    def expanse_scripts(self, expanse_scripts):
        """ Sets the number of files per folder.
        """

        # Validate that the scripts variable is a list and contains the proper element.
        if not isinstance(expanse_scripts, (list,)):
            raise ValueError(f"The expanse scripts list must be a list.")

        elif len(expanse_scripts) > 0 and not all(map(lambda x: isinstance(x, (Es,)), expanse_scripts)):
            raise ValueError(f"The expanse scripts list must be a list of Expanse scripts.")

        self.__expanse_scripts = expanse_scripts

    @expanse_scripts.deleter
    def expanse_scripts(self):
        """ Deletes this parameter.
        """
        del self.__expanse_scripts

    # --------------------------------------------------------------------------

    @property
    def files_copy_list(self):
        """ Returns the files copy list, i.e., the list of extra files that must
            be copied to the simulation folders.
        """
        return self.__files_copy_list

    @files_copy_list.setter
    def files_copy_list(self, files_copy_list):
        """ Sets the number of files per folder.
        """

        # Validate that the number of files per folder is greater than zero.
        if not isinstance(files_copy_list, (list,)) or not all(map(lambda x: isinstance(x, (str,)), files_copy_list)):
            raise ValueError(f"The files copy list must be a list of strings.")

        self.__files_copy_list = files_copy_list

    @files_copy_list.deleter
    def files_copy_list(self):
        """ Deletes this parameter.
        """
        del self.__files_copy_list

    # --------------------------------------------------------------------------

    @property
    def files_per_folder(self):
        """ Returns the number of files per folder.
        """
        return self.__files_per_folder

    @files_per_folder.setter
    def files_per_folder(self, files_per_folder):
        """ Sets the number of files per folder.
        """

        # Validate that the number of files per folder is greater than zero.
        if files_per_folder < 0:
            raise ValueError(f"The number of files per folder must be greater than or"
                             f" equal to zero. Current files_per_folder: {files_per_folder}.")

        self.__files_per_folder = files_per_folder

    @files_per_folder.deleter
    def files_per_folder(self):
        """ Deletes this parameter.
        """
        del self.__files_per_folder

    # --------------------------------------------------------------------------

    @property
    def folder_origin(self):
        """ Returns the name of the origin folder.
        """
        return self.__folder_origin

    @folder_origin.setter
    def folder_origin(self, folder_origin):
        """ The name of the folder that is to be cleaned.
        """

        # Check that this is a string.
        if not isinstance(folder_origin, (str,)):
            raise ValueError(f"The origin folder must be a string. Current: {type(folder_origin)}")

        self.__folder_origin = folder_origin

    @folder_origin.deleter
    def folder_origin(self):
        """ Deletes this parameter.
        """
        del self.__folder_origin

    # --------------------------------------------------------------------------

    @property
    def folder_simulations(self):
        """ Returns the name of the simulations folder.
        """
        return self.__folder_simulations

    @folder_simulations.setter
    def folder_simulations(self, folder_simulations):
        """ The name of the folder that is to be cleaned.
        """

        # Check that this is a string.
        if not isinstance(folder_simulations, (str,)):
            raise ValueError(f"The simulation group folders must be a string. Current: {type(folder_simulations)}")

        self.__folder_simulations = folder_simulations

    @folder_simulations.deleter
    def folder_simulations(self):
        """ Deletes this parameter.
        """
        del self.__folder_simulations

    # --------------------------------------------------------------------------

    @property
    def folder_zip(self):
        """ Returns the name of the zip folder.
        """
        return self.__folder_zip

    @folder_zip.setter
    def folder_zip(self, folder_zip):
        """ The name of the zip folder where to store the remaining files
            after moving the lattice related files.
        """

        # Check that this is a string.
        if not isinstance(folder_zip, (str,)):
            raise ValueError(f"The zip folder name must be a string. Current: {type(folder_zip)}")

        self.__folder_zip = folder_zip

    @folder_zip.deleter
    def folder_zip(self):
        """ Deletes this parameter.
        """
        del self.__folder_zip

    # --------------------------------------------------------------------------

    @property
    def offset_folder_simulations(self):
        """ Returns the number by which the folder simulation numbering must be
            offset.
        """
        return self.__offset_folder_simulations

    @offset_folder_simulations.setter
    def offset_folder_simulations(self, offset_folder_simulations):
        """ The name of the folder that is to be cleaned.
        """

        # Check that this is a string.
        if not isinstance(offset_folder_simulations, (int,)):
            raise ValueError(f"The offset must be an integer number."
                             f" Current type: {type(offset_folder_simulations)}")

        elif offset_folder_simulations < 0:
            raise ValueError(f"The offset must be an integer number greater than or equal to zero."
                             f" Current value: {offset_folder_simulations}")

        self.__offset_folder_simulations = offset_folder_simulations

    @offset_folder_simulations.deleter
    def offset_folder_simulations(self):
        """ Deletes this parameter.
        """
        del self.__offset_folder_simulations

    # --------------------------------------------------------------------------
    # Constructor(s)
    # --------------------------------------------------------------------------

    def __init__(self):
        """ Initializes the parameters to their default value.
        """

        # ----------------------------------------------------------------------
        # Expanse scripts are initially empty.
        # ----------------------------------------------------------------------

        self.expanse_scripts = []

        # ----------------------------------------------------------------------
        # Number of files per folder.
        # ----------------------------------------------------------------------

        # Set the number of files per folder.
        self.files_per_folder = 4

        # ----------------------------------------------------------------------
        # Folder parameters.
        # ----------------------------------------------------------------------

        # Set the origin folder to be cleaned.
        self.folder_origin = ""

        # Set the base name for the folders where the simulations will be distributed.
        self.folder_simulations = "Simulations"

        # Set the name of the zip folder where the documents that were moved will be stored.
        self.folder_zip = "zipFolder"

        # ----------------------------------------------------------------------
        # Other folder lists.
        # ----------------------------------------------------------------------

        # Set the other files to be copied into the simulation folders.
        self.files_copy_list = [""]

        # ----------------------------------------------------------------------
        # Numbering Offset Parameters.
        # ----------------------------------------------------------------------

        # The numbering offset for the simulations folders.
        self.offset_folder_simulations = 0

    # --------------------------------------------------------------------------
    # Dunder Methods(s)
    # --------------------------------------------------------------------------

    def __repr__(self):
        """ Returns string representation with the parameter values.

            :return self.__str__(): A string with the values of the current
            parameters.
        """
        return self.__str__()

    def __str__(self):
        """ Returns string representation with the parameter values.

            :return strng: A string with the values of the current parameters.
        """

        # Get the list of names of the expanse scripts.
        expanse_script_list = [script.script_name for script in self.expanse_scripts]

        # Message to the user.
        strng = f"FileCleaningParameters values:\n\n"

        # Append the differente parameters with their values.
        strng += f"\tself.expanse_scripts: {expanse_script_list}\n"

        # Append the differente parameters with their values.
        strng += f"\tself.files_copy_list: {self.files_copy_list}\n"

        # Append the differente parameters with their values.
        strng += f"\tself.files_per_folder: {self.files_per_folder}\n"

        # Append the name of the folder where the cleaning will take place.
        strng += f"\tself.folder_origin: {self.folder_origin}\n"

        # Append the name of the base name for the group folders the simulations will be moved.
        strng += f"\tself.folder_simulations: {self.folder_simulations}\n"

        # Append the name of the base name for the zip folder where the remaining files will be moved.
        strng += f"\tself.folder_zip: {self.folder_zip}\n"

        # Append the umbering offset for the simulations folders.
        strng += f"\tself.offset_folder_simulations: {self.offset_folder_simulations}\n"

        return strng


class FileCleaning:
    """ Several useful cleaning utility functions.
    """

    @staticmethod
    def clean_files(parameters):
        """ The parameters required for the files to be cleaned.

            :param parameters: The parameters required to perform a proper
            cleaning. Must be an instance of FileCleaningParameters.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def create_expanse_scripts(folder_names0, script_list0):
            """ Creates the given expanse scripts in the given folder names.

                :param folder_names0: The names of the folders where the expanse
                scripts must be created.

                :param script_list0: The list of Expanse scripts to be
                generated.
            """

            # Create the given scripts in each folder.
            for folder0_1 in folder_names0:
                for script0 in script_list0:
                    script0.generate_script(folder0_1)

        def create_folders(destination_path0, number_of_files0, files_per_folder0, offset0=0):
            """ Creates the folders where the groups of files will be stored,
                according to a determined amount of files to be stored in these
                folders.

                :param destination_path0: The path where the folders are to be
                created.

                :param number_of_files0: The number of files that serve as an
                indication of how many folders must be created.

                :param files_per_folder0: How many files can be stored in a
                folder.

                :param offset0: The offset of the numbering for the folders;
                i.e., the number at which the numbering must start.

                :return folder_names0: The names of the created folders.
            """

            # Auxiliary variables.
            cntr0 = offset0

            # Create a list where to store the file names.
            folder_names0 = []

            # Determine the numbering and how many files must be created.
            for _ in range(offset0, offset0 + number_of_files0, files_per_folder0):
                folder_name_string0 = destination_path0 + f"{cntr0}"

                # Search for the next valid folder name.
                while os.path.isdir(folder_name_string0):
                    folder_name_string0 = destination_path0 + f"{cntr0}"
                    cntr0 += 1

                # Add the folder name to the list.
                folder_names0 += [folder_name_string0]

                # Create the folder.
                os.mkdir(folder_name_string0)

                # Remember to create just one folder if the number of files per
                # folder is zero.
                if files_per_folder0 == 0:
                    break

            return folder_names0

        def remove_files(origin_path0):
            """ Removes the *.log, *.gsd and *.pickle files from the folder,
                leaving it clean.

                :param origin_path0: The path of the folder where the files to
                be deleted are.
            """

            # Get the proper path of the origin file.
            origin_path0_1 = origin_path0.strip()
            origin_path0_1 += "" if origin_path0_1[-1] == os.sep else os.sep

            # Define the search patterns for the files to be deleted.
            search_patterns0 = ("*.log", "*.gsd", "*.pickle")

            # Get the list of files to be deleted.
            files_to_be_deleted0 = []
            for pattern0 in search_patterns0:
                files_to_be_deleted0 += glob.glob(origin_path0_1 + pattern0)

            # Remove the files.
            for files0 in files_to_be_deleted0:
                os.remove(files0)

        def sort_files(files_to_sort0, folder_names0, files_per_folder0):
            """ Sorts the given files in the given folders, given the number of
                files per folder.

                :param files_to_sort0: The files to sort.

                :param folder_names0: The name of the folders where the files
                will be sorted.

                :param files_per_folder0: The number of files per folder.
            """

            # Auxiliary variables.
            cntr0 = 0
            cntr1 = 0
            cntr2 = files_per_folder0

            # Add the os separator if needed.
            folder_names0 = list(map(lambda x: x + "" if x[-1] == os.sep else x + os.sep, folder_names0))

            # Get the raw list of files.
            file_names0 = list(map(lambda x: x.split(os.sep)[-1], files_to_sort0))

            # Number of files.
            file_number0 = len(file_names0)

            # Move all the files to a single folder if it is the case.
            if len(folder_names0) == 1:
                # Get the destination path.
                destination0 = [folder_names0[0] + file_name0 for file_name0 in file_names0]

                # Move the files.
                for i0, file_name0 in enumerate(file_names0):
                    shutil.move(file_name0, destination0[i0])

                return

            # Move the files to the corresponding folders.
            for i in range(0, file_number0, files_per_folder0):
                # Get the origin paths.
                tmp_list0 = files_to_sort0[cntr1: cntr2]

                # Get the destination paths.
                tmp_list1 = file_names0[cntr1: cntr2]
                tmp_list1 = list(map(lambda x:  folder_names0[cntr0] + x, tmp_list1))

                # Move the files accordingly.
                list(map(lambda x, y: shutil.move(x, y), tmp_list0, tmp_list1))

                # Increase the counters
                cntr0 += 1
                cntr1 += files_per_folder0
                cntr2 += files_per_folder0 if cntr1 <= file_number0 else file_number0

        def validate_files(files_names0, lattice_constants0):
            """ Validate that the number of gsd and pickle files are the same,
                as well as the lattice constants.

                :param files_names0: The file names with lattice constant
                indicators.

                :param lattice_constants0: The values of the lattice constants.
            """

            # The number of files must be the same.
            expected_length = len(files_names0[0])

            # All lengths must be the same.
            valid = all(map(lambda x: len(x) == expected_length, files_names0))
            valid = valid and all(map(lambda x: len(x) == expected_length, lattice_constants0))

            if not valid:
                raise ValueError(f"The number of lattice files must be the same for each file type, as"
                                 f" well as the number of lattices. Currently, this is not the case.")

        def validate_parameters():
            """ Validates if the parameter variable is the correct type.
            """

            # Validate that the passed parameters are valid.
            if not isinstance(parameters, (FileCleaningParameters,)):
                raise ValueError(f"The parameters to clean the file must be of"
                                 f"type {FileCleaningParameters}. Current type: {type(parameters)}.")

        def zip_files(origin_path0, zip_file_name0):
            """ Zips ALL the files in the folder.
            """

            # Get the proper folder.
            origin_path0 = origin_path0.strip()
            origin_path0 += "" if origin_path0[-1] == os.sep else os.sep

            # Use glob to look for all the FILES in the folder and filters them.
            files_in_folder0 = list(filter(os.path.isfile, glob.glob(origin_path0 + "*")))

            # Zip the files.
            FileCleaning.zip_files(files_in_folder0, origin_path0, zip_file_name0)

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Validate the parameters.
        validate_parameters()

        # ----------------------------------------------------------------------
        # Zip all the FILES for convenience.
        # ----------------------------------------------------------------------

        # Zip the files.
        zip_files(parameters.folder_origin, parameters.folder_zip)

        # ----------------------------------------------------------------------
        # Obtain the lattice file names.
        # ----------------------------------------------------------------------

        # Get the lattice related files to move.
        function0 = FileCleaning.get_lattice_files_multiple_extensions
        file_names, lattice_constants = function0(parameters.folder_origin, ["gsd", "pickle"])

        # Validate that the number of gsd and pickle files are the same, as well as the lattice constants.
        validate_files(file_names, lattice_constants)

        # Get the specific file names for each extension.
        gsd_files, pickle_files = file_names[0], file_names[1]

        # ----------------------------------------------------------------------
        # Create the folders for the different sets of simulations.
        # ----------------------------------------------------------------------

        # Get the base path.
        folder_path = parameters.folder_origin.strip()
        folder_path += "" if folder_path[-1] == os.sep else os.sep
        folder_path += parameters.folder_simulations

        # Get the offset.
        offset = parameters.offset_folder_simulations

        # Get number of files per folder.
        files_per_folder = parameters.files_per_folder

        # Create the number of required folders.
        folder_names = create_folders(folder_path, len(gsd_files), files_per_folder, offset)

        # ----------------------------------------------------------------------
        # Sort the files to their destination folders.
        # ----------------------------------------------------------------------

        # Sort the files to the given folders.
        sort_files(gsd_files, folder_names, files_per_folder)
        sort_files(pickle_files, folder_names, files_per_folder)

        # Copy the other required scripts to the folders.
        for folder0 in folder_names:
            FileCleaning.copy_files(parameters.folder_origin, folder0, parameters.files_copy_list, new_folder=False)

        # ----------------------------------------------------------------------
        # Create the expanse scripts.
        # ----------------------------------------------------------------------

        # Save the expanse files.
        create_expanse_scripts(folder_names, parameters.expanse_scripts)

        # ----------------------------------------------------------------------
        # Delete the *.gsd, *.log and *.pickle files.
        # ----------------------------------------------------------------------

        # Delete files.
        remove_files(parameters.folder_origin)

    @staticmethod
    def copy_files(origin_path, destination_path, search_patterns=(" ",), new_folder=False):
        """ Copies the requested files to the destination path. If the
            destination path does NOT exist, it tries to create the directory
            using a consecutively numbered pattern.

            :param origin_path: The absolute path of where the files to be moved
            are.

            :param destination_path: The absolute path of where the files are
            to be moved.

            :param search_patterns: The search patterns to be used to search for
            the files using glob.glob. An empty string means that no search will
            be performed.

            :param new_folder: If a new folder must be created in case the
            folder already exists.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def get_files_to_copy(origin_path0, search_patterns0):
            """ Gets the files given the implemented search patterns.

                :param origin_path0: The absolute path of where the files to
                be moved are located.

                :param search_patterns0: The search patterns to be used to
                search for the files using glob.glob.

                :return to_move_names0: The list of the files to be moved.
            """

            # Get the simulation files to move to the temporary folder.
            to_move_names0 = []
            for search_pattern0 in search_patterns0:
                to_move_names0 += glob.glob(origin_path0 + search_pattern0)

            # Make sure the file names are not repeated.
            to_move_names0 = list(np.unique(to_move_names0))

            return to_move_names0

        def normalize_folder_name(folder_name1):
            """ Trims the string and removes operating system path separators at
                the end of the string.

                :param folder_name1: The name of the folder.

                :return folder_name10: The "normalized folder name."
            """

            # Remove leading and trailing white spaces.
            folder_name10 = folder_name1.strip()

            # Remove operating system path separators.
            while folder_name10[-1] == os.sep:
                folder_name10 = folder_name10[:-1]

            return folder_name10.strip()

        def validate_folder_name(destination_path0):
            """ Accordingly, re-numbers the directory where the temporary files
                are to be saved,0 if the passed directory exists.

                :param destination_path0: The path of the directory to created.

                :return folder_name1: The final name the folder will have.
            """

            # Normalize the folder name.s
            destination_path_name1 = normalize_folder_name(destination_path0)

            # Set the folder without numbering, initially.
            folder_name1 = destination_path_name1

            i0 = 0
            while os.path.isdir(folder_name1):
                folder_name1 = destination_path_name1 + str(i0)
                i0 += 1

            # Include the path separator if needed.
            folder_name1 += "" if folder_name1[-1] == os.sep else os.sep

            # Create the folder.
            os.mkdir(folder_name1)

            return folder_name1

        def validate_path(destination_path0):
            """ Validates if the folder path exists for files to be saved.

                :param destination_path0: The path of the directory where the
                files are to be moved.

                :return destination_path0: The validated path where the file is
                to be saved.
            """

            # Remove trailing and leading white spaces.
            destination_path0 = destination_path0.strip()

            # Validate the name.
            if not os.path.isdir(destination_path0):
                raise ValueError("The folder where to save the file does not exist.")

            # Include the separator at the end if needed.
            destination_path0 += "" if destination_path0[-1] == os.sep else os.sep

            return destination_path0

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # No need to go forth if all the search patterns are empty strings.
        if all(map(lambda x: x.strip() == "", list(search_patterns))):
            return

        # Set the proper origin path.
        origin_path = normalize_folder_name(origin_path)
        origin_path += "" if origin_path[-1] == os.sep else os.sep

        # Get the files to move.
        to_move_names = get_files_to_copy(origin_path, list(search_patterns))

        # No point in creating a new folder if there are no files to copy.
        if len(to_move_names) == 0:
            return

        # Get the full path of the destination folder.
        destination_path_name0 = validate_folder_name(destination_path) if new_folder else validate_path(destination_path)

        # Get the origin and destination file names.
        origin_path1 = [origin_path + x.split(os.sep)[-1] for x in to_move_names]
        destination_path1 = [destination_path_name0 + x.split(os.sep)[-1] for x in origin_path1]

        # Copy the files to final directory.
        for i, _ in enumerate(origin_path1):
            shutil.copy(origin_path1[i], destination_path1[i])

    @staticmethod
    def get_lattice_files_single_extension(directory_path, extension, prestring="", poststring=""):
        """ Gets the names of the files that contain a lattice parameter
            specification, sorted in order of increasing lattice constant, along
            with the lattice constants. A prestring and poststring to the file
            search can be included in case a more specific search is required.

            :param directory_path: The path of the directory where the search is
            to be made.

            :param extension: The extension of the requested files; without the
            final dot separator.

            :param prestring: Restrics the search to be more specific by adding
            extra text to the first part of the file name.

            :param poststring: Restrics the search to be more specific by adding
            extra text to the second part of the file name before the extension.

            :return: An unpacked 2-dimensional tuple with the list of organized
            file names and lattice constants.
        """

        # Normalize the directory path.
        directory_path0 = directory_path.strip()
        directory_path0 += "" if directory_path0[-1] == os.sep else os.sep

        # Use glob to obtain all the file names with the requested extensions.
        file_names = [name for name in glob.glob(directory_path0 + prestring + "*" + poststring + "." + extension)]

        # Auxiliary arrays.
        tmp_files0 = []
        lattice_constants0 = []

        # Get the file names with the specific extension.
        file_names0 = [name for name in file_names if name.split(".")[-1] == extension]

        # --------------------------------------------------------------------------
        # Obtain the file names and lattices.
        # --------------------------------------------------------------------------

        # Set the search/filter pattern.
        pattern0 = re.compile(r"_a\d[0-9]+_")
        matches0 = list(map(pattern0.finditer, file_names0))

        # Iterate over the file names.
        for i0, match0 in enumerate(matches0):
            # Iterate over the solutions and store them; one pass is enough.
            for item0 in match0:
                tmp_files0 += [file_names[i0]]
                lattice_constants0 += [float(file_names[i0][item0.start() + 2: item0.end() - 1])]
                break

        # Get the sorted arguments.
        sorted_arguments0 = np.argsort(lattice_constants0)

        # Sort the file names and arguments.
        file_names1 = [tmp_files0[i0] for i0 in sorted_arguments0]
        lattice_constants0 = [lattice_constants0[i0] for i0 in sorted_arguments0]

        return file_names1, lattice_constants0

    @staticmethod
    def get_lattice_files_multiple_extensions(directory_path, extensions_list=(), prestring="", poststring=""):
        """ Gets the names of the files that contain a lattice parameter
            specification, sorted in order of increasing lattice constant, along
            with the lattice constants. A prestring and poststring to the file
            search can be included in case a more specific search is required.

            :param directory_path: The path of the directory where the search is
            to be made.

            :param extensions_list: The extension of the requested files; without
            the final dot separator.

            :param prestring: Restrics the search to be more specific by adding
            extra text to the first part of the file name.

            :param poststring: Restrics the search to be more specific by adding
            extra text to the second part of the file name before the extension.

            :return: An unpacked 2-dimensional tuple with the list of organized
            file names and lattice constants.
        """

        # Variables where the results will be stored.
        file_list = []
        lattice_constant_list = []

        # Get the requested file names and lattice constants for each extension.
        for extension in extensions_list:
            # Get the list for the extension.
            tmp0, tmp1 = FileCleaning.get_lattice_files_single_extension(directory_path, extension, prestring, poststring)

            # Add them to the lists.
            file_list += [tmp0]
            lattice_constant_list += [tmp1]

        return file_list, lattice_constant_list

    @staticmethod
    def log_files_to_data_frame(file_paths, destination_path):
        """ Re-writes the given log files from HOODLT output such that the
            generated files are valid to generate a consistent data frame.

            :param file_paths: The path of the data files.

            :param destination_path: The destination path of the directory where
            the files are to be generated. If the folder does not exist, an
            attempt will be made to create it and write the files there.
        """

        # Do this for each file.
        for _, file in enumerate(file_paths):
            FileCleaning.log_file_to_data_frame(file, destination_path)

    @staticmethod
    def log_file_to_data_frame(file_path, destination_path):
        """ Re-writes the given log file from HOODLT output such that the
            generated files are valid to generate a consistent data frame.

            :param file_path: The path of the data file.

            :param destination_path: The destination path of the directory where
            the files are to be generated. If the folder does not exist, an
            attempt will be made to create it and write the files there.
        """

        # Get the list of data.
        with open(file_path, "r") as fl:
            file_data = list(csv.reader(fl, delimiter="\t"))

        # Get the indexes where the different data frame headers are.
        data_headers = [j for j, name in enumerate(file_data) if name[0].strip()[0] == "#"]

        # Only consider files that have more than one data frame to normalize.
        if len(data_headers) <= 1:
            return

        # Extract the different data frames according to the location of the headers.
        data_frames = []
        for j, index in enumerate(data_headers):
            data_frames += [file_data[index:data_headers[j + 1]] if j < len(data_headers) - 1 else file_data[index:]]

        # Filter out data frames that only have one entry (i.e., only the header).
        data_frames = [data_frame for data_frame in data_frames if len(data_frame) > 1]

        # Remove the pound sign from the data frame headers.
        for j, data in enumerate(data_frames):
            data[0][0] = data[0][0].strip()[1:] if data[0][0][0] == "#" else data[0][0].strip()
            data_frames[j] = data

        # Set the destination path properly.
        destination_path0 = destination_path.strip()
        destination_path0 += "" if destination_path0[-1] == os.sep else os.sep

        # Create the directory if needed.
        if not os.path.isdir(destination_path0):
            os.mkdir(destination_path0)

        # Get the file name and remove the extension.
        destination_path0 += file_path.split(os.sep)[-1]
        destination_path0 = ".".join(destination_path0.split(".")[:-1])

        # Write to numbered files, overwrite.
        for k, data_frame in enumerate(data_frames):
            # Rename the file with a new extension.
            new_file_name = f"{destination_path0}({k}).log"

            # Save the data frame.
            with open(new_file_name, "w") as fl:
                for line in data_frame:
                    # Remove all the leading and trailing spaces in the entries.
                    line = [x.strip() for x in line if not x.strip() == ""]

                    # Write the line to the file.
                    fl.write(("\t".join(line)).strip() + "\n")

    @staticmethod
    def move_files(origin_path, destination_path, search_patterns=(" ",), new_folder=False):
        """ Moves the requested files to the destination path. If the
            destination path does NOT exist, it tries to create the directory
            using a consecutively numbered pattern.

            :param origin_path: The absolute path of where the files to be moved
            are.

            :param destination_path: The absolute path of where the files are
            to be moved.

            :param search_patterns: The search patterns to be used to search for
            the files using glob.glob. An empty string means that no search will
            be performed.

            :param new_folder: If a new folder must be created in case the
            folder already exists.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def get_moving_files(origin_path0, search_patterns0):
            """ Gets the files given the implemented search patterns.

                :param origin_path0: The absolute path of where the files to
                be moved are located.

                :param search_patterns0: The search patterns to be used to
                search for the files using glob.glob.

                :return to_move_names0: The list of the files to be moved.
            """

            # Get the simulation files to move to the temporary folder.
            to_move_names0 = []
            for search_pattern0 in search_patterns0:
                to_move_names0 += glob.glob(origin_path0 + search_pattern0)

            # Make sure the file names are not repeated.
            to_move_names0 = list(np.unique(to_move_names0))

            return to_move_names0

        def normalize_folder_name(folder_name1):
            """ Trims the string and removes operating system path separators at
                the end of the string.

                :param folder_name1: The name of the folder.

                :return folder_name10: The "normalized folder name."
            """

            # Remove leading and trailing white spaces.
            folder_name10 = folder_name1.strip()

            # Remove operating system path separators.
            while folder_name10[-1] == os.sep:
                folder_name10 = folder_name10[:-1]

            return folder_name10.strip()

        def validate_folder_name(destination_path0):
            """ Accordingly, re-numbers the directory where the temporary files
                are to be saved,0 if the passed directory exists.

                :param destination_path0: The path of the directory to created.

                :return folder_name1: The final name the folder will have.
            """

            # Normalize the folder name.s
            destination_path_name1 = normalize_folder_name(destination_path0)

            # Set the folder without numbering, initially.
            folder_name1 = destination_path_name1

            i0 = 0
            while os.path.isdir(folder_name1):
                folder_name1 = destination_path_name1 + str(i0)
                i0 += 1

            # Include the path separator if needed.
            folder_name1 += "" if folder_name1[-1] == os.sep else os.sep

            # Create the directory.
            os.mkdir(folder_name1)

            return folder_name1

        def validate_path(destination_path0):
            """ Validates if the folder path exists for files to be saved.

                :param destination_path0: The path of the directory where the
                files are to be moved.

                :return destination_path0: The validated path where the file is
                to be saved.
            """

            # Remove trailing and leading white spaces.
            destination_path0 = destination_path0.strip()

            # Validate the name.
            if not os.path.isdir(destination_path0):
                raise ValueError("The folder where to save the file does not exist.")

            # Include the separator at the end if needed.
            destination_path0 += "" if destination_path0[-1] == os.sep else os.sep

            return destination_path0

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # No need to go forth if all the search patterns are empty strings.
        if all(map(lambda x: x.strip() == "", list(search_patterns))):
            return

        # Set the proper origin path.
        origin_path = normalize_folder_name(origin_path)
        origin_path += "" if origin_path[-1] == os.sep else os.sep

        # Get the files to move.
        to_move_names = get_moving_files(origin_path, list(search_patterns))

        # No point in creating a new folder if there are no files to move.
        if len(to_move_names) == 0:
            return

        # Get the full path of the destination folder.
        destination_path_name0 = validate_folder_name(destination_path) if new_folder else validate_path(destination_path)

        # Get the origin and destination file names.
        origin_path1 = [origin_path + x.split(os.sep)[-1] for x in to_move_names]
        destination_path1 = [destination_path_name0 + x.split(os.sep)[-1] for x in origin_path1]

        # Move the files to final directory.
        for i, _ in enumerate(origin_path1):
            shutil.move(origin_path1[i], destination_path1[i])

    @staticmethod
    def normalize_path(origin_path):
        """ Formats the path correctly, i.e., removes leading and trailing
            spaces and adds an operating system path-separator at the end, if
            needed.

            :param origin_path: The path to be normalized.

            :return origin_path_0: The normalized path, i.e., the path with leading and trailing
            spaces removed and adds an operating system path-separator at the
            end.
        """

        # Take the trailing and leading spaces off.
        origin_path_0 = origin_path.strip()

        # Add an operating system folder path-separator if needed.
        origin_path_0 += "" if origin_path_0[-1] == os.sep else os.sep

        return origin_path_0

    @staticmethod
    def plot_data_frame(file_path, columns=3, independent_variable=None, save_path=None, show_scatter=True):
        """ Plots the requested data frame. This function is NOT intended to
            customize the plots fully, but to give a quick overview of how the
            dependent variables quantities relate to the independent variable.

            The data frame must be in a tab separated format, i.e., a data frame
            whose column entries are separated by the "tab" character (\t).

            :param file_path: The path of the file the contains the data frame.

            :param columns: The number of columns per row the figure will have.

            :param independent_variable: The name of the independent variable.
            If it is not given, it is assumed to be the first column of the data
            frame. None, by default.

            :param save_path: The directory path where the generated plots will
            be saved. If it is not specified, the image will be saved in the
            directory of origin of the data set.

            :param show_scatter: If the markers of the points must be shown.
        """

        # ----------------------------------------------------------------------
        # Auxiliary Functions.
        # ----------------------------------------------------------------------

        def format_axis(axis0, data_to_plot0, show_scatter0):
            """ Plots the specific quantities, and formats the axis.

                :param axis0: The axis to be formatted.

                :param data_to_plot0: The data to be plotted. Must be a list
                with two entries.

                :param show_scatter0: True if the points should be shown. False,
                otherwise.

                :return axis0: The formatted axis.
            """

            # Plot the data.
            axis0.plot(data_to_plot0[0], data_to_plot0[1])
            if show_scatter0:
                axis0.scatter(data_to_plot0[0], data_to_plot0[1])

            # Customize the axis labels.
            axis0.set_xlabel(independent_variable)
            axis0.set_ylabel(data_to_plot[2])

            # Get the minimum and maximum data points for each axis.
            minx = min(0, min(data_to_plot0[0])) * 1.10
            maxx = max(0, max(data_to_plot0[0])) * 1.10

            # Set the position of the horizontal axes.
            axis0.set_xlim(minx, maxx)
            axis0.spines['left'].set_position(('data', minx))
            axis0.spines['right'].set_position(('data', maxx))

            # Get the minimum and maximum data points for each axis.
            miny = min(0, min(data_to_plot0[1])) * 1.10
            maxy = max(0, max(data_to_plot0[1])) * 1.10

            # Set the position of the vertical axes.
            axis0.set_ylim(miny, maxy)
            axis0.spines['bottom'].set_position(('data', miny))
            axis0.spines['top'].set_position(('data', maxy))

            # Horizontal line from the left-most part to the right most part.
            axis0.plot([minx, maxx], [0, 0], 'black', linewidth=0.25)
            axis0.plot([0, 0], [miny, maxy], 'black', linewidth=0.25)

        def validate_independent_variable(data_columns0, independent_variable0):
            """ Return Determines if the independent variable exists in the data frame
                columns. If it does not exist, it raises an error.

                :param data_columns0: The names of the columns of the
                data frame.

                :param independent_variable0: The name of the independent
                variable.
            """

            # If the independent is the default, no need to validate.
            if independent_variable0 is None:
                return

            # Get the column names.
            if independent_variable0 not in data_columns0:
                raise ValueError(f"The intended data frame column name is not in {data_columns0}."
                                 f"Current independent variable: {independent_variable0}")

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Auxiliary variables.
        cntr = 0

        # ----------------------------------------------------------------------
        # Get the data frame information.
        # ----------------------------------------------------------------------

        # Import the pandas data frame.
        data_frame = pandas.read_csv(file_path, sep="\t")

        # Get the columns from the data frame.
        data_columns = data_frame.columns

        # Validate the independent variable name.
        validate_independent_variable(data_columns, independent_variable)

        # Set the proper indenpendent variable.
        independent_variable = copy.deepcopy(data_columns[0]) if independent_variable is None else independent_variable

        # Get the names of the columns, except for the independent variable.
        data_columns = [column for column in data_columns if not column == independent_variable]

        # --------------------------------------------------------------------------
        # Define the size of the grid.
        # --------------------------------------------------------------------------

        # Calculate the number of rows and columns for the plot.
        rows = len(data_columns) // columns
        rows += int(1) if not len(data_columns) % columns == 0 else int(0)

        # Get the main figure in mid resolution...no need for high resolution.
        fi = pyplot.figure(constrained_layout=True, dpi=600)

        # Set the figure title.
        title = file_path.split(os.sep)[-1]
        fi.suptitle(title, fontsize=10)

        # For the plots to be centered the rows must have an even number of
        # remaining slots.
        remaining_slots = columns - len(data_columns) % columns
        remaining_slots = remaining_slots if not remaining_slots == columns else 0

        # Get the grid specifications according to the dimensions.
        if remaining_slots % 2 == 0:
            gs = GridSpec(rows, columns, figure=fi)
            last_row_offset = remaining_slots // 2

            # Add the graphs.
            for i, data_column in enumerate(data_columns):
                # Keep a count of where we are in terms of rows.
                cntr += 1 if i % columns == 0 and i > 0 else 0

                # Counter to set the axis.
                j = 0 if cntr < rows - 1 else last_row_offset

                axis = fi.add_subplot(gs[cntr, i % columns + j: i % columns + j + 1])

                # Data to plot.
                data_to_plot = (data_frame[independent_variable], data_frame[data_column], data_column)

                # Format the given axis and retrieve its value.
                format_axis(axis, data_to_plot, show_scatter0=show_scatter)

        else:
            gs = GridSpec(rows, 2 * columns, figure=fi)
            last_row_offset = remaining_slots

            # Add the graphs.
            for i, data_column in enumerate(data_columns):
                # Keep a count of where we are in terms of rows.
                cntr += 1 if i % columns == 0 and i > 0 else 0

                # Counter to set the axis.
                j = 0 if cntr < rows - 1 else last_row_offset

                axis = fi.add_subplot(gs[cntr, 2 * (i % columns) + j: 2 * (i % columns) + j + 2])

                # Data to plot.
                data_to_plot = (data_frame[independent_variable], data_frame[data_column], data_column)

                # Format the given axis and retrieve its value.
                format_axis(axis, data_to_plot, show_scatter0=show_scatter)

        # Get the file name without any extension.
        image_name = ".".join((file_path.split(os.sep)[-1]).split(".")[:-1])

        # Set the correct file path in case it is needed.
        file_path = os.sep.join(file_path.strip().split(os.sep)[:-1])
        file_path += "" if file_path == os.sep else os.sep

        # Normalize the save path.
        save_path = save_path.strip() if save_path is not None else file_path
        save_path += "" if save_path == os.sep else os.sep

        # Save the image to the given folder.
        pyplot.savefig(save_path + image_name + ".png")

    @staticmethod
    def zip_files(files_to_zip, destination_path, zip_folder_name="tmp_zip"):
        """ Given a list of files, full path, and the destination folder, it
            zips the given files to in the given zip folder. If the zip folder
            exists, adequate numbering will be added to the folder.
        """

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Check if the destination path is valid.
        destination_path = FileCleaning.normalize_path(destination_path)
        destination_path += "" if destination_path[-1] == os.sep else os.sep

        # Check if the zip FILE is a zip FILE, otherwise increase the numbering.
        j = 0
        base_name = f"{destination_path}{zip_folder_name}.zip"
        while os.path.isfile(base_name):
            base_name = f"{destination_path}{zip_folder_name}{j}.zip"
            j += 1

        # Zip the files.
        with zipfile.ZipFile(base_name, "w") as fl:
            # Remove the directory structure.
            arcnames = list(map(lambda x: x.split(os.sep)[-1], files_to_zip))

            # Zip the files.
            for i, file in enumerate(files_to_zip):
                fl.write(file, arcname=arcnames[i], compress_type=zipfile.ZIP_DEFLATED)


# Main program.
if __name__ == "__main__":
    pass
