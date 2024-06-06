"""
:module: Squeeze
:platform: Unix, Windows
:synopsis: Generates the configurations of DDAB ligands for a Perovskite face.

.. moduleauthor:: Jonas Hallstrom <jlhallst@asu.edu> August 2021
.. history:
..                Jonas Hallstrom <jlhallst@asu.edu> August 2021
..                  -Created the scripts to generate the DDAB ligand
..                  configuration on a Perovskite cube.
..
"""

# Imports: General imports.
import random
import copy as cp
import numpy as np


class PerovskiteDDABLigandGenerator:
    """ Static class that provides static methods to distribute different
        ligands in a Perovskite cube.
    """

    @staticmethod
    def arrange_ligands(l_size, ligand_number, iterations_max, site_trials_max=10, weights=(1.0, 1.0 / 2.5), extra_detail=False):
        """ Initializes the random distribution of ligands by moving them around
            the face to minimize the coordination number of nearest neighboring
            neighbors, constrained by the given statistical weights.

            :param l_size: The number of cells along the length (or width) of
            the cross-sectional dimension of the cube face.

            :param ligand_number: The total number of ligands to be placed in
            the cube.

            :param iterations_max: The maximum number of attempts to place the
            ligands in the square face.

            :param site_trials_max: The number of attempts to move a ligand to a
            new position; i.e., safeguard to avoid getting stuck in loops; the
            default is 10 trials.

            :param weights: Statistical values to weight the probability of
            placing a ligand at a specific site, given its position within the
            face; tuple with default values (1.0, 1.0 / 2.5).

            :param extra_detail: If more information should be printed to the
            user; False, by default.

            :return ligand_position: The positions of the ligands withing a
            very specific face.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def get_coordination_number(ligand_map0, identifier0, l_size0):
            """ Gets a list with the number of nearest neighbors and next
                nearest neighbors around a site with a given identifier0, based
                on the ligand map. A one represents the presence of a ligand, a
                zero represents an absence of the ligand.

                :param ligand_map0: The map of ligands of the current face. A
                one represents the presence of a ligand, a zero represents an
                absence of the ligand.

                :param identifier0: The identifier of the site.

                :param l_size0: The number of cells along the length (or width)
                of the square face.

                :returns coordination_number0: A list with the number of nearest
                and next nearest neighbor ligands of a cell, based on the ligand
                map.
            """

            # Get the neighbor sites.
            neighbor_sites0 = PerovskiteDDABLigandGenerator.get_neighbors(identifier0, l_size0)

            # Get the number of neighbor and next neighbor ligands to the current ligand.
            neighbor_ligands0 = np.sum(ligand_map0[neighbor_sites0[0]])
            next_neighbor_ligands0 = np.sum(ligand_map0[neighbor_sites0[1]])

            # Get the coordination numbers.
            coordination_number0 = [neighbor_ligands0, next_neighbor_ligands0]
            coordination_number0 = np.array(coordination_number0, dtype=np.double)

            return coordination_number0

        def get_empty_neighbors(ligand_map0, identifier0, l_size0):
            """ Gets a list with the number of nearest neighbors and next
                nearest neighbors around a site with a given identifier0, based
                on the ligand map. A one represents the presence of a ligand, a
                zero represents an absence of the ligand.

                :param ligand_map0: The map of ligands of the current face. A
                one represents the presence of a ligand, a zero represents an
                absence of the ligand.

                :param identifier0: The identifier of the site.

                :param l_size0: The number of cells along the length (or width)
                of the square face.

                :returns coordination_number0: A list with the number of nearest
                and next nearest neighbor ligands of a cell, based on the ligand
                map.
            """

            # Get the neighbor sites.
            neighbor_sites0 = PerovskiteDDABLigandGenerator.get_neighbors(identifier0, l_size0)
            neighbor_sites0 = np.concatenate(neighbor_sites0)

            # Get the considered neighboring sites that have no ligands.
            empty_neighbors0 = [x for x in neighbor_sites0 if ligand_map0[x] == 0]

            return empty_neighbors0

        def get_face_density(ligand_map0, ligand_positions0, weights0, l_size0):
            """ Gets the local density around a site with a given identifier0,
                based on the ligand map. A one represents the presence of a
                ligand, a zero represents an absence of the ligand.

                :param ligand_map0: The map of ligands of the current face. A
                one represents the presence of a ligand, a zero represents an
                absence of the ligand.

                :param ligand_positions0: The list of the sites that contain the
                ligands.

                :param weights0: The statistical weights to calculate the local
                density.

                :param l_size0: The number of cells along the length (or width)
                of the square face.

                :returns The face density, i.e., the sum of the local densities
                where ligands are present.
            """

            # Calculate the local density for each site with a ligand.
            face_density0 = map(lambda x: get_local_density(ligand_map0, x, weights0, l_size0), ligand_positions0)

            # Add the densities.
            face_density0 = np.sum(list(face_density0))

            return face_density0

        def get_local_density(ligand_map0, identifier0, weights0, l_size0):
            """ Gets the local density around a site with a given identifier0,
                based on the ligand map. A one represents the presence of a
                ligand, a zero represents an absence of the ligand.

                :param ligand_map0: The map of ligands of the current face. A
                one represents the presence of a ligand, a zero represents an
                absence of the ligand.

                :param identifier0: The identifier of the site.

                :param weights0: The statistical weights to calculate the local
                density.

                :param l_size0: The number of cells along the length (or width)
                of the square face.

                :returns coordination_number0: A list with the number of nearest
                and next nearest neighbor ligands of a cell, based on the ligand
                map.
            """

            # Get the neighbor sites.
            coordination_number0 = get_coordination_number(ligand_map0, identifier0, l_size0)

            # Get the local density by taking the dot product.
            local_density0 = np.dot(coordination_number0, weights0)

            return local_density0

        def print_message(index, extra_detail0, *args):
            """ Prints the requested message, if extra detail is requested.

                :param index: The index of the message requested to be printed.

                :param extra_detail0: If extra detail is needed.
            """

            def print_matrix(matrix0):
                """ Function that prints a nxm array. If the matrix has further
                    dimensions these are ignored.

                    :param matrix0: The nxm matrix to print.
                """

                # Take each row and format it.
                for i0, submatrix in enumerate(matrix0):
                    strng0 = f"\t" if i0 == 0 else f""

                    # Collect the items on each row (ones or zeros) and print
                    # them.
                    strng1 = f"\t"
                    for j0, element in enumerate(submatrix):
                        strng2 = f"    " if j0 < len(submatrix) - 1 else f""
                        strng1 += f"{int(element)}" + strng2

                    print(strng1)

                # Make sure to leave a blank space.
                print("")

            # No worries if extra details are not needed.
            if not extra_detail0:
                return

            if index == 0:
                print(f"\tLooking for the least dense {args[0]}x{args[0]} face with {args[1]} ligands.\n")

            elif index == 1:
                config = np.reshape(args[0], (args[1], args[1]))
                print(f"\tInitial random configuration is {args[2]}.\n\tConfiguration map:\n")
                print_matrix(config)

            elif index == 2:
                print(f"\t----- Iteration {args[0]} -----\n")

            elif index == 3:
                print(f"\tLigand positions: {args[0]}")
                config = np.reshape(args[1], (args[2], args[2]))
                print(f"\tConfiguration map:\n")
                print_matrix(config)
                print(f"\tFace density: {args[3]}\n")

            elif index == 4:
                print("\tFound the lowest density with the current configuration; breaking the loop.")

            elif index == 5:
                print(f"\tThe same density, {args[0]}, has happened for {args[1]} iterations. Breaking loop.")

            elif index == 6:
                print(f"\tNothing changed this iteration, allowing lateral move next iteration.\n")

            elif index == 7:
                print(f"\tNo ligands have moved for 3 iterations, breaking loop with final density {args[0]}.")

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # Message to the user.
        print_message(0, extra_detail, l_size, ligand_number)

        # ----------------------------------------------------------------------
        # Define the general variables.
        # ----------------------------------------------------------------------

        # Define a matrix that maps the Cs atom ligand occupation: zero is
        # available, one is not available.
        ligand_map = np.zeros(l_size ** 2, dtype=np.int64)

        # Define the maximum number of ligands per face.
        maximum_face_ligands = l_size ** 2

        # Obtain random positions on the face where to place the ligands.
        ligand_position = random.sample(range(maximum_face_ligands), ligand_number)

        # These are the places where the ligands will be placed, initially.
        ligand_map[ligand_position] = int(1)

        # ----------------------------------------------------------------------
        # Further auxiliary variables.
        # ----------------------------------------------------------------------

        # Stores the latest ligand position before a move.
        last_ligand_position = cp.deepcopy(ligand_position)

        # Set the initial density to "infinity".
        last_ligand_density = np.inf

        # If a lateral move is allowed.
        allow_lateral_move = False
        exact_repeats = 1
        density_repeats = 1

        # Message to the user.
        print_message(1, extra_detail, ligand_map, l_size, sorted(cp.deepcopy(ligand_position)))

        # Move the ligands around to minimize the density.
        for i in range(iterations_max):
            # Message to the user.
            print_message(2, extra_detail, i)

            # Shift the ligands arround.
            for current in ligand_position:

                # Obtain the local density from the given site.
                current_density = get_local_density(ligand_map, current, weights, l_size)

                # Obtain the empty neighbors around the specific site.
                empty_neighbors = get_empty_neighbors(ligand_map, current, l_size)

                # Look for lower local densities by trying neighboring empty sites.
                ligand_map[current] = 0
                leader_position = current
                leader_density = current_density

                # Check for prospect moves.
                for prospect in empty_neighbors:
                    # Calculate the density for each possible shift.
                    prospect_density = get_local_density(ligand_map, prospect, weights, l_size)

                    # The prospect becomes the new leader based on whether it is less dense than the previous leader
                    if allow_lateral_move and round(prospect_density, 3) <= leader_density:
                        leader_position = prospect
                        leader_density = prospect_density
                    elif not allow_lateral_move and round(prospect_density, 3) < leader_density:
                        leader_position = prospect
                        leader_density = prospect_density

                # After considering all prospects, the remaining leader will now be the new position of the ligand
                ligand_map[leader_position] = 1
                ligand_position = [j for j, x in enumerate(ligand_map) if x == 1]

            # Calculate the face density.
            face_density = get_face_density(ligand_map, ligand_position, weights, l_size)

            # Print a message to the user.
            print_message(3, extra_detail, ligand_position, ligand_map, l_size, face_density)

            # The face density has achieved a minimum, no need to keep looking.
            if face_density == 0:
                # Message to the user.
                print_message(4, extra_detail)
                break

            # Check the density repeats.
            density_repeats = (density_repeats + 1) if round(face_density, 3) == round(last_ligand_density, 3) else 1
            if density_repeats == site_trials_max:
                # Message to the user.
                print_message(5, extra_detail, face_density, site_trials_max)
                break

            # Determine if lateral moves are allowed.
            allow_lateral_move = sorted(last_ligand_position) == sorted(ligand_position)
            if allow_lateral_move:
                print_message(6, extra_detail)

            # Update the attempt-repeat counter if needed.
            exact_repeats = (exact_repeats + 1) if allow_lateral_move else 1
            if exact_repeats == 3:
                print_message(7, extra_detail, face_density)
                break

            # Set the temporary variables.
            last_ligand_density = cp.deepcopy(face_density)
            last_ligand_position = cp.deepcopy(ligand_position)

            # Shuffles the order of the ligands so that they will move in a different order.
            random.shuffle(ligand_position)

        return ligand_position

    @staticmethod
    def get_faces(sigma, l_size, iterations_max=100, extra_detail=False):
        """ Iterates over the six faces of the cube, to individually arrange the
            ligands, according to the provided density. Since this is a random
            process, the function will only attempt the process a maximum
            number of iterations.

            :param sigma: The requested density of the ligands.

            :param l_size: The number of cells along the length (or width) of
            the cross-sectional dimension of the cube face.

            :param iterations_max: The maximum number of attempts to be made at
            setting up the ligands in the lattice.

            :param extra_detail: If more information should be printed to the
            user; False, by default.

            :return positions: The indexes that represent the positions of the
            ligands on each of the faces.
        """

        # ----------------------------------------------------------------------
        # Parameter definitions.
        # ----------------------------------------------------------------------

        # Cs-Cs nearest neighbor distance (nanometers).
        cs_cs_nn_distance = 0.587

        # Set the area of the square where the ligands are to be placed.
        face_area = (cs_cs_nn_distance * l_size) ** 2

        # Total number of ligands to be placed on ALL the faces; must be an
        # integer number.
        total_ligands = round(face_area * sigma * 6)

        # Minimum number of ligands per faced.
        minimum_ligands = total_ligands // 6

        # Leftover ligands to be distributed among the faces.
        left_over_ligands = total_ligands % 6

        # Initialize the array where the positions are to be stored.
        positions = []

        # ----------------------------------------------------------------------
        # Parameter definitions.
        # ----------------------------------------------------------------------

        # Get the position for each face.
        for i in range(6):
            # Message to the user.
            if extra_detail:
                print(f"\nStarting ligand walker for face {i + 1}:\n")

            # Determine if more ligands must be placed.
            a = 1 if i < left_over_ligands else 0

            # Obtain the positions.
            face_position = PerovskiteDDABLigandGenerator.arrange_ligands(l_size, minimum_ligands + a, iterations_max, extra_detail=extra_detail)
            positions.extend([x * 6 + i for x in face_position])

        # ----------------------------------------------------------------------
        # Print useful information to the user.
        # ----------------------------------------------------------------------

        # Print the extra information to the user.
        if extra_detail:
            print(f"\nTotal number of ligands on the cube is {len(positions)}.")
            print(f"HOODLT-ready positions of the ligands are: {sorted(cp.deepcopy(positions))}.")

        return positions

    @staticmethod
    def get_face_labels(l_size):
        """ Returns a 6 x l_size x l_size matrix that represents the six faces
            of the cube of cross-sectional dimension l_size x l_size.

            :param l_size: The number of cells along the length (or width) of
            the cross-sectional dimension of the cube face.

            :return face_list: The list of faces with the labels of the sites,
            as a numpy 6 x l_size x l_size.
        """

        # Get a cube with the proper dimensions, we will change it later.
        face_list = [0 for _ in range(6)]

        # Loop through the six faces.
        for i in range(6):
            tmp_list = np.array([i + j * 6 for j in range(l_size * l_size)])
            face_list[i] = np.transpose(tmp_list.reshape((l_size, l_size)))

        # Turn the list into a numpy array.
        face_list = np.array(face_list)

        return face_list

    @staticmethod
    def get_neighbors(identifier, l_size):
        """ From the numerical identifier of a cell on a square-face-cell on a
            cube, it returns its neighbor's positions and types.

            :param identifier: The numerical identifier of cell on a
            square-face-cell.

            :param l_size: The length (or width) face-size of a cell.

            :return The site identifier with its neighbor identifiers.
        """

        # ----------------------------------------------------------------------
        # Auxiliary functions.
        # ----------------------------------------------------------------------

        def get_fc_neighbors(matrix, coordinates, l_size0):
            """ Given a set of coordinates, it gets the face cardinal indexes
                from the matrix.

                :param matrix: The matrix from the which the indexes will be
                obtained.

                :param coordinates: The matrix coordinates.

                :param l_size0: The cross-sectional lenght (or width) of the
                square matrix.

                :return fc_list: The list of face cardinals of the given site.
            """

            # Get the additive indexes.
            a_indexes = [1, -1]

            # Get the coordinates.
            crds0 = [[coordinates[0] + i, coordinates[1]] for i in a_indexes]
            crds0 += [[coordinates[0], coordinates[1] + i] for i in a_indexes]

            # Generate the list.
            fc_list = []
            for crd in crds0:
                if PerovskiteDDABLigandGenerator.validate_index(crd, l_size0):
                    fc_list += [matrix[crd[0], crd[1]]]

            # Sort the list.
            fc_list = sorted(fc_list)

            return fc_list

        def get_fd_neighbors(matrix, coordinates, l_size0):
            """ Given a set of coordinates, it gets the face diagonal indexes
                from the matrix.

                :param matrix: The matrix from the which the indexes will be
                obtained.

                :param coordinates: The matrix coordinates.

                :param l_size0: The cross-sectional lenght (or width) of the
                square matrix.

                :return fd_list: The list of face face diagonal indexes of the
                given site.
            """

            # Get the additive indexes.
            a_indexes = [1, -1]

            # Get the coordinates.
            crds1 = []
            for i in a_indexes:
                for j in a_indexes:
                    crds1 += [[coordinates[0] + i, coordinates[1] + j]]

            # Generate the list.
            fd_list = []
            for crd in crds1:
                if PerovskiteDDABLigandGenerator.validate_index(crd, l_size0):
                    fd_list += [matrix[crd[0], crd[1]]]

            # Sort the list.
            fd_list = sorted(fd_list)

            return fd_list

        def generate_labels(l_size0):
            """ Generates an l_size0 x l_size0 matrix with the corresponing
                labels.

                :param l_size0: The cross-sectional lenght (or width) of the
                square matrix.

                :returns square_matrix: The square matrix with the corresponding
                labels of the sites.
            """

            # Obtain the linear array and reshape to the proper size.
            square_matrix = np.array([i for i in range(l_size0 * l_size0)])
            square_matrix = np.transpose(square_matrix.reshape((l_size0, l_size0)))

            return square_matrix

        def validate_identifer(identifier0, l_size0):
            """ Raises an exception if the size of the lattice is not correct.

                :param identifier0: The numerical identifier of cell on a
                square-face-cell.

                :param l_size0: The cross-sectional lenght (or width) of the
                square matrix.
            """
            # Possible identifier types.
            int_types = (int, np.int32, np.int64)

            # Raise an error if an invalid identifier is passed.
            if not 0 <= identifier0 < l_size0 * l_size0 or not isinstance(identifier0, int_types):
                raise ValueError(f"A wrong identifier was passed. It must be in {int_types} between 0 and "
                                 + f"{l_size0 ** 2 - 1}. Identifier is currently: {identifier0}, "
                                 + f"identifier type is {type(identifier0)}.")

        # ----------------------------------------------------------------------
        # Implementation.
        # ----------------------------------------------------------------------

        # First validate the identifier.
        validate_identifer(identifier, l_size)

        # Get the square matrix.
        mat = generate_labels(l_size)

        # Get the coordinate of the system.
        base_crd = [identifier % l_size, identifier // l_size]

        # Remember to sort the lists before.
        fc_neighbors = get_fc_neighbors(mat, base_crd, l_size)
        fd_neighbors = get_fd_neighbors(mat, base_crd, l_size)

        return [fc_neighbors, fd_neighbors]

    @staticmethod
    def validate_index(coordinate, l_size):
        """ Determines if the given matrix cell indexes are valid, given the
            cell indexes and the dimensions of the matrix.

            :param coordinate: The coordinate index to be determined if it
            is valid.

            :param l_size: The size of the matrix.

            :return valid: If the ith coordinate index is in the range
            0 <= index[i] < l_size.
        """

        # Check ALL indexes fall within the proper range.
        valid = True
        for coord in coordinate:
            valid = valid and 0 <= coord < l_size

        return valid
