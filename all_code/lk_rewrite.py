############
############ ALTHOUGH I GIVE YOU THE 'BARE BONES' OF THIS PROGRAM WITH THE NAME
############ 'skeleton.py', YOU CAN RENAME IT TO ANYTHING YOU LIKE. HOWEVER, FOR
############ THE PURPOSES OF THE EXPLANATION IN THESE COMMENTS, I ASSUME THAT
############ THIS PROGRAM IS STILL CALLED 'skeleton.py'.
############
############ IF YOU WISH TO IMPORT STANDARD MODULES, YOU CAN ADD THEM AFTER THOSE BELOW.
############ NOTE THAT YOU ARE NOT ALLOWED TO IMPORT ANY NON-STANDARD MODULES!
############

import os
import sys
import time
import random

############
############ NOW PLEASE SCROLL DOWN UNTIL THE NEXT BLOCK OF CAPITALIZED COMMENTS.
############
############ DO NOT TOUCH OR ALTER THE CODE IN BETWEEN! YOU HAVE BEEN WARNED!
############

def read_file_into_string(input_file, ord_range):
    the_file = open(input_file, 'r')
    current_char = the_file.read(1)
    file_string = ""
    length = len(ord_range)
    while current_char != "":
        i = 0
        while i < length:
            if ord(current_char) >= ord_range[i][0] and ord(current_char) <= ord_range[i][1]:
                file_string = file_string + current_char
                i = length
            else:
                i = i + 1
        current_char = the_file.read(1)
    the_file.close()
    return file_string

def remove_all_spaces(the_string):
    length = len(the_string)
    new_string = ""
    for i in range(length):
        if the_string[i] != " ":
            new_string = new_string + the_string[i]
    return new_string

def integerize(the_string):
    length = len(the_string)
    stripped_string = "0"
    for i in range(0, length):
        if ord(the_string[i]) >= 48 and ord(the_string[i]) <= 57:
            stripped_string = stripped_string + the_string[i]
    resulting_int = int(stripped_string)
    return resulting_int

def convert_to_list_of_int(the_string):
    list_of_integers = []
    location = 0
    finished = False
    while finished == False:
        found_comma = the_string.find(',', location)
        if found_comma == -1:
            finished = True
        else:
            list_of_integers.append(integerize(the_string[location:found_comma]))
            location = found_comma + 1
            if the_string[location:location + 5] == "NOTE=":
                finished = True
    return list_of_integers

def build_distance_matrix(num_cities, distances, city_format):
    dist_matrix = []
    i = 0
    if city_format == "full":
        for j in range(num_cities):
            row = []
            for k in range(0, num_cities):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    elif city_format == "upper_tri":
        for j in range(0, num_cities):
            row = []
            for k in range(j):
                row.append(0)
            for k in range(num_cities - j):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    else:
        for j in range(0, num_cities):
            row = []
            for k in range(j + 1):
                row.append(0)
            for k in range(0, num_cities - (j + 1)):
                row.append(distances[i])
                i = i + 1
            dist_matrix.append(row)
    if city_format == "upper_tri" or city_format == "strict_upper_tri":
        for i in range(0, num_cities):
            for j in range(0, num_cities):
                if i > j:
                    dist_matrix[i][j] = dist_matrix[j][i]
    return dist_matrix

def read_in_algorithm_codes_and_tariffs(alg_codes_file):
    flag = "good"
    code_dictionary = {}
    tariff_dictionary = {}
    if not os.path.exists(alg_codes_file):
        flag = "not_exist"
        return code_dictionary, tariff_dictionary, flag
    ord_range = [[32, 126]]
    file_string = read_file_into_string(alg_codes_file, ord_range)
    location = 0
    EOF = False
    list_of_items = []
    while EOF == False:
        found_comma = file_string.find(",", location)
        if found_comma == -1:
            EOF = True
            sandwich = file_string[location:]
        else:
            sandwich = file_string[location:found_comma]
            location = found_comma + 1
        list_of_items.append(sandwich)
    third_length = int(len(list_of_items)/3)
    for i in range(third_length):
        code_dictionary[list_of_items[3 * i]] = list_of_items[3 * i + 1]
        tariff_dictionary[list_of_items[3 * i]] = int(list_of_items[3 * i + 2])
    return code_dictionary, tariff_dictionary, flag

############
############ THE RESERVED VARIABLE 'input_file' IS THE CITY FILE UNDER CONSIDERATION.
############
############ IT CAN BE SUPPLIED BY SETTING THE VARIABLE BELOW OR VIA A COMMAND-LINE
############ EXECUTION OF THE FORM 'python skeleton.py city_file.txt'. WHEN SUPPLYING
############ THE CITY FILE VIA A COMMAND-LINE EXECUTION, ANY ASSIGNMENT OF THE VARIABLE
############ 'input_file' IN THE LINE BELOW iS SUPPRESSED.
############
############ IT IS ASSUMED THAT THIS PROGRAM 'skeleton.py' SITS IN A FOLDER THE NAME OF
############ WHICH IS YOUR USER-NAME, E.G., 'abcd12', WHICH IN TURN SITS IN ANOTHER
############ FOLDER. IN THIS OTHER FOLDER IS THE FOLDER 'city-files' AND NO MATTER HOW
############ THE NAME OF THE CITY FILE IS SUPPLIED TO THIS PROGRAM, IT IS ASSUMED THAT
############ THE CITY FILE IS IN THE FOLDER 'city-files'.
############

input_file = "AISearchfile012.txt"

############
############ PLEASE SCROLL DOWN UNTIL THE NEXT BLOCK OF CAPITALIZED COMMENTS.
############
############ DO NOT TOUCH OR ALTER THE CODE IN BETWEEN! YOU HAVE BEEN WARNED!
############

if len(sys.argv) > 1:
    input_file = sys.argv[1]

the_particular_city_file_folder = "city-files"

if os.path.isfile("../" + the_particular_city_file_folder + "/" + input_file):
    ord_range = [[32, 126]]
    file_string = read_file_into_string("../" + the_particular_city_file_folder + "/" + input_file, ord_range)
    file_string = remove_all_spaces(file_string)
    print("I have found and read the input file " + input_file + ":")
else:
    print("*** error: The city file " + input_file + " does not exist in the folder '" + the_particular_city_file_folder + "'.")
    sys.exit()

location = file_string.find("SIZE=")
if location == -1:
    print("*** error: The city file " + input_file + " is incorrectly formatted.")
    sys.exit()

comma = file_string.find(",", location)
if comma == -1:
    print("*** error: The city file " + input_file + " is incorrectly formatted.")
    sys.exit()

num_cities_as_string = file_string[location + 5:comma]
num_cities = integerize(num_cities_as_string)
print("   the number of cities is stored in 'num_cities' and is " + str(num_cities))

comma = comma + 1
stripped_file_string = file_string[comma:]
distances = convert_to_list_of_int(stripped_file_string)

counted_distances = len(distances)
if counted_distances == num_cities * num_cities:
    city_format = "full"
elif counted_distances == (num_cities * (num_cities + 1))/2:
    city_format = "upper_tri"
elif counted_distances == (num_cities * (num_cities - 1))/2:
    city_format = "strict_upper_tri"
else:
    print("*** error: The city file " + input_file + " is incorrectly formatted.")
    sys.exit()

dist_matrix = build_distance_matrix(num_cities, distances, city_format)
print("   the distance matrix 'dist_matrix' has been built.")

############
############ YOU NOW HAVE THE NUMBER OF CITIES STORED IN THE INTEGER VARIABLE 'num_cities'
############ AND THE TWO_DIMENSIONAL MATRIX 'dist_matrix' HOLDS THE INTEGER CITY-TO-CITY
############ DISTANCES SO THAT 'dist_matrix[i][j]' IS THE DISTANCE FROM CITY 'i' TO CITY 'j'.
############ BOTH 'num_cities' AND 'dist_matrix' ARE RESERVED VARIABLES AND SHOULD FEED
############ INTO YOUR IMPLEMENTATIONS.
############

############
############ THERE NOW FOLLOWS CODE THAT READS THE ALGORITHM CODES AND TARIFFS FROM
############ THE TEXT-FILE 'alg_codes_and_tariffs.txt' INTO THE RESERVED DICTIONARIES
############ 'code_dictionary' AND 'tariff_dictionary'. DO NOT AMEND THIS CODE!
############ THE TEXT FILE 'alg_codes_and_tariffs.txt' SHOULD BE IN THE SAME FOLDER AS
############ THE FOLDER 'city-files' AND THE FOLDER WHOSE NAME IS YOUR USER-NAME, E.G., 'abcd12'.
############

code_dictionary, tariff_dictionary, flag = read_in_algorithm_codes_and_tariffs("../alg_codes_and_tariffs.txt")

if flag != "good":
    print("*** error: The text file 'alg_codes_and_tariffs.txt' does not exist.")
    sys.exit()

print("The codes and tariffs have been read from 'alg_codes_and_tariffs.txt':")

############
############ YOU NOW NEED TO SUPPLY SOME PARAMETERS.
############
############ THE RESERVED STRING VARIABLE 'my_user_name' SHOULD BE SET AT YOUR USER-NAME, E.G., "abcd12"
############

my_user_name = "chpf93"

############
############ YOU CAN SUPPLY, IF YOU WANT, YOUR FULL NAME. THIS IS NOT USED AT ALL BUT SERVES AS
############ AN EXTRA CHECK THAT THIS FILE BELONGS TO YOU. IF YOU DO NOT WANT TO SUPPLY YOUR
############ NAME THEN EITHER SET THE STRING VARIABLES 'my_first_name' AND 'my_last_name' AT
############ SOMETHING LIKE "Mickey" AND "Mouse" OR AS THE EMPTY STRING (AS THEY ARE NOW;
############ BUT PLEASE ENSURE THAT THE RESERVED VARIABLES 'my_first_name' AND 'my_last_name'
############ ARE SET AT SOMETHING).
############

my_first_name = "Finlay"
my_last_name = "Boyle"

############
############ YOU NEED TO SUPPLY THE ALGORITHM CODE IN THE RESERVED STRING VARIABLE 'algorithm_code'
############ FOR THE ALGORITHM YOU ARE IMPLEMENTING. IT NEEDS TO BE A LEGAL CODE FROM THE TEXT-FILE
############ 'alg_codes_and_tariffs.txt' (READ THIS FILE TO SEE THE CODES).
############

algorithm_code = "LK"

############
############ DO NOT TOUCH OR ALTER THE CODE BELOW! YOU HAVE BEEN WARNED!
############

if not algorithm_code in code_dictionary:
    print("*** error: the agorithm code " + algorithm_code + " is illegal")
    sys.exit()
print("   your algorithm code is legal and is " + algorithm_code + " -" + code_dictionary[algorithm_code] + ".")

############
############ YOU CAN ADD A NOTE THAT WILL BE ADDED AT THE END OF THE RESULTING TOUR FILE IF YOU LIKE,
############ E.G., "in my basic greedy search, I broke ties by always visiting the first
############ city found" BY USING THE RESERVED STRING VARIABLE 'added_note' OR LEAVE IT EMPTY
############ IF YOU WISH. THIS HAS NO EFFECT ON MARKS BUT HELPS YOU TO REMEMBER THINGS ABOUT
############ YOUR TOUR THAT YOU MIGHT BE INTERESTED IN LATER.
############

added_note = ""

############
############ NOW YOUR CODE SHOULD BEGIN.
############

from queue import Queue

"""
Calculates the length of a tour
"""
def tour_length_calc(tour):
    dist = dist_matrix[tour[0]][tour[len(tour) - 1]]

    for i in range(0, len(tour) - 1):
        dist += dist_matrix[tour[i]][tour[i + 1]]

    return dist

"""
Converts a tour to an edge list
"""
def get_edge_list_from_tour(tour, list=True):
    # Since we are dealing with the symmetric case for LK we sort the contents of the tuples
    return [tuple(sorted((tour[i], tour[(i + 1) % len(tour)]))) for i in range(len(tour))]

"""
Returns the two adjacent nodes to a node
Index 0 has prev
Index 1 has next
"""
def get_adjacent_nodes(tour, node):
    position = tour.index(node)
    return [tour[position - 1], tour[(position + 1) % len(tour)]]

"""
Returns the two adjacent edges to a node
Index 0 has (prev, node)
Index 1 has (node, next)
"""
def get_adjacent_edges(tour, node):
    position = tour.index(node)
    return [tuple(sorted((tour[position - 1], tour[position]))), tuple(sorted((tour[position], tour[(position + 1) % len(tour)])))]

"""
Returns the 'max_count' best neighbours of 'node' in the 'tour'
(where best is the lowest weight edge)
Can skip nodes if they are in an iterable in excluded_nodes
A list of tuples. tuple[0] is the neighbour, tuple[1] is the weight
We also don't want the pair to exist in the tour already
"""
def get_best_neighbours(tour, node, max_count, excluded_nodes = None):
    # Don't use the default parameter excluded_nodes = [] since the list
    # never resets itself
    if excluded_nodes == None:
        excluded_nodes = set()
    else:
        excluded_nodes = set(excluded_nodes)

    # Don't want to match to itself
    if node not in excluded_nodes:
        excluded_nodes.add(node)

    weights = []

    # Index of the weight is also the city itself
    for neighbour, weight in enumerate(dist_matrix[node]):
        if neighbour in excluded_nodes:
            continue

        weights.append((neighbour, weight))

    # Sort based on the weight with the smallest at index 0
    return sorted(weights, key=lambda k: k[1], reverse=False)[:max_count]

"""
Returns the weight of a single edge
"""
def get_edge_weight(edge):
    return dist_matrix[edge[0]][edge[1]]

"""
When picking y_i's we require that:
1. y_i is not in X
2. G_i - y_i_weight is positive
3. If y_i is selected then there exists an x_(i + 1) that can be broken
Helsgaun (5) states that we should only search for the 5 nearest
This applies the restricted lookahead
"""
def get_y_i_candidates(X, Y, G_i_with_x_i, t_last, original_tour, max_count = 5, excluded_nodes = None):
    # Don't include the current node itself
    if excluded_nodes == None:
        excluded_nodes = set()
    else:
        excluded_nodes = set(excluded_nodes)

    if t_last not in excluded_nodes:
        excluded_nodes.add(t_last)

    potential_gains = []
    original_tour_edges = set(get_edge_list_from_tour(original_tour))

    # Check the neighbours of the endpoint
    for neighbour, weight in enumerate(dist_matrix[t_last]):
        if neighbour in excluded_nodes:
            continue

        # Construct the potential y_i and consider its effect on G_i
        y_i_candidate = tuple(sorted((t_last, neighbour)))
        gain = G_i_with_x_i - weight

        # 4d - Gain Criterion (total gain > 0)
        if gain <= 0:
            continue

        # 4c (y_i cannot be an edge previously broken)
        if y_i_candidate in X:
            continue

        # We don't want to repeat the edge otherwise we will end up with len(X) != len(Y)
        # then we will never be able to make a tour
        if y_i_candidate in Y or y_i_candidate in original_tour_edges:
            continue

        # 4e (check that the feasibility criterion is satisfied for i + 1)
        adjacent_edges = get_adjacent_edges(original_tour, neighbour)
        possible = False

        # The original paper states to try and maximise |x_(i +1)| - |y_i|
        # however in my own testing this negatively impacted the testing and performance
        # as such I decided to use |x_i| - |y_i| instead
        for edge in adjacent_edges:
            if edge in X or edge in Y:
                continue

            possible = True
            break

        if possible:
            potential_gains.append((neighbour, gain))

    # Return the maximum amount requested sorted by their gain
    # At this point, all edges in the list should be valid y_i's
    return sorted(potential_gains, key=lambda k: k[1], reverse=False)[:max_count]

"""
Given the existing tour as a list of nodes and the sets X and Y
this will construct a tour or tell you it is invalid
**OLD: This was previously used for asymmetric version but the results were much worse**
"""
def construct_tour(tour, X, Y):
    edges = set(get_edge_list_from_tour(tour))
    # Use set operations to find the new edge list
    new_edges = (edges - X) | Y

    node_sequence = dict()

    for edge in new_edges:
        start, to = edge
        # Start is already pointing to a node
        if start in node_sequence.keys():
            return False, []

        node_sequence[start] = to

    assert len(tour) == num_cities

    # Don't have all of the original nodes
    if len(node_sequence.keys()) != len(tour):
        return False, []

    # Now we can reconstruct the tour
    # We need to make sure we don't have any loops

    new_tour = []
    current_node = 0

    while len(new_tour) != len(tour):
        # We have a loop
        if current_node in new_tour:
            return False, []

        new_tour.append(current_node)
        current_node = node_sequence[current_node]

    # Check that the end node points to the start
    if node_sequence[new_tour[-1]] != new_tour[0]:
        return False, []

    return True, new_tour

"""
Given a tour of nodes and the sets X and Y this constructs a tour
The edges in X, Y satisify (a, b): b > a
The tour will always start at 0
"""
def symmetric_construct_tour(tour, X, Y):
    edges = set(get_edge_list_from_tour(tour))
    # Compute the new edge list using set operations
    # edges - X = remove the broken edges
    # | Y = add the replacement edges
    new_edges = (edges - X) | Y

    node_sequence = dict()
    occurences = {node: 0 for node in tour}

    for edge in new_edges:
        start, end = edge

        # We have a singular loop
        if start == end:
            return False, []

        occurences[start] += 1
        occurences[end] += 1

        # Start is already pointing to a node so add it to the list
        if start in node_sequence.keys():
            node_sequence[start].add(end)
            continue

        # Create a fresh list
        node_sequence[start] = set([end])

    # We need to make sure that each node appears exactly twice
    for node in occurences:
        if occurences[node] != 2:
            return False, []

    new_tour = []
    current_node = 0

    while len(new_tour) != len(tour) - 1:
        new_tour.append(current_node)

        # Check if the current node has adjacent nodes
        if current_node in node_sequence.keys():
            # If it does then take the lowest one as the next node
            if len(node_sequence[current_node]) != 0:
                next = min(node_sequence[current_node])
                node_sequence[current_node].remove(next)
                current_node = next
                continue
            else:
                # If we've used all adjacent nodes then remove it
                del node_sequence[current_node]

        found = False

        # If it was not a dict key then it might be in a set within the dict
        for next in node_sequence.keys():
            if current_node in node_sequence[next]:
                node_sequence[next].remove(current_node)
                current_node = next
                found = True
                break

        # If we couldn't find it then we have a problem
        if not found:
            return False, []

    # Final node added to the tour
    new_tour.append(current_node)
    return True, new_tour

# Seen hashes is used to avoid checkout time (Helsgaun, Rule 7)
seen_hashes = set()
# Stores up to the 5 most recent optima to compute the common edges
recent_optima = Queue()
# For i >= 4 we don't break any edges in this set
common_edges = set()

def pick_x_i(oldX, oldY, G_i, t_1, t_last, original_tour):
    global seen_hashes, common_edges

    i = len(oldX) + 1
    # Get the edges connected to the node, should always be 2
    t_last_adjacent_edges = get_adjacent_edges(original_tour, t_last)

    # When i >= 4 we apply reduction (Helsgaun, Rule 6)
    # This means that we don't break seemingly important edges for i >= 4
    # because it is expensive
    if i >= 4:
        allowed_edges = []

        left_edge = t_last_adjacent_edges[0]
        left_edge_breakable = not(left_edge in common_edges)

        right_edge = t_last_adjacent_edges[1]
        right_edge_breakable = not(right_edge in common_edges)

        # Only allow if they are not common edges
        if left_edge_breakable:
            allowed_edges.append(left_edge)

        if right_edge_breakable:
            allowed_edges.append(right_edge)

        # For i = 4 (and when we have two choices)
        # we only take the highest weight edge
        # (Helsgaun, Rule 9)
        if len(allowed_edges) == 2 and i == 4:
            allowed_edges = []

            left_weight = get_edge_weight(left_edge)
            right_weight = get_edge_weight(right_edge)

            if left_weight > right_weight:
                allowed_edges.append(left_edge)
            else:
                allowed_edges.append(right_edge)

        t_last_adjacent_edges = allowed_edges

    # Step 6b says we should try the alternative edge if we can't improve
    # with our initial selection
    for adjacent_edge in t_last_adjacent_edges:
        # Duplicate the sets so we don't unintentionally alter something
        # we don't mean to
        X = set(oldX)
        Y = set(oldY)

        # Get t_2i from the edge
        t_2i = adjacent_edge[0] if adjacent_edge[0] != t_last else adjacent_edge[1]
        x_i = adjacent_edge[:]

        # We must maintain that X and Y are disjoint
        # Do not break links that we have added in Y
        if x_i in Y:
            continue

        X.add(x_i)

        # Helsgaun (Rule 2) check feasibility
        y_i_star = tuple(sorted((t_2i, t_1)))

        if y_i_star in Y:
            continue

        # Maintain disjointness
        # If we already broke this edge then we can't re-add it
        if y_i_star in X:
            continue

        # Also checks that x_i != y_i_star since x_i is a member of X
        Y.add(y_i_star)
        is_closed, closed_tour = symmetric_construct_tour(original_tour, X, Y)

        # This tour cannot be closed so move on
        # Only for i >= 3 (Helsgaun, Rule 3)
        if not is_closed and i >= 3:
            continue

        # This avoids checkout time
        # No point searching further if we already explored this path
        if is_closed:
            # Do this via hashing
            # Checking for hashes in a set is quick
            solution_hash = hash(tuple(closed_tour))

            if solution_hash in seen_hashes:
                continue

            seen_hashes.add(solution_hash)

        x_i_weight = get_edge_weight(x_i)
        y_i_star_weight = get_edge_weight(y_i_star)

        # Compute the total gain. It must be > 0 by the Gain Criterion
        y_i_star_G_i = G_i + x_i_weight - y_i_star_weight

        # We reconstruct the tour on return instead of using closed_tour
        # this is because it's possible that for i = 2 we have a closed tour too
        if y_i_star_G_i > 0 and is_closed:
            Y.add(y_i_star)
            # We have a valid tour
            return symmetric_construct_tour(original_tour, X, Y)

        # If we are going to keep going remove the edge that closes up the tour
        Y.remove(y_i_star)
        # We will now choose the y_i
        # We pass the current X, Y and the gain with the effect of x_i
        # Then we continue passing the first node, the most recent node and the original tour we are optimising
        improved, new_tour = pick_y_i(X, Y, G_i + x_i_weight, t_1, t_2i, original_tour)

        return improved, new_tour

    # If we couldn't find a good x_i to break then we didn't improve the solution
    return False, []

def pick_y_i(oldX, oldY, G_i_with_x_i, t_1, t_last, original_tour):
    # This was suggested by the Arthur Maheo article
    # Altering these values will affect the quality of solutions
    # We consider more for the lower values of i as step 6a says
    # that we consider them in order of increasing length if we don't improve
    max_candidates = 5 if len(oldX) <= 3 else 1
    candidates = get_y_i_candidates(oldX, oldY, G_i_with_x_i, t_last, original_tour, max_count = max_candidates)

    for candidate, gain in candidates:
        # Create the y_i edge
        y_i = tuple(sorted((t_last, candidate)))

        X = set(oldX)
        Y = set(oldY)

        # No need to check for disjointness etc since we do all of that in the
        # get_y_i_candidates function
        Y.add(y_i)

        # Once we have a y_i we go back to pick x_(i + 1) or see if we can close the tour
        improved, new_tour = pick_x_i(X, Y, gain, t_1, candidate, original_tour)

        if improved:
            return improved, new_tour

    return False, []

def iterate_lk(tour):
    # Step 7 says to try all t_1 (step 6e is also relevant)
    for t_1 in tour:
        t_1_adjacent_edges = get_adjacent_edges(tour, t_1)

        # Step 6d says we try the alternate t_2 if the original didn't lead to improvement
        for adjacent_edge in t_1_adjacent_edges:
            t_2 = adjacent_edge[0] if adjacent_edge[0] != t_1 else adjacent_edge[1]
            x_1 = adjacent_edge[:]
            x_1_weight = get_edge_weight(x_1)

            # If we didn't improve with the original best then we try the next best
            # and so on
            for t_3, y_1_weight in get_best_neighbours(tour, t_2, 5):
                y_1 = tuple(sorted((t_2, t_3)))
                g_1 = x_1_weight - y_1_weight

                # Step 3 requires this
                if g_1 <= 0:
                    continue

                X = set()
                X.add(x_1)
                Y = set()
                Y.add(y_1)

                # Pick an x_2 and start the whole process of optimising
                improved, new_tour = pick_x_i(X, Y, g_1, t_1, t_3, tour[:])

                if improved:
                    return improved, new_tour

    return False, []

def update_common_edges(new_solution):
    global common_edges, recent_optima
    # This is used to apply the reduction heuristic

    # We only want distinct optima
    if new_solution in recent_optima.queue:
        return

    # We require at least 2 optima to do this with
    # so for the first one we can't check the intersection
    if recent_optima.qsize() == 0:
        recent_optima.put(new_solution)
        return

    # We limit ourselves to the 5 most recent distinct optima
    if recent_optima.qsize() >= 5:
        recent_optima.get(block = False)

    recent_optima.put(new_solution)

    solutions = list(recent_optima.queue)
    shared = solutions[0]

    # & is the intersection operator for sets
    for solution in solutions[1:]:
        shared &= solution

    common_edges = shared

"""
References:
1) Lin, Shen, and Brian W. Kernighan. "An effective heuristic algorithm for the traveling-salesman problem." Operations research 21.2 (1973): 498-516. (https://pdfs.semanticscholar.org/88c3/ae44f61301aa2974f4e65f73d17f5944c0bb.pdf)
2) Helsgaun, Keld. "An effective implementation of the Lin–Kernighan traveling salesman heuristic." European Journal of Operational Research 126.1 (2000): 106-130. (https://homes.di.unimi.it/righini/Didattica/AlgoritmiEuristici/MaterialeAE/Helsgaun.pdf)
3) Arthur Maheo https://arthur.maheo.net/implementing-lin-kernighan-in-python/

In my comments, 'Steps' refer to the original paper (ref 1)
'Rules' are from Helsgaun's paper (ref 2)
"""
def start_lk(tour):
    improved = True
    best_solution = tour

    while improved:
        # When we improve the solution we restart the iteration with the improvement
        improved, possible_solution = iterate_lk(best_solution)

        if improved:
            update_common_edges(set(get_edge_list_from_tour(possible_solution)))
            best_solution = possible_solution

    return best_solution, tour_length_calc(best_solution)

def two_opt(start_solution):
    best_solution = start_solution[:]
    best_score = tour_length_calc(best_solution)
    improved = True

    while improved:
        improved = False

        for i in range(0, num_cities):
            for k in range(i + 1, num_cities):
                new_solution = best_solution[:i] + best_solution[i:k + 1][::-1] + best_solution[k + 1:]
                new_score = tour_length_calc(new_solution)

                if new_score < best_score:
                    improved = True
                    best_solution = new_solution
                    best_score = new_score

    return best_solution, best_score

base_tour = [x for x in range(num_cities)]
#base_tour = [220,521,97,131,403,410,204,499,207,217,5,294,301,430,235,309,332,413,10,257,282,65,3,339,9,492,359,42,104,440,212,69,72,200,457,498,515,476,78,125,458,73,512,363,87,24,378,157,234,213,265,534,305,98,321,297,351,358,163,325,241,340,51,356,120,348,166,269,424,336,46,165,95,88,32,429,431,448,134,158,194,379,530,449,176,441,243,21,118,361,156,447,231,380,0,67,376,227,382,407,102,354,38,327,360,510,489,258,54,478,278,421,151,398,189,195,143,178,43,443,455,337,259,460,14,428,94,240,395,292,188,190,250,185,419,412,468,288,58,145,383,394,36,475,179,40,437,409,517,81,329,426,254,526,70,436,106,59,318,525,255,514,52,101,2,390,245,304,293,153,408,459,7,502,528,507,400,298,100,119,272,275,276,82,296,404,37,501,286,141,152,35,84,444,183,30,113,267,22,450,364,283,186,191,225,238,25,150,374,438,393,268,389,174,420,375,132,485,114,367,175,6,77,148,161,427,506,232,423,312,405,343,371,197,532,402,322,162,181,471,17,320,456,111,260,347,505,491,529,352,500,1,117,467,509,334,344,103,435,33,274,414,422,55,19,126,50,472,99,406,392,279,503,520,316,110,328,388,47,252,149,142,196,396,261,154,246,342,170,115,432,108,211,357,333,262,139,366,159,511,18,480,192,133,62,284,300,13,522,496,28,417,249,416,533,56,239,128,8,26,222,397,487,137,138,433,488,373,29,135,483,484,53,20,497,123,387,187,229,214,461,228,216,86,495,41,355,313,266,130,68,271,122,140,264,504,182,109,253,147,465,518,263,4,353,307,236,173,331,247,202,401,399,287,486,215,218,112,474,146,345,469,27,311,49,451,323,66,425,45,513,386,219,303,324,96,16,92,418,233,116,237,290,201,168,167,338,302,76,481,350,446,93,462,85,107,519,75,121,494,277,164,439,477,370,90,74,464,63,346,299,490,209,223,124,11,89,206,248,527,60,369,15,289,335,180,372,224,415,326,244,34,144,39,91,445,61,341,452,48,208,210,280,524,127,198,310,71,129,64,508,454,184,230,330,411,319,31,493,251,171,12,368,57,381,453,169,44,79,160,155,470,315,466,83,242,23,270,205,136,391,177,193,349,516,377,172,256,531,482,473,221,226,105,384,314,442,281,463,199,434,273,285,479,308,295,385,306,291,523,365,317,80,362,203]

print("Starting...")
# print(base_tour)
# print(get_edge_list_from_tour(base_tour))
# print(get_adjacent_edges(base_tour, 8))

# print(dist_matrix[5])
# print(get_best_neighbours(base_tour, 5, 7))

start_time = time.time()
tour, tour_length = start_lk(base_tour)
end_time = time.time()

print(tour, tour_length)
print("Took", end_time - start_time, "s")

import sys
sys.exit(0)

############
############ YOUR CODE SHOULD NOW BE COMPLETE AND WHEN EXECUTION OF THIS PROGRAM 'skeleton.py'
############ REACHES THIS POINT, YOU SHOULD HAVE COMPUTED A TOUR IN THE RESERVED LIST VARIABLE 'tour',
############ WHICH HOLDS A LIST OF THE INTEGERS FROM {0, 1, ..., 'num_cities' - 1}, AND YOU SHOULD ALSO
############ HOLD THE LENGTH OF THIS TOUR IN THE RESERVED INTEGER VARIABLE 'tour_length'.
############

############
############ YOUR TOUR WILL BE PACKAGED IN A TOUR FILE OF THE APPROPRIATE FORMAT AND THIS TOUR FILE,
############ WHOSE NAME WILL BE A MIX OF THE NAME OF THE CITY FILE, THE NAME OF THIS PROGRAM AND THE
############ CURRENT DATA AND TIME. SO, EVERY SUCCESSFUL EXECUTION GIVES A TOUR FILE WITH A UNIQUE
############ NAME AND YOU CAN RENAME THE ONES YOU WANT TO KEEP LATER.
############

############
############ DO NOT TOUCH OR ALTER THE CODE BELOW THIS POINT! YOU HAVE BEEN WARNED!
############

flag = "good"
length = len(tour)
for i in range(0, length):
    if isinstance(tour[i], int) == False:
        flag = "bad"
    else:
        tour[i] = int(tour[i])
if flag == "bad":
    print("*** error: Your tour contains non-integer values.")
    sys.exit()
if isinstance(tour_length, int) == False:
    print("*** error: The tour-length is a non-integer value.")
    sys.exit()
tour_length = int(tour_length)
if len(tour) != num_cities:
    print("*** error: The tour does not consist of " + str(num_cities) + " cities as there are, in fact, " + str(len(tour)) + ".")
    sys.exit()
flag = "good"
for i in range(0, num_cities):
    if not i in tour:
        flag = "bad"
if flag == "bad":
    print("*** error: Your tour has illegal or repeated city names.")
    sys.exit()
check_tour_length = 0
for i in range(0, num_cities - 1):
    check_tour_length = check_tour_length + dist_matrix[tour[i]][tour[i + 1]]
check_tour_length = check_tour_length + dist_matrix[tour[num_cities - 1]][tour[0]]
if tour_length != check_tour_length:
    flag = print("*** error: The length of your tour is not " + str(tour_length) + "; it is actually " + str(check_tour_length) + ".")
    sys.exit()
print("You, user " + my_user_name + ", have successfully built a tour of length " + str(tour_length) + "!")

local_time = time.asctime(time.localtime(time.time()))
output_file_time = local_time[4:7] + local_time[8:10] + local_time[11:13] + local_time[14:16] + local_time[17:19]
output_file_time = output_file_time.replace(" ", "0")
script_name = os.path.basename(sys.argv[0])
if len(sys.argv) > 2:
    output_file_time = sys.argv[2]
output_file_name = script_name[0:len(script_name) - 3] + "_" + input_file[0:len(input_file) - 4] + "_" + output_file_time + ".txt"

f = open(output_file_name,'w')
f.write("USER = " + my_user_name + " (" + my_first_name + " " + my_last_name + "),\n")
f.write("ALGORITHM CODE = " + algorithm_code + ", NAME OF CITY-FILE = " + input_file + ",\n")
f.write("SIZE = " + str(num_cities) + ", TOUR LENGTH = " + str(tour_length) + ",\n")
f.write(str(tour[0]))
for i in range(1,num_cities):
    f.write("," + str(tour[i]))
f.write(",\nNOTE = " + added_note)
f.close()
print("I have successfully written your tour to the tour file:\n   " + output_file_name + ".")
