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
    return [(tour[i], tour[(i + 1) % len(tour)]) for i in range(len(tour))]

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
    return [(tour[position - 1], tour[position]), (tour[position], tour[(position + 1) % len(tour)])]

"""
Returns the 'max_count' best neighbours of 'node' in the 'tour'
Can skip nodes if they are in an iterable in excluded_nodes
A list of tuples. tuple[0] is the neighbour, tuple[1] is the weight
We also don't want the pair to exist in the tour already
"""
def get_best_neighbours(tour, node, max_count, excluded_nodes = None):
    if excluded_nodes == None:
        excluded_nodes = set()
    else:
        excluded_nodes = set(excluded_nodes)

    for adjacent in get_adjacent_nodes(tour, node):
        excluded_nodes.add(adjacent)

    # Don't want to match to itself
    if node not in excluded_nodes:
        excluded_nodes.add(node)

    weights = []

    for neighbour, weight in enumerate(dist_matrix[node]):
        if neighbour in excluded_nodes:
            continue

        weights.append((neighbour, weight))

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
"""
def get_y_i_candidates(X, Y, G_i_with_x_i, t_last, original_tour, max_count = 5, excluded_nodes = None):
    if excluded_nodes == None:
        excluded_nodes = set()
    else:
        excluded_nodes = set(excluded_nodes)

    if t_last not in excluded_nodes:
        excluded_nodes.add(t_last)

    potential_gains = []
    original_tour_edges = set(get_edge_list_from_tour(original_tour))

    # TODO: We may need to check that y_i is not in the original tour edges or in Y

    for neighbour, weight in enumerate(dist_matrix[t_last]):
        if neighbour in excluded_nodes:
            continue

        y_i_candidate = (t_last, neighbour)
        gain = G_i_with_x_i - weight

        # 4d
        if gain <= 0:
            continue

        # 4c
        if y_i_candidate in X:
            continue

        if y_i_candidate in Y or y_i_candidate in (original_tour_edges - X):
            continue

        # 4e
        adjacent_edges = get_adjacent_edges(original_tour, neighbour)
        possible = False

        for edge in adjacent_edges:
            if edge in X or edge in Y:
                continue

            possible = True
            break

        if possible:
            potential_gains.append((neighbour, gain))

    return sorted(potential_gains, key=lambda k: k[1], reverse=False)[:max_count]

"""
Given the existing tour as a list of nodes abd the sets X and Y
this will construct a tour or tell you it is invalid
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

def pick_x_i(oldX, oldY, G_i, t_1, t_last, original_tour):
    assert len(oldX) == len(oldY)

    #print(oldX, oldY)
    i = len(oldX) + 1
    # t_last is t_(2i-1)

    t_last_adjacent_edges = get_adjacent_edges(original_tour, t_last)

    if i == 4:
        left_weight = get_edge_weight(t_last_adjacent_edges[0])
        right_weight = get_edge_weight(t_last_adjacent_edges[1])

        if left_weight > right_weight:
            t_last_adjacent_edges = [t_last_adjacent_edges[0]]
        else:
            t_last_adjacent_edges = [t_last_adjacent_edges[1]]

    for is_next, adjacent_edge in enumerate(t_last_adjacent_edges):
        X = set(oldX)
        Y = set(oldY)

        t_2i = adjacent_edge[is_next]
        x_i = adjacent_edge[:]

        # We must maintain that X and Y are disjoint
        if x_i in Y:
            continue

        X.add(x_i)

        # Helsgaun (2) check feasibility
        y_i_star = (t_2i, t_1)

        if y_i_star in Y:
            continue

        # Maintain disjointness

        if y_i_star in X:
            continue

        # Also checks that x_i != feasibility_criterion_edge
        Y.add(y_i_star)
        is_closed, closed_tour = construct_tour(original_tour, X, Y)

        # This tour cannot be closed so move on
        if not is_closed and i >= 3:
            if i > 3:
                print(":(")
                print(original_tour)
                print(t_1)
                print(X)
                print(Y)
            continue

        if i>3 and is_closed:
            print("!!!")
            # print(original_tour, X, Y)

        x_i_weight = get_edge_weight(x_i)
        y_i_star_weight = get_edge_weight(y_i_star)

        y_i_star_G_i = G_i + x_i_weight - y_i_star_weight

        if y_i_star_G_i > 0 and is_closed:
            Y.remove(y_i_star)
            improved, new_tour = pick_y_i(X, Y, G_i + x_i_weight, t_1, t_2i, original_tour)

            if improved:
                return improved, new_tour

            Y.add(y_i_star)
            # We have a valid tour
            return construct_tour(original_tour, X, Y)

        # If we are going to keep going remove the finishing edge
        Y.remove(y_i_star)
        # We will now choose the y_i
        # We may need to impose additional conditions here
        #print(X, Y)
        improved, new_tour = pick_y_i(X, Y, G_i + x_i_weight, t_1, t_2i, original_tour)

        return improved, new_tour

    return False, []

def pick_y_i(oldX, oldY, G_i_with_x_i, t_1, t_last, original_tour):
    assert len(oldX) == len(oldY) + 1
    candidates = get_y_i_candidates(oldX, oldY, G_i_with_x_i, t_last, original_tour)
    #print(candidates)

    for candidate, gain in candidates:
        y_i = (t_last, candidate)

        X = set(oldX)
        Y = set(oldY)

        r = len(Y)
        Y.add(y_i)
        assert len(Y) != r

        improved, new_tour = pick_x_i(X, Y, gain, t_1, candidate, original_tour)

        if improved:
            return improved, new_tour

    return False, []

def iterate_lk(tour):
    # Step 7 says to try all t_1
    for t_1 in tour:
        t_1_adjacent_edges = get_adjacent_edges(tour, t_1)
        t_1_adjacent_nodes = get_adjacent_nodes(tour, t_1)
        assert len(t_1_adjacent_edges) == 2

        for is_next, adjacent_edge in enumerate(t_1_adjacent_edges):
            # is_next will be 0 if we have (prev, t_1) and 1 if we have (t_1, next)
            t_2 = adjacent_edge[is_next]
            x_1 = adjacent_edge[:]
            x_1_weight = get_edge_weight(x_1)

            # We may need to increase this or change how it works
            for t_3, y_1_weight in get_best_neighbours(tour, t_2, 5):
                y_1 = (t_2, t_3)
                g_1 = x_1_weight - y_1_weight

                # Step 3 requires this
                if g_1 <= 0:
                    continue

                X = set()
                X.add(x_1)
                Y = set()
                Y.add(y_1)

                improved, new_tour = pick_x_i(X, Y, g_1, t_1, t_3, tour[:])

                if improved:
                    return improved, new_tour

    return False, []

def start_lk(tour):
    # TODO: May need to zero the tour
    improved = True
    best_solution = tour#two_opt(tour)[0]

    while improved:
        improved, possible_solution = iterate_lk(best_solution)

        if improved:
            #print(improved, best_solution, tour_length_calc(best_solution))
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

# random.shuffle(base_tour)
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
