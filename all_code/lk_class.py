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

def tour_length_calc(state):
    dist = dist_matrix[state[0]][state[len(state) - 1]]

    for i in range(0, len(state) - 1):
        dist += dist_matrix[state[i]][state[i + 1]]

    return dist

class LKUtil:
    g_star = 0
    g_star_X = None
    g_star_Y = None
    g_star_tour = None
    seen_hashes = set()

    @staticmethod
    def convert_tour_to_edges(tour):
        """
        Returns a tour of edges from a tour of nodes
        Should have the same length as the original tour
        We are using sets so we can do set operations and they are fast
        """
        return set([(tour[i], tour[(i + 1) % len(tour)]) for i in range(len(tour))])

    @staticmethod
    def get_length_of_edge_tour(edge_tour):
        """
        Gets the length of a tour of edges
        """
        return sum([dist_matrix[a][b] for a, b in edge_tour])

    @staticmethod
    def convert_edge_tour_to_nodes(edge_tour):
        """
        Gets the tour as an array of nodes
        It also returns whether the tour is valid
        If it is not valid the returned array is empty
        returns in the order <valid>, <tour>
        Tour will always start at 0 if it's valid
        """

        node_sequence = {}

        for a, b in edge_tour:
            if a in node_sequence.keys():
                return False, []

            node_sequence[a] = b

        if len(node_sequence.keys()) != num_cities:
            return False, []

        node_tour = []
        current_node = 0

        while len(node_tour) != num_cities:
            node_tour.append(current_node)
            current_node = node_sequence[current_node]

        if len(set(node_tour)) != num_cities:
            return False, []

        return True, node_tour

    @staticmethod
    def get_best_neighbours(node, max_count, exclude_nodes = None):
        """
        Returns the max_count best neighbours of node excluding exclude_nodes
        Returns tuples an array (in sorted order) of tuples of the form (<node>, <weight>)
        """
        if exclude_nodes == None:
            exclude_nodes = set()

        # Don't return to itself
        if node not in exclude_nodes:
            exclude_nodes.add(node)

        weights = dict()

        for neighbour, weight in enumerate(dist_matrix[node]):
            if neighbour in exclude_nodes:
                continue

            weights[neighbour] = weight

        return [(neighbour, weights[neighbour]) for neighbour in sorted(weights, key=weights.get, reverse=False)[:max_count]]

    @staticmethod
    def apply_lookahead(instance, max_count, exclude_nodes = None):
        if exclude_nodes == None:
            exclude_nodes = set()

        if instance.t_last not in exclude_nodes:
            exclude_nodes.add(instance.t_last)

        weights = dict()
        current_index = instance.original_nodes.index(instance.t_last)

        for neighbour, weight in enumerate(dist_matrix[instance.t_last]):
            if len(weights.keys()) == max_count:
                break

            if neighbour in exclude_nodes:
                continue

            neighbour_index = instance.original_nodes.index(neighbour)
            y_i = (neighbour, instance.t_last) if neighbour_index < current_index else (instance.t_last, neighbour)
            y_i_weight = dist_matrix[y_i[0]][y_i[1]]

            if instance.G_i - y_i_weight <= 0 or y_i in instance.X or y_i in instance.Y or y_i in ((instance.original_T - instance.X) | instance.Y):
                continue

            neighbour_adjacent_nodes = [instance.original_nodes[neighbour_index - 1], instance.original_nodes[(neighbour_index + 1) % num_cities]]

            for t_next_is_after, t_next in enumerate(neighbour_adjacent_nodes):
                next_x_i = (t_next, neighbour) if t_next_is_after == 0 else (neighbour, t_next)
                next_x_i_weight = dist_matrix[next_x_i[0]][next_x_i[1]]

                if next_x_i in instance.X or next_x_i in instance.Y:
                    continue

                # May need to reverse?
                weights[neighbour] = next_x_i_weight - y_i_weight

        #print([(neighbour, weights[neighbour]) for neighbour in sorted(weights, key=weights.get, reverse=True)[:max_count]])
        return [(neighbour, weights[neighbour]) for neighbour in sorted(weights, key=weights.get, reverse=True)[:max_count]]

# Will act as a single 'thread' of the LK chain we form
# Each branch will form a new thread and thus a new LK instance
class LKInstance:
    def __init__(self, X, Y, G_i, original_nodes, original_T, t_1, t_last, last_weight):
        """
        Create a new instance
        Shallow copies X and Y so we can make changes
        We also store t_1 and most recent end point t_last as well as the original edge and node tour
        i is deduced as len(X)
        G_i is the sum of the gains so far (G_3 = g_1 + g_2 + g_3 etc.)
        """
        self.X = set(X)
        self.Y = set(Y)
        self.G_i = G_i

        self.original_nodes = original_nodes
        self.original_T = original_T
        self.t_1 = t_1
        self.t_last = t_last
        self.last_weight = last_weight

        self.i = len(self.X)

    def solve(self):
        pass

    def pick_x_i(self):
        """
        When we pick x_i we need to look at the adjacent nodes of t_last
        We require that:
        1. 'x_i is chosen so that if t_2i is joined to t_1, the resulting configuration is a tour' (S.Lin pg502)#
        2. x_i is not in X or Y
        """

        if len(self.X) != len(self.Y):
            print(self.X, self.Y)

        # x_i consists of (t_(2i - 1), t_2i), y_i consists of (t_2i, t_(2i + 1))

        # t_last is actually t_(2i - 1)
        t_last_position = self.original_nodes.index(self.t_last)
        t_last_adjacent_nodes = [self.original_nodes[t_last_position - 1], self.original_nodes[(t_last_position + 1) % num_cities]]
        #print("l", len(t_last_adjacent_nodes))

        # For x_4 we always take the longest edge to remove
        if self.i == 4:
            #print(self.X, self.Y)
            left_edge = (t_last_adjacent_nodes[0], self.t_last)
            left_weight = dist_matrix[left_edge[0]][left_edge[1]]

            right_edge = (self.t_last, t_last_adjacent_nodes[1])
            right_weight = dist_matrix[right_edge[0]][right_edge[1]]

            if right_weight > left_weight:
                t_last_adjacent_nodes = [t_last_adjacent_nodes[1]]
            else:
                t_last_adjacent_nodes = [t_last_adjacent_nodes[0]]

        for t_2i_is_after, t_2i in enumerate(t_last_adjacent_nodes):
            x_i = (t_2i, self.t_last) if t_2i_is_after == 0 else (self.t_last, t_2i)

            # We can't break edges we have already broke or have added ourselves
            if x_i in self.X or x_i in self.Y:
                continue

            #print(self.i, x_i)
            x_i_weight = dist_matrix[x_i[0]][x_i[1]]
            G_i_no_y_i = self.G_i + x_i_weight

            # Before we go on to construct y_i step 4f says we check if closing up
            # will give us a better gain value than the best we have already seen

            feasibility_test_edge = (t_2i, self.t_1)

            if feasibility_test_edge in self.X or feasibility_test_edge in self.Y:
                # We've already got this edge so it can't be put in again
                continue

            # We need to make sure we can make a tour if we remove x_i
            # This uses set operations to construct the tour
            # When i == 2 we have special behaviour which is from Helsgaun (3)
            feasibility_tour = (self.original_T - self.X - set([x_i])) | (self.Y | set([feasibility_test_edge]))

            tour_hash = hash(frozenset(feasibility_tour))

            if tour_hash in LKUtil.seen_hashes:
                return False

            LKUtil.seen_hashes.add(tour_hash)

            is_valid_tour = LKUtil.convert_edge_tour_to_nodes(feasibility_tour)[0]

            # looking_for = {(0, 10), (10, 5), (5, 6), (6, 8), (8, 2), (2, 1), (1, 7), (7, 3), (3, 11), (11, 4), (4, 9), (9, 0)}

            # if looking_for == feasibility_tour:
            #     print("IT'S HERE")

            # The first part returned is a boolean stating if the tour is valid
            if not is_valid_tour and self.i > 2:
                #print(feasibility_tour)
                #print("Not valid :(", self.i)
                #print("Nope", len(self.X), len(self.Y))
                continue
            # TODO: Backtracking see step 6

            if self.i > 2:
                print("Valid! :)", self.i, len(self.X))
                print(feasibility_tour)

            #print(len(self.X), len(self.Y))


            # We've found a better tour than we have previously so lets save it
            if is_valid_tour:
                feasibility_test_edge_weight = dist_matrix[feasibility_test_edge[0]][feasibility_test_edge[1]]
                if G_i_no_y_i - feasibility_test_edge_weight > 0:
                    self.X.add(x_i)
                    self.Y.add(feasibility_test_edge)

                    if(sum([dist_matrix[a[0]][a[1]] for a in self.X]) < sum([dist_matrix[a[0]][a[1]] for a in self.Y])):
                        continue

                    # print("New tour")
                    # print(self.X, [dist_matrix[a[0]][a[1]] for a in self.X], "-", sum([dist_matrix[a[0]][a[1]] for a in self.X]))
                    # print(self.Y, [dist_matrix[a[0]][a[1]] for a in self.Y], "+", sum([dist_matrix[a[0]][a[1]] for a in self.Y]))
                    # LKUtil.g_star = G_i_no_y_i - feasibility_test_edge_weight
                    LKUtil.g_star_tour = feasibility_tour
                    # print(LKUtil.g_star_tour, LKUtil.g_star, "Actual:", LKUtil.get_length_of_edge_tour(LKUtil.g_star_tour))
                    return True

            #print("Go to y_i")
            new_X = set(self.X)
            new_X.add(x_i)
            y_i_search = LKInstance(new_X, self.Y, G_i_no_y_i, LKUtil.convert_edge_tour_to_nodes(self.original_T)[1], self.original_T, self.t_1, t_2i, x_i_weight).pick_y_i()

            if self.i == 2 and y_i_search:
                return True

            return y_i_search

        return False

    def pick_y_i(self):
        """
        When we pick y_i we require that:
        1. X and Y are disjoint
        2. G_i is positive (note that self.G_i consists of G_(i - 1) + x_i_weight so need to subtract y_i_weight)
        3. At i + 1, the y_i chosen must permit the breaking of an x_(i + 1)
        We already check condition 4f in pick_x_i

        According to Helsgaun (5) we should restrict to the nearest 5 neighbours of t_2i (which is self.t_last here)
        TODO: Implement the lookahead specified in Helsgaun (8)
        """

        #print("Last", self.t_last)

        #t_last_nearest_neighbours = LKUtil.get_best_neighbours(self.t_last, 5 if self.i == 2 else 1)
        t_last_nearest_neighbours = LKUtil.apply_lookahead(self, 5 if self.i == 2 else 1)
        #print(t_last_nearest_neighbours)

        for neighbour, y_i_weight in t_last_nearest_neighbours:
            y_i = (self.t_last, neighbour)

            # Already seen it
            if y_i in self.X or y_i in self.Y:
                #print("Seen it")
                continue

            if self.last_weight - y_i_weight <= 0:
                #print("s")
                continue

            # We require the gain to be positive
            if self.G_i - y_i_weight <= 0:
                #print("Not good")
                continue

            # TODO: Check that at i + 1 we can permit the breaking of x_(i + 1)

            new_Y = set(self.Y)
            new_Y.add(y_i)

            if(LKInstance(self.X, new_Y, self.G_i - y_i_weight, LKUtil.convert_edge_tour_to_nodes(self.original_T)[1], self.original_T, self.t_1, neighbour, None).pick_x_i()):
                #print("Yep")
                return True

        return False


def start_lin_kernighan(original):
    # We first need to turn this tour consisting of n cities in an array into
    # a tour of edges which is needed for LK
    original_tour_edges = LKUtil.convert_tour_to_edges(original)
    LKUtil.g_star_tour = original_tour_edges
    last_g_star_tour = None
    its = 0
    improved = True

    while improved:
        improved = improve()



        #print(LKUtil.g_star)

def improve():
    #print(LKUtil.g_star_tour, LKUtil.get_length_of_edge_tour(LKUtil.g_star_tour))
    LKUtil.g_star = 0
    last_g_star_tour = LKUtil.g_star_tour
    #print(original)
    #print(LKUtil.convert_edge_tour_to_nodes(last_g_star_tour)[1])
    last_g_star_nodes = LKUtil.convert_edge_tour_to_nodes(last_g_star_tour)[1]#

    # Step 7 says we should try all nodes in the original
    for t_1_position, t_1 in enumerate(last_g_star_nodes):
        #print(t_1_position)
        # We need to look at any edge adjacent to t_1 so we can use t_1_position for this
        t_1_adjacent_nodes = [last_g_star_nodes[t_1_position - 1], last_g_star_nodes[(t_1_position + 1) % num_cities]]

        # Step 3
        # Now get t_2 which is the adjacent node
        # Step 6d says we should try both x_1's
        for t_2_is_after, t_2 in enumerate(t_1_adjacent_nodes):
            # t_2_is_after == 0 means it points to t_1, otherwise t_1 points to t_2
            x_1 = (t_2, t_1) if t_2_is_after == 0 else (t_1, t_2)
            x_1_weight = dist_matrix[x_1[0]][x_1[1]]

            attempts = 0

            # g_i = x_i_weight - y_i_weight
            # Step 6d says we should try all y_1's that have g_1 > 0
            # Look at the neighbours of t_2

            for t_3, y_1_weight in LKUtil.get_best_neighbours(t_2, 5):
                if attempts == 5: break
                # Q: Do I need to impose t_3 == t_1 here??
                # We can't have t_3 being t_2 or an adjacent node to t_1 or t_1 itself
                if t_3 in t_1_adjacent_nodes or t_3 == t_2:
                    continue

                g_1 = x_1_weight - y_1_weight

                # <= 0 so skip it
                if g_1 <= 0:
                    continue

                y_1 = (t_2, t_3)

                if x_1 == y_1:
                    continue

                # We now have a gain so lets do LK on it
                # This finishes the setup for the algorithm now we move into the more
                # in depth steps
                X = set()
                X.add(x_1)
                Y = set()
                Y.add(y_1)

                lk_instance = LKInstance(X, Y, g_1, last_g_star_nodes[:], last_g_star_tour, t_1, t_3, None)
                improved = lk_instance.pick_x_i()

                if improved:
                    print("Improved", LKUtil.g_star, LKUtil.get_length_of_edge_tour(LKUtil.g_star_tour))
                    return True

                attempts += 1

    return False

#base = [2, 8, 6, 5, 10, 0, 9, 4, 11, 1, 3, 7]
base = [x for x in range(num_cities)]
#random.shuffle(base)
s = time.time()
start_lin_kernighan(base)
q = time.time() - s
print()
print("T", q, "s")
print(LKUtil.g_star)
print(LKUtil.g_star_tour)
node_tour = LKUtil.convert_edge_tour_to_nodes(LKUtil.g_star_tour)[1]
print(node_tour, tour_length_calc(node_tour))

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
