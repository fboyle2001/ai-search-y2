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

def tour_length_from_edges(T):
    dist = 0

    for pair in T:
        dist += dist_matrix[pair[0]][pair[1]]

    return dist

def tour_to_edges(T):
    return set([(T[i], T[i + 1]) for i in range(len(T) - 1)] + [(T[-1], T[0])])

def convert_edge_list_to_tour(edge_set):
    #if len(edge_set) != num_cities:
    #    #print("Not enough edges for a full tour", len(edge_set), num_cities)
    #    return None

    after = dict()

    for start, end in edge_set:
        if start in after.keys():
            return None

        after[start] = end

    if len(after.keys()) != len(edge_set):
        #print("Rep edges not enough")
        return None

    T = []
    current = 0

    #print(edge_set, after)

    while len(T) != len(edge_set):
        T.append(current)
        current = after[current]

    return T

def is_valid_tour(T):
    return T != None and len(set(T)) == num_cities and len(T) == num_cities

# T = tour, t_1 is the start node, t_last = t_3/t_5/...
# can deduce i from len(X)
# X = edges we broke, Y = edges we are adding
def pick_x_i(nodes, T, t_1, t_last, X, Y, total_gain):
    i = len(X)
    #print(i)
    t_last_index = nodes.index(t_last)
    #print(t_last_index)
    direct_nodes = (nodes[t_last_index - 1], nodes[(t_last_index + 1) % num_cities])

    # x_i consists of (t_(2i - 1), t_2i), y_i consists of (t_2i, t_(2i + 1))

    #print(T, tour_length_from_edges(T))

    for is_after, t_2i in enumerate(direct_nodes):
        if t_2i == t_last:
            continue

        x_i = (t_last, t_2i) if is_after == 1 else (t_2i, t_last)
        # The new x_i must satisfy:
        # 1. Connecting t_2i to t_1 results in a tour (for i >= 3) (https://homes.di.unimi.it/righini/Didattica/AlgoritmiEuristici/MaterialeAE/Helsgaun.pdf, pg 4)
        # 2. Not already in X or Y
        # Check (2)
        if x_i in X or x_i in Y:
            continue

        cost_x_i = dist_matrix[x_i[0]][x_i[1]]
        temp_total_gain = total_gain + cost_x_i

        complete_tour_edge = (t_2i, t_1)

        if complete_tour_edge in X or complete_tour_edge in Y:
            continue

        complete_tour_edge_cost = dist_matrix[t_2i][t_1]

        #print("Adding", (Y | set([complete_tour_edge])), len(Y | set([complete_tour_edge])))
        #print("Removing", (X | set([x_i])), len(X | set([x_i])))

        complete_tour = (T - X - set([x_i])) | (Y | set([complete_tour_edge]))
        #print(nodes_complete_tour)
        #print(complete_tour)

        #print("pre", T, X, Y)
        nodes_complete_tour = convert_edge_list_to_tour(complete_tour)
        valid_tour = is_valid_tour(nodes_complete_tour)

        #if(valid_tour):
            #print("ct", complete_tour, len(complete_tour), valid_tour)

        # Check (1)
        # if i >= 3 and not valid_tour:
        #     # Try the other node instead
        #     continue

        X.add(x_i)

        if temp_total_gain - complete_tour_edge_cost > 0 and valid_tour:
            Y.add(complete_tour_edge)
            # complete tour is better so use it

            #print("Returning", complete_tour)
            return complete_tour, tour_length_from_edges(complete_tour)

        # otherwise find an edge to replace x_i we just removed
        #print("yes")
        return pick_y_i(nodes, T, t_1, t_2i, X, Y, temp_total_gain)

        if q != None:
            print(q[1], "\t" * len(X), q)

        return q

def pick_y_i(nodes, T, t_1, t_last, X, Y, total_gain):
    # Must satisify:
    # 1. Not in X or Y (c)
    # 2. total_gain + |x_i| - |y_i| > 0 (d)
    # 3. y_i must permit a breaking of an x_(i+1) (e)
    # We also need to check (f)

    # We will use nearest neighbour to select the y_i edge
    neighbours = dist_matrix[t_last]

    complete_tour_edge = (t_last, t_1)
    #complete_tour_edge_cost = dist_matrix[t_2i][t_1]
    complete_tour = (T - X) | (Y | set([complete_tour_edge]))

    candidates = []
    a = 0

    for neighbour, distance in enumerate(neighbours):
        y_i = (t_last, neighbour)

        if t_last == neighbour:
            continue

        # Check (1)
        if y_i in X or y_i in Y:
            continue

        # Check (2)
        if total_gain - distance <= 0:
            continue

        # Check (f)

        dupe_X = set([x for x in X])
        dupe_Y = set([y for y in Y])

        dupe_Y.add(y_i)

        q = pick_x_i(nodes, T, t_1, neighbour, dupe_X, dupe_Y, total_gain - distance)

        return q

        if q != None:
            print(a, q[0], q[1], len(X), is_valid_tour(q[0]))

        if a == 1 or (a == 0 and q != None and q[1] == 69):
            return q

        a += 1

    #print("Nothing", complete_tour)
    #print("Rem", X)
    #print("Add", Y)
    #print("Additional Edge", complete_tour_edge)
    #return complete_tour, tour_length_from_edges(complete_tour)

def lk(start_solution):
    # Step 1
    nodes = start_solution[:]
    #print(nodes)
    # Write the tour as a set of its edges instead
    T = tour_to_edges(nodes)
    #print(start_solution, T, convert_edge_list_to_tour(T), is_valid_tour(convert_edge_list_to_tour(T)))

    # Step 2
    # Adapt the G* = 0 to be a measure of the best solution found so far instead
    f_T = tour_length_from_edges(T)
    improved = True
    its = 0

    while improved:
        its += 1
        improved = False
        improved_tour, improved_score = T, f_T
        #print(improved_tour)
        nodes = convert_edge_list_to_tour(improved_tour)
        #print(nodes)
        #print("Starting with", improved_tour, improved_score)

        # Pick a t_1 - Links with Step 7
        for index, t_1 in enumerate(nodes):
            #print(index)
            # t_1 is the first node
            # t_1 has a node before and after it
            direct_nodes = (nodes[index - 1], nodes[(index + 1) % num_cities])

            # t_2 comes from Step 3
            for is_after, t_2 in enumerate(direct_nodes):
                # x_1 is the edge we will break first
                # We need to determine which side we came from
                # If we don't we can expect to get no improvements
                x_1 = (t_2, t_1) if is_after == 0 else (t_1, t_2)
                cost_x_1 = dist_matrix[x_1[0]][x_1[1]]

                # now we find y_1
                distances_from_t_2 = dist_matrix[t_2]

                for t_3, cost_y_1 in enumerate(distances_from_t_2):
                    Y = set()
                    X = set([x_1])
                    # Don't go to the other neighbour of t_1
                    # Doing it by hand it broke everytime
                    if t_3 in direct_nodes or t_3 == t_2:
                        continue

                    y_1 = (t_2, t_3)

                    # X and Y must remain disjoint
                    if y_1 in X:
                        continue

                    g_1 = cost_x_1 - cost_y_1

                    if g_1 > 0:
                        Y.add(y_1)
                        #print(t_1, t_2, t_3)
                        # By this point we have our initial x_1, y_1, t_1, t_2 and t_3
                        # now we enter the recursive steps to try and improve the solution
                        # these ones are worth going further with
                        res = pick_x_i(nodes, T, t_1, t_3, X, Y, g_1)

                        if res != None:
                            #print(res)
                            temp_tour, temp_score = res

                            #print(temp_score)

                            if improved_score == None or improved_score > temp_score:
                                improved = True
                                improved_tour, improved_score = temp_tour, temp_score


        if improved:
            #print("Org:", T, f_T)
            #print("Better:", improved_tour, improved_score)
            T, f_T = improved_tour, improved_score

    print("ITS", its)

    return T, f_T


base = [2, 8, 6, 5, 10, 0, 9, 4, 11, 1, 3, 7]
#random.shuffle(base)

a = time.time()
tour, tour_length = lk(base)
b = time.time() - a

print(base, tour_length_calc(base))
print(convert_edge_list_to_tour(tour), tour_length_calc(convert_edge_list_to_tour(tour)), tour_length)

print("Took", b, "s")

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
