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

input_file = "AISearchfile021.txt"

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

algorithm_code = "WF"

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

import math

# Calculates the length of a tour
def tour_length_calc(state):
    dist = dist_matrix[state[len(state) - 1]][state[0]]

    for i in range(0, len(state) - 1):
        dist += dist_matrix[state[i]][state[i + 1]]

    return dist

# Lin-K-H heuristic instead of 2-opt
def lkh(start_solution):
    # Copy the tour and find its length
    best_solution = start_solution[:]
    best_score = tour_length_calc(best_solution)

    # We stop when this is false
    improved = True

    while improved:
        improved = False

    return { "solution": best_solution, "score": best_score }

def alter_local_solution(current_solution):
    altered = current_solution[:]

    first_i = random.randint(0, num_cities - 1)
    second_i = random.randint(0, num_cities - 1)

    # Start by randomly swapping two cities
    altered[first_i], altered[second_i] = altered[second_i], altered[first_i]

    # Then apply 2-opt
    return lkh(altered)

# This represents a single flow in the algorithm
# References:
# (1) Ayman Srour, Zulaiha Ali Othman, Abdul Razak Hamdan, "A Water Flow-Like Algorithm for the Travelling Salesman Problem", Advances in Computer Engineering, vol. 2014, Article ID 436312, 14 pages, 2014. https://doi.org/10.1155/2014/436312 (https://www.hindawi.com/journals/aceng/2014/436312/)
# (2) Feng-Cheng Yang & Yuan-Peng Wang (2007) WATER FLOW-LIKE ALGORITHM FOR OBJECT GROUPING PROBLEMS, Journal of the Chinese Institute of Industrial Engineers, 24:6, 475-488, DOI: 10.1080/10170660709509062 (http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.120.117&rep=rep1&type=pdf)
class WaterFlow:
    # Parameters from the original papers
    T = 20
    SPLIT_UPPER_LIMIT = 3
    g = 9.81
    EVAPORATION_RATE = 1 / 5

    # w, v are weight and velocity
    def __init__(self, solution, w, v, score=None):
        self.solution = solution
        self.w = w
        self.v = v
        self.score = score if score != None else tour_length_calc(solution)
        self.comparable_hash = hash(tuple(solution))

    # Used to determine if a flow will split
    def get_momentum(self):
        return self.w * self.v

    def does_split(self):
        return not(self.is_regular_flow()) and not(self.is_stagnant())

    # Stagnant solutions represent a local minimum
    def is_stagnant(self):
        return self.get_momentum() == 0

    # Flows in this range aren't going to split but they also aren't stagnant
    def is_regular_flow(self):
        return 0 < self.get_momentum() < WaterFlow.T

    # Determines the next step for this flow
    # Either splits, returns itself (stagnant) or flows to a new location
    def make_move(self):
        if self.does_split():
            return self.split_into_subflows()

        if self.is_stagnant():
            return [self]

        return self.move_to_new_location()

    def move_to_new_location(self):
        # Alter the solution
        altered = alter_local_solution(self.solution)
        improved_solution = altered["solution"]
        improved_score = altered["score"]

        # Calculate its velocity
        score_improvement = self.score - improved_score
        square_vel = pow(self.v, 2) + 2 * WaterFlow.g * score_improvement
        new_vel = 0

        # A negative square_vel implies that this solution isn't worth pursuing
        if square_vel > 0:
            new_vel = math.sqrt(square_vel)

        # Return in an array so it is the same format as subflows
        return [WaterFlow(improved_solution, self.w, new_vel)]

    def split_into_subflows(self):
        # This works by making small changes in the neighbourhood
        # Then calculate their change in obj score from original
        # eqX references paper (1)
        # Their velocity is given by sqrt(V_i^2 + 2*g*(obj score change)) if it is > 0 (see eq3)
        # Their mass is given by their relative rankings (see eq2)

        # The number of subflows in given by eq1
        number_of_subflows = int(min(max(1, self.get_momentum() // WaterFlow.T), WaterFlow.SPLIT_UPPER_LIMIT))

        # Get the altered solutions and sort by the score so the best is at index 0
        improved_solutions = [alter_local_solution(self.solution) for x in range(number_of_subflows)]
        sorted_solutions = sorted(improved_solutions, key=lambda obj: obj["score"], reverse=False)

        # Now calculate the weights for the new solutions
        rank_total = sum(range(number_of_subflows + 1))
        new_flows = []

        for k, obj in enumerate(sorted_solutions):
            # Same principle as move_to_new_location now
            square_vel = pow(self.v, 2) + 2 * WaterFlow.g * (self.score - obj["score"])
            vel = 0

            if square_vel > 0:
                vel = math.sqrt(square_vel)

            new_flows.append(WaterFlow(obj["solution"], ((number_of_subflows + 1 - (k + 1)) / rank_total) * self.w, vel, score=obj["score"]))

        return new_flows

# Basic greedy implementation of nearest_neighbour to give a base solution to the TSP for this city set
def nearest_neighbour(source):
    tour = []
    current = source

    while True:
        neighbours = dist_matrix[current]
        nearest = None
        nearest_dist = max(neighbours) + 1

        for i, dist in enumerate(neighbours):
            if i == current:
                continue

            if i in tour:
                continue

            if nearest_dist > dist:
                nearest_dist = dist
                nearest = i

        tour.append(current)

        if nearest == None:
            break

        current = nearest

    return tour

def water_flow_optimise(initial_solution, max_it, w_nought, v_nought):
    best_solution = initial_solution
    best_score = tour_length_calc(best_solution)
    # Keep track of the flows that are active in this iteration
    active_flows = [WaterFlow(initial_solution, w_nought, v_nought, score=best_score)]

    for iteration in range(max_it):
        # Debugging information
        # print("Iteration:", iteration)
        # print("Best sol:", best_solution)
        # print("Best score:", best_score)
        # print("Total Active Flows:", len(active_flows))
        # print()

        # Use a dict so we can merge the flows easily
        new_flow_dict = {}

        # Find the best solution from the current flows and
        # find the next flows including the subflows
        for flow in active_flows:
            if flow.score < best_score:
                best_score = flow.score
                best_solution = flow.solution

            new_flow_arr = flow.make_move()

            # We use the hash of the new flows to begin the merging process
            for nf in new_flow_arr:
                if nf.comparable_hash in new_flow_dict.keys():
                    new_flow_dict[nf.comparable_hash]["quantity"] += 1
                else:
                    new_flow_dict[nf.comparable_hash] = {
                        "flow": nf,
                        "quantity": 1
                    }

        # We'll squash the dictionary into this array
        new_flows = []

        # Merge the flows
        for key in new_flow_dict:
            obj = new_flow_dict[key]
            if obj["quantity"] == 1:
                new_flows.append(obj["flow"])
            else:
                # Merge the flows if there are multiple with the same hash (and thus the same solution)
                flow = obj["flow"]
                combined_weight = flow.w
                combined_velocity = flow.v

                # Probably doable without the loop but obj["quantity"] <= WaterFlow.SPLIT_UPPER_LIMIT (= 3)
                # so it is insignificant
                for _ in range(1, obj["quantity"]):
                    combined_weight += flow.w
                    # Note that eq5 uses W_i + W_j but combined_weight is already that
                    combined_velocity = (combined_weight * combined_velocity + flow.w * flow.v) / (combined_weight)

                new_flows.append(flow)

        # We'll use these to do precipitation after the evaporation
        velocity_sum = 0
        mass_sum = 0

        # Now apply evaporation
        for flow in new_flows:
            flow.w *= (1 - WaterFlow.EVAPORATION_RATE)
            velocity_sum += flow.v
            mass_sum += flow.w

        # Enforced precipitation
        # Occurs when all velocities are 0
        # Thinking about floats here, we may have rounding issues so go for below boundary instead
        # Prevents the whole system from stagnating
        if(velocity_sum < 0.01):
            for flow in new_flows:
                flow.w = (flow.w / mass_sum) * w_nought
                flow.v = v_nought

        # Regular precipitation
        # Redistribute the evaporated water
        for flow in new_flows:
            flow.w = (flow.w / mass_sum) * w_nought - mass_sum

        # Set the active flows for the next iteration
        active_flows = new_flows

    # Finally retrieve the best from the final iteration
    # Otherwise we did all that work for nothing
    for flow in active_flows:
        if flow.score < best_score:
            best_score = flow.score
            best_solution = flow.solution

    return best_solution, best_score

initial_solution = nearest_neighbour(0)
max_it = 2
w_nought = 8
v_nought = 5
# Additional Parameters are set in Waterflow:
# T = 20
# SPLIT_UPPER_LIMIT = 3
# g = 9.81
# EVAPORATION_RATE = 1 / 5

tour, tour_length = water_flow_optimise(initial_solution, max_it, w_nought, v_nought)

added_note = f"max_it = {max_it}, w_nought = {w_nought}, v_nought = {v_nought}, T = {WaterFlow.T}, SPLIT_UPPER_LIMIT = {WaterFlow.SPLIT_UPPER_LIMIT}, g = {WaterFlow.g}, EVAPORATION_RATE = {WaterFlow.EVAPORATION_RATE}"

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
