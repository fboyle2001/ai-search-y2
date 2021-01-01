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

import random
import math
from pprint import pprint

def tour_length_calc(state):
    dist = None

    try:
        dist = dist_matrix[state[0]][state[len(state) - 1]]
    except:
        print("State:", state)
        print("L:", len(state))
        sys.exit(1)

    for i in range(0, len(state) - 1):
        dist += dist_matrix[state[i]][state[i + 1]]

    return dist

seen_hashes = set()
stored_scores = dict()

def two_opt(start_solution):
    best_solution = start_solution[:]
    best_score = tour_length_calc(best_solution)
    improved = True
    hits = 0
    scorings = 0

    while improved:
        improved = False

        for i in range(0, num_cities):
            for k in range(i + 1, num_cities):
                new_solution = best_solution[:i] + best_solution[i:k + 1][::-1] + best_solution[k + 1:]
                new_score = 0

                h = hash(tuple(new_solution))
                scorings += 1

                if h in seen_hashes:
                    new_score = stored_scores[h]
                    hits += 1
                else:
                    new_score = tour_length_calc(new_solution)
                    seen_hashes.add(h)
                    stored_scores[h] = new_score

                if new_score < best_score:
                    improved = True
                    best_solution = new_solution
                    best_score = new_score

    print("Hits:", hits, "/", scorings, (hits/scorings) * 100, "%")
    return { "solution": best_solution, "score": best_score }

def alter_local_solution(current_solution):
    # Start by randomly swapping two cities
    altered = current_solution[:]

    first_i = random.randint(0, num_cities - 1)
    second_i = random.randint(0, num_cities - 1)

    altered[first_i], altered[second_i] = altered[second_i], altered[first_i]

    # Then apply 2-opt
    return two_opt(altered)

class WaterFlow:
    # Params from the original papers
    T = 20
    SPLIT_UPPER_LIMIT = 3
    g = 9.81
    EVAPORATION_RATE = 1 / 5

    def __init__(self, solution, w, v, score=None):
        self.solution = solution
        self.w = w
        self.v = v
        self.score = score if score != None else tour_length_calc(solution)
        self.comparable_hash = hash(tuple(solution))

    def get_momentum(self):
        return self.w * self.v

    def does_split(self):
        return not(self.is_regular_flow()) and not(self.is_stagnant())

    def is_stagnant(self):
        return self.get_momentum() == 0

    def is_regular_flow(self):
        return 0 < self.get_momentum() < WaterFlow.T

    def make_move(self):
        if self.does_split():
            return self.split_into_subflows()

        if self.is_stagnant():
            return [self]

        return self.move_to_new_location()

    def move_to_new_location(self):
        print("In move_to_new_location")
        altered = alter_local_solution(self.solution)
        improved_solution = altered["solution"]
        improved_score = altered["score"]
        score_improvement = self.score - improved_score
        square_vel = pow(self.v, 2) + 2 * WaterFlow.g * score_improvement
        new_vel = 0

        if square_vel > 0:
            new_vel = math.sqrt(square_vel)

        return [WaterFlow(improved_solution, self.w, new_vel)]

    def split_into_subflows(self):
        # This works by making small changes in the neighbourhood
        # Then calculate their change in obj score from original
        # Their velocity is given by sqrt(V_i^2 + 2*g*(obj score change)) if it is > 0 (see eq3)
        # Their mass is given by their relative rankings (see eq2)

        # The number of subflows in given by eq1
        number_of_subflows = int(min(max(1, self.get_momentum() // WaterFlow.T), WaterFlow.SPLIT_UPPER_LIMIT))

        if not isinstance(number_of_subflows, int):
            print("Non int??")
            print(number_of_subflows)
            print(int(number_of_subflows))
            sys.exit(1)

        improved_solutions = [alter_local_solution(self.solution) for x in range(number_of_subflows)]
        sorted_solutions = sorted(improved_solutions, key=lambda obj: obj["score"], reverse=False)

        # Now calculate the weights for the new solutions
        rank_total = sum(range(number_of_subflows + 1))
        new_flows = []

        for k, obj in enumerate(sorted_solutions):
            square_vel = pow(self.v, 2) + 2 * WaterFlow.g * (self.score - obj["score"])
            vel = 0

            if square_vel > 0:
                vel = math.sqrt(square_vel)

            new_flows.append(WaterFlow(obj["solution"], ((number_of_subflows + 1 - (k + 1)) / rank_total) * self.w, vel, score=obj["score"]))

        return new_flows

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
    active_flows = [WaterFlow(initial_solution, w_nought, v_nought, score=best_score)]

    for iteration in range(max_it):
        print("Iteration:", iteration)
        print("Best sol:", best_solution)
        print("Best score:", best_score)
        print("Total Active Flows:", len(active_flows))
        print()

        new_flow_dict = {}

        # Find the best solution from the current flows and
        # find the next flows including the subflows
        for flow in active_flows:
            if flow.score < best_score:
                best_score = flow.score
                best_solution = flow.solution

            new_flow_arr = flow.make_move()

            for nf in new_flow_arr:
                if nf.comparable_hash in new_flow_dict.keys():
                    new_flow_dict[nf.comparable_hash]["quantity"] += 1
                else:
                    new_flow_dict[nf.comparable_hash] = {
                        "flow": nf,
                        "quantity": 1
                    }

        new_flows = []

        # Merge the flows
        for key in new_flow_dict:
            obj = new_flow_dict[key]
            if obj["quantity"] == 1:
                new_flows.append(obj["flow"])
            else:
                flow = obj["flow"]
                combined_weight = flow.w
                combined_velocity = flow.v

                for _ in range(1, obj["quantity"]):
                    combined_weight += flow.w
                    # Note that eq5 uses W_i + W_j but combined_weight is already that
                    combined_velocity = (combined_weight * combined_velocity + flow.w * flow.v) / (combined_weight)

                new_flows.append(flow)

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
        # ** MAY NEED TO REVISIT SEE EQ8 **
        if(velocity_sum < 0.01):
            for flow in new_flows:
                flow.w = (flow.w / mass_sum) * w_nought
                flow.v = v_nought

        # Regular precipitation
        for flow in new_flows:
            flow.w = (flow.w / mass_sum) * w_nought - mass_sum

        active_flows = new_flows

    # Finally retrieve the best from the final iteration
    for flow in active_flows:
        if flow.score < best_score:
            best_score = flow.score
            best_solution = flow.solution

    return best_solution, best_score

initial_solution = [220, 521, 97, 131, 403, 410, 204, 499, 207, 217, 5, 294, 301, 430, 235, 309, 332, 413, 10, 257, 282, 65, 3, 339, 9, 492, 359, 42, 104, 440, 212, 69, 72, 200, 457, 498, 515, 476, 78, 125, 458, 73, 512, 363, 87, 24, 378, 157, 234, 213, 265, 534, 305, 98, 321, 297, 351, 358, 163, 325, 241, 340, 51, 356, 120, 348, 166, 269, 424, 336, 46, 165, 95, 88, 32, 429, 431, 448, 134, 158, 194, 379, 530, 449, 176, 441, 243, 21, 118, 361, 156, 447, 231, 380, 0, 67, 376, 227, 382, 407, 102, 354, 38, 327, 360, 510, 489, 258, 54, 478, 278, 421, 151, 398, 189, 195, 143, 178, 43, 443, 455, 337, 259, 460, 14, 428, 94, 240, 395, 292, 188, 190, 250, 185, 419, 412, 468, 288, 58, 145, 383, 394, 36, 475, 179, 40, 437, 409, 517, 81, 329, 426, 254, 526, 70, 436, 106, 59, 318, 525, 255, 514, 52, 101, 2, 390, 245, 304, 293, 153, 408, 459, 7, 502, 528, 507, 400, 298, 100, 119, 272, 275, 276, 82, 296, 404, 37, 501, 286, 141, 152, 35, 84, 444, 183, 30, 113, 267, 22, 450, 364, 283, 186, 191, 225, 238, 25, 150, 374, 438, 393, 268, 389, 174, 420, 375, 132, 485, 114, 367, 175, 6, 77, 148, 161, 427, 506, 232, 423, 312, 405, 343, 371, 197, 532, 402, 322, 162, 181, 471, 17, 320, 456, 111, 260, 347, 505, 491, 529, 352, 500, 1, 117, 467, 509, 334, 344, 103, 435, 33, 274, 414, 422, 55, 19, 126, 50, 472, 99, 406, 392, 279, 503, 520, 316, 110, 328, 388, 47, 252, 149, 142, 196, 396, 261, 154, 246, 342, 170, 115, 432, 108, 211, 357, 333, 262, 139, 366, 159, 511, 18, 480, 192, 133, 62, 284, 300, 13, 522, 496, 28, 417, 249, 416, 533, 56, 239, 128, 8, 26, 222, 397, 487, 137, 138, 433, 488, 373, 29, 135, 483, 484, 53, 20, 497, 123, 387, 187, 229, 214, 461, 228, 216, 86, 495, 41, 355, 313, 266, 130, 68, 271, 122, 140, 264, 504, 182, 109, 253, 147, 465, 518, 263, 4, 353, 307, 236, 173, 331, 247, 202, 401, 399, 287, 486, 215, 218, 112, 474, 146, 345, 469, 27, 311, 49, 451, 323, 66, 425, 45, 513, 386, 219, 303, 324, 96, 16, 92, 418, 233, 116, 237, 290, 201, 168, 167, 338, 302, 76, 481, 350, 446, 93, 462, 85, 107, 519, 75, 121, 494, 277, 164, 439, 477, 370, 90, 74, 464, 63, 346, 299, 490, 209, 223, 124, 11, 89, 206, 248, 527, 60, 369, 15, 289, 335, 180, 372, 224, 415, 326, 244, 34, 144, 39, 91, 445, 61, 341, 452, 48, 208, 210, 280, 524, 127, 198, 310, 71, 129, 64, 508, 454, 184, 230, 330, 411, 319, 31, 493, 251, 171, 12, 368, 57, 381, 453, 169, 44, 79, 160, 155, 470, 315, 466, 83, 242, 23, 270, 205, 136, 391, 177, 193, 349, 516, 377, 172, 256, 531, 482, 473, 221, 226, 105, 384, 314, 442, 281, 463, 199, 434, 273, 285, 479, 308, 295, 385, 306, 291, 523, 365, 317, 80, 362, 203]
print(len(initial_solution))
max_it = 6

s = time.time()
tour, tour_length = water_flow_optimise(initial_solution, max_it, 8, 5)
p = time.time() - s

print(tour)
print("Score:", tour_length)
print("Took", p, "s")













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
