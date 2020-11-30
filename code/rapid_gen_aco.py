import aco
import signal
import time

active = True

def sig_intercept(sig, frame):
    global active
    active = False
    print("Finishing iteration...")

signal.signal(signal.SIGINT, sig_intercept)

nn_tour = aco.nearest_neighbour(0)
nn_tour_length = aco.tour_length_calc(nn_tour)
ants = aco.num_cities
max_it = 100
tau_nought = ants / nn_tour_length
alpha = 1
beta = 3
rho = 0.5

best_tour = []
best_tour_length = 20000000000
total_time = 0
total_runs = 0

worst_tour = []
worst_tour_length = 0

total_tour_length = 0

#while active:
time_start = time.time()
tour, tour_length = aco.ant_colony_optimise(tau_nought, nn_tour, ants, max_it, alpha, beta, rho)
time_end = time.time()

if tour_length < best_tour_length:
    best_tour = tour
    best_tour_length = tour_length

if tour_length > worst_tour_length:
    worst_tour = tour
    worst_tour_length = tour_length

time_taken = time_end - time_start
total_time += time_taken
total_runs += 1
total_tour_length += tour_length

print(total_runs, time_taken)

print()
print("Best")
print("Best Tour:", best_tour)
print("Best Tour Length:", best_tour_length)
print()

print("Worst")
print("Worst Tour:", worst_tour)
print("Worst Tour Length:", worst_tour_length)
print()

print("Average")
print("Average Tour Length:", total_tour_length / total_runs)
print()

print("Timing")
print("Total Time:", total_time)
print("Total Runs:", total_runs)
print("Avg Time:", total_time / total_runs)
print()
