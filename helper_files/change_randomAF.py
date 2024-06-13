import random
import sys

angsd_freq = sys.argv[1]
random_freq =sys.argv[2]

# file= '/projects/korneliussen/people/vlr112/final_simulations/af_copy.txt'
def read_float_vector_file(file_path):
    with open(file_path, 'r') as file:
        vector = [float(line.strip()) for line in file]
    return vector

af_vector = read_float_vector_file(angsd_freq)

with open(random_freq, 'w') as f:
    for value in af_vector:
        change = random.uniform(-0.15, 0.10) * value
        randomized_value = value + change
        if randomized_value < 0:
            randomized_value == 0
        if randomized_value > 1:
            randomized_value == 1
        f.writelines(f"{randomized_value}")
        f.writelines("\n")