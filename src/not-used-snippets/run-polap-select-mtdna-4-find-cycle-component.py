import sys

from postman_problems.solver import cpp
from postman_problems.stats import calculate_postman_solution_stats


def write_array_of_arrays_to_file(array_of_arrays, file_name):
    with open(file_name, "w") as f:
        for pair in array_of_arrays:
            # Join the two strings with a comma or any delimiter, and write each pair to a new line
            f.write(f"{pair[0]},{pair[1]}\n")


# Example usage:
# array_of_arrays = [["apple", "banana"], ["cat", "dog"], ["red", "blue"]]
# file_name = "output_file.txt"
# write_array_of_arrays_to_file(array_of_arrays, file_name)

if __name__ == "__main__":

    # Create a graph given in the above diagram
    # 5 vertices numbered from 0 to 4

    # g = Graph(5)
    if len(sys.argv) != 3:
        print("Usage: python script.py <path_to_file> <outfile>")
        sys.exit(1)

    file_path = sys.argv[1]
    file_path_out = sys.argv[2]

    # find CPP solution
    circuit, graph = cpp(
        edgelist_filename=file_path,
    )

    # print solution route
    # for e in circuit:
    #     print(e[0], ",", e[1])

    write_array_of_arrays_to_file(circuit, file_path_out)

    # print solution summary stats
    # for k, v in calculate_postman_solution_stats(circuit).items():
    # print(k, v)
