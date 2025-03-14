import os
import argparse


def read_file_to_set(filepath):
    """Read a file and return a sorted set of its lines."""
    with open(filepath, "r") as f:
        return sorted(line.strip() for line in f if line.strip())


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Process mtcontigs files to find unique sets."
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_file",
        required=True,
        help="Path to mtcontig_files.txt",
    )
    parser.add_argument(
        "-o", "--out", required=True, help="Output index file (e.g., mtcontig-set.txt)"
    )
    parser.add_argument(
        "-p",
        "--prefix",
        required=True,
        help="Prefix for the new set files (e.g., mtcontig)",
    )
    args = parser.parse_args()

    mtcontigs_file = args.input_file
    output_index_file = args.out
    prefix = args.prefix
    # output_dir = os.path.dirname(output_index_file)
    # if output_dir:
    #     os.makedirs(output_dir, exist_ok=True)

    # Step 1: Read mtcontigs_files.txt and load file paths with indices
    file_paths = {}
    with open(mtcontigs_file, "r") as f:
        for line in f:
            index, path = line.strip().split()
            file_paths[int(index)] = path

    # Step 2: Load each file as a sorted set and index by the line index
    indexed_sets = {index: read_file_to_set(path) for index, path in file_paths.items()}

    # Step 3: Find unique sets
    unique_sets = {}
    set_to_index = {}
    counter = 1

    for index, current_set in indexed_sets.items():
        current_set_tuple = tuple(current_set)  # Convert list to tuple for hashing
        if current_set_tuple not in set_to_index:
            unique_filename = f"{prefix}-{counter}.txt"
            unique_sets[unique_filename] = current_set
            set_to_index[current_set_tuple] = counter
            counter += 1

    # Step 4: Write unique sets to files
    for filename, unique_set in unique_sets.items():
        with open(filename, "w") as f:
            f.write("\n".join(unique_set) + "\n")

    # Step 5: Create the output index file
    with open(output_index_file, "w") as f:
        for unique_set_tuple, unique_index in set_to_index.items():
            original_indices = [
                str(idx)
                for idx, s in indexed_sets.items()
                if tuple(s) == unique_set_tuple
            ]
            f.write(f"{prefix}-{unique_index}.txt {' '.join(original_indices)}\n")


if __name__ == "__main__":
    main()
