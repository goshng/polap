################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

# polaplib/polap-py-unique-mtcontigs.py

################################################################################
# This script is very complicated script, which used to be written in
# bash shell script.
#
# command-line options:
# '--input', 'input/run-polap-py-unique-mtcontigs.input1.2.tmp',
# '--out', 'output/run-polap-py-unique-mtcontigs.output1.2.tmp',
# '--prefix', 'output/run-polap-py-unique-mtcontigs.output1.2/mt.contig.name'
#
# Files
# =====
#
# input1:
# input/run-polap-py-unique-mtcontigs.input1.2.tmp
# ------------------------------------------------
# 2 input/2/8-mt.contig.name.txt
#
# file: input/2/8-mt.contig.name.txt
# ----------------------------------
# edge_1
#
# output/run-polap-py-unique-mtcontigs.output1.2.tmp
# -------------------------------------------------- 
# output/run-polap-py-unique-mtcontigs.output1.2/mt.contig.name-1.txt 2
# 
# prefix: output/run-polap-py-unique-mtcontigs.output1.2/mt.contig.name
#
# actual output file:
# output/run-polap-py-unique-mtcontigs.output1.2/mt.contig.name-1.txt
# -------------------------------------------------------------------
# edge_1
#
# Example:
# python polap-py-unique-mtcontigs.py 
#   --input input/run-polap-py-unique-mtcontigs.input1.2.tmp 
#   --out output/run-polap-py-unique-mtcontigs.output1.2.tmp 
#   --prefix output/mt.contig.name
#
# This script is used to finalize the mt.contig.name after selecting seeds.
#
# Used by:
# polap_disassemble-seeds() {
#
# Check:
################################################################################

# %%
import os
import argparse

# %%
debug = os.getenv("_POLAP_DEBUG", "0")

# %%
# Read mt.contig.name file to sort the edge number strings.
# We call this function for each of the list of mt.contig.name files.
def read_file_to_set(filepath):
    with open(filepath, "r") as f:
        return sorted(set(line.strip() for line in f if line.strip()))

# %%
def main():

# %%
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Process mtcontigs files to find unique sets."
    )

    # Parse command-line arguments
    parser.add_argument(
        "-i",
        "--input",
        dest="input_file",
        required=True,
        help="Path to mtcontig_files.txt",
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Output index file (e.g., mtcontig-set.txt)"
    )
    parser.add_argument(
        "-p",
        "--prefix",
        required=True,
        help="Prefix for the new set files (e.g., mtcontig)",
    )

#%%
    # Test arguments
    # custom_args = ['--input', 'input/run-polap-py-unique-mtcontigs.input1.2.tmp',
    #     '--out', 'output/run-polap-py-unique-mtcontigs.output1.2.tmp',
    #     '--prefix', 'output/run-polap-py-unique-mtcontigs.output1.2/mt.contig.name']
    # args = parser.parse_args(custom_args)
    # print(args)

# %%
    # Actual arguments from the command-line
    args = parser.parse_args()

# %%
    mtcontigs_file = args.input_file
    output_index_file = args.out
    prefix = args.prefix

# %%
    # An empty file path dictionary for index and its file path
    file_paths = {}

# %%
    # file_paths contains these two columns: index number and path string
    # The index is the plastid seed contig selection scheme number.
    # file_paths has all the mt.contig.name candidate files.
    # We could have multiple seed selection schemes for mtDNA case.
    # For ptDNA case, we have only one, so we might not need this, though.
    # We want to find a unique set of such mt.contig.name files based on
    # the contents of the mt.contig.name files.
    # Example: mtcontigs_file
    # 0 /path/to/file0.txt
    # 1 /path/to/file1.txt
    with open(mtcontigs_file, "r") as f:
        for line in f:
            index, path = line.strip().split()
            file_paths[int(index)] = path

# %%
    # mt.contig.name file is read to create a set of mt.contig.name content
    # We want to find a unique set of mt.contig.name files based on the 
    # content.
    indexed_sets = {index: read_file_to_set(path) for index, path in file_paths.items()}

# %%
    # The indexed_sets has the contig selection scheme index and its
    # seed contigs. Two schemes can generate two seemingly different
    # sets of the seed contigs if they are not properly ordered.
    # Tuple can be used as an index to find unique items, we
    # use set_tuple for that.
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

# %%
    # Write out each mt.contig.name file with its edge numbers.
    # Edge numbers are printed one per line using:
    # "\n".join(unique_set)
    for filename, unique_set in unique_sets.items():
        with open(filename, "w") as f:
            f.write("\n".join(unique_set) + "\n")

# %%
    # Now, after having all mt.contig.name files, we record the new
    # index for each contig seed selection scheme, so that we can
    # find which mt.contig.name corresponds to which selection scheme.
    with open(output_index_file, "w") as f:
        for unique_set_tuple, unique_index in set_to_index.items():
            original_indices = [
                str(idx)
                for idx, s in indexed_sets.items()
                if tuple(s) == unique_set_tuple
            ]
            f.write(f"{prefix}-{unique_index}.txt {' '.join(original_indices)}\n")

# %%

if __name__ == "__main__":
    main()
