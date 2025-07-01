#!/usr/bin/env bash
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

################################################################################
# Create a bash script for building (or installing) a conda package.
# 1. A bioconda package needs a file named build.sh. It copies files
# from a conda package to the destination folder. This script creates build.sh
# that is used by a bioconda package recipes.
# 2. We also use the build.sh file to update the conda environment with a local
# source file or github repository.
#
# FIXME:
# We may not need this script much because we used this script before we have
# polaplib folder.
#
# Check: 2025-06-16
################################################################################

_head_message=$(
	cat <<HEREDOC
#!/usr/bin/env bash

mkdir -p \$PREFIX/bin

files=(polap
	polap-ncbitools
  polap-batch-v0.sh
  polap-batch-v2.sh
  polap-data-v0.sh
  polap-data-v1.sh
  polap-data-v2.sh
  polap-data-v3.sh
  polap-data-v4.sh
	polap-data-aflye
	polap-data-cflye
	polap-data-dflye
	polap-data-taxon
  bolap
  bolap.sh
HEREDOC
)

_tail_message=$(
	cat <<HEREDOC
	polap.sh)

for i in "\${files[@]}"; do
	cp src/\$i \$PREFIX/bin
  chmod +x \$PREFIX/bin/\$i
done
cp -pr src/polaplib \$PREFIX/bin

chmod +x \$PREFIX/bin/polap
chmod +x \$PREFIX/bin/polap.sh
chmod +x \$PREFIX/bin/polap-ncbitools
chmod +x \$PREFIX/bin/polap-data-aflye
chmod +x \$PREFIX/bin/polap-data-cflye
chmod +x \$PREFIX/bin/polap-data-dflye
chmod +x \$PREFIX/bin/polap-data-taxon
chmod +x \$PREFIX/bin/bolap
chmod +x \$PREFIX/bin/bolap.sh
HEREDOC
)

echo "${_head_message}"

ls -1 *.sh |
	grep -v '^build\.sh$' |
	grep -v '^polap-build\.sh$' |
	grep -v '^detect_plastome_structure\.sh$' |
	grep -v '^ptGAUL2\.sh$' |
	grep -v '^1\.sh$' |
	grep -v '^2\.sh$' |
	grep -v '^x\.sh$' |
	grep -v '^polap\.sh$' |
	sed 's/^/\t/'
# ls -1 *.py | sed 's/^/\t/'
# ls -1 *.R | sed 's/^/\t/'
# ls -1 *.csv | sed 's/^/\t/'
# ls -1 polap-template-*.txt | sed 's/^/\t/'

echo "${_tail_message}"

exit

files=(polap
	polap-conda-environment-fmlrc.yaml
	polap-conda-environment.yaml
	polap-mt.1.c70.3.faa
	polap-pt.2.c70.3.faa
	run-polap-ncbitools
	polap-command-completion.sh
	polap-constants.sh
	polap-function-set-variables.sh
	polap-git-hash-version.sh
	polap-package-common.sh
	polap-package-mtcontigs.sh
	polap-parsing.sh
	polap-test-variables.sh
	polap-variables-common.sh
	polap-variables-main.sh
	polap-variables-mtcontigs.sh
	polap-version.sh
	run-polap-function-annotate-contig.sh
	run-polap-function-annotate.sh
	run-polap-function-archive.sh
	run-polap-function-assemble.sh
	run-polap-function-bioproject.sh
	run-polap-function-errors.sh
	run-polap-function-include.sh
	run-polap-function-log.sh
	run-polap-function-menus.sh
	run-polap-function-miscellaneous.sh
	run-polap-function-mtdna.sh
	run-polap-function-oga.sh
	run-polap-function-polishing.sh
	run-polap-function-seeds.sh
	run-polap-function-template.sh
	run-polap-function-test.sh
	run-polap-function-utilities.sh
	run-polap-function-wga.sh
	polap-bash-create-depth-file.sh
	polap-bash-half-cut.sh
	polap-bash-minimap2-paf2tab.sh
	polap-py-find-cc.py
	run-polap-py-mtdna-find-cycle-with-node-revisits.py
	run-polap-py-select-mtdna-2-nx-find-circular-path.py
	run-polap-py-select-mtdna-2-nx-simple-cycles.py
	run-polap-pairs.R
	run-polap-r-assemble-bioproject-3-length-match.R
	run-polap-r-blast-mtdna-1-determine-gene.R
	run-polap-r-bridge.R
	run-polap-r-cc2mtcontig.R
	run-polap-r-contig2edge.R
	run-polap-r-depth-distribution.R
	run-polap-r-depthfilter-gfa.R
	run-polap-r-determine-depth-range.R
	run-polap-r-determine-depth-range_1.R
	run-polap-r-determine-depth-range_2.R
	run-polap-r-determine-depth-range_3.R
	run-polap-r-determine-depth-range_4.R
	run-polap-r-determine-depth-range_5.R
	run-polap-r-determine-depth-range_6.R
	run-polap-r-edges-stats.R
	run-polap-r-final-filter-mtcontig.R
	run-polap-r-final-mtcontig.R
	run-polap-r-final-seed-mtcontig.R
	run-polap-r-genes-bed4.R
	run-polap-r-genes.R
	run-polap-r-get-bioproject-1.R
	run-polap-r-jellyfish.R
	run-polap-r-mtcontig-contig.R
	run-polap-r-mtcontig.R
	run-polap-r-pairs.R
	run-polap-r-plastid-determine-depth-range.R
	run-polap-r-plastid-determine-depth-range_1.R
	run-polap-r-plastid-determine-depth-range_6.R
	run-polap-r-plot-mtdna.R
	run-polap-r-prepare-cc.R
	run-polap-r-preselect-annotation.R
	run-polap-r-select-mtdna-1-nx-gfa-links.R
	run-polap-r-template.R
	run-polap-r-test-reads-bar-graph.R
	polap.sh)

for i in "${files[@]}"; do
	cp src/$i $PREFIX/bin
done

chmod +x $PREFIX/bin/polap
chmod +x $PREFIX/bin/polap.sh
