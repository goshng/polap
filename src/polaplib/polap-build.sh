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
	polap-cmd-annotate.sh
	polap-cmd-archive.sh
	polap-cmd-assemble.sh
	run-polap-function-bioproject.sh
	run-polap-function-errors.sh
	run-polap-function-include.sh
	run-polap-function-log.sh
	run-polap-function-menus.sh
	run-polap-function-miscellaneous.sh
	polap-cmd-mtdna.sh
	polap-cmd-oga.sh
	polap-cmd-polishing.sh
	polap-cmd-seeds.sh
	run-polap-function-template.sh
	run-polap-function-test.sh
	run-polap-function-utilities.sh
	polap-cmd-wga.sh
	polap-bash-create-depth-file.sh
	polap-bash-half-cut.sh
	polap-bash-minimap2-paf2tab.sh
	polap-py-find-cc.py
	run-polap-py-mtdna-find-cycle-with-node-revisits.py
	run-polap-py-select-mtdna-2-nx-find-circular-path.py
	run-polap-py-select-mtdna-2-nx-simple-cycles.py
	run-polap-pairs.R
	polap-r-assemble-bioproject-3-length-match.R
	run-polap-r-blast-mtdna-1-determine-gene.R
	polap-r-bridge.R
	polap-r-cc2mtcontig.R
	run-polap-r-contig2edge.R
	polap-r-depth-distribution.R
	polap-r-depthfilter-gfa.R
	polap-r-determine-depth-range.R
	polap-r-determine-depth-range_1.R
	polap-r-determine-depth-range_2.R
	polap-r-determine-depth-range_3.R
	polap-r-determine-depth-range_4.R
	polap-r-determine-depth-range_5.R
	polap-r-determine-depth-range_6.R
	polap-r-edges-stats.R
	polap-r-final-filter-mtcontig.R
	polap-r-final-mtcontig.R
	polap-r-final-seed-mtcontig.R
	run-polap-r-genes-bed4.R
	polap-r-genes.R
	polap-r-get-bioproject.R
	polap-r-jellyfish.R
	polap-r-mtcontig-contig.R
	polap-r-mtcontig.R
	polap-r-pairs.R
	polap-r-plastid-determine-depth-range.R
	polap-r-plastid-determine-depth-range_1.R
	polap-r-plastid-determine-depth-range_6.R
	run-polap-r-plot-mtdna.R
	polap-r-prepare-cc.R
	polap-r-preselect-annotation.R
	run-polap-r-select-mtdna-1-nx-gfa-links.R
	run-polap-r-template.R
	polap-r-test-reads-bar-graph.R
	polap.sh)

for i in "${files[@]}"; do
	cp src/$i $PREFIX/bin
done

chmod +x $PREFIX/bin/polap
chmod +x $PREFIX/bin/polap.sh
