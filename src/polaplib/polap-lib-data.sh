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
# Ensure that the current script is sourced only once
source "${_POLAPLIB_DIR}/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u
if [[ -n "${!_POLAP_INCLUDE_}" ]]; then
	set -u
	return 0
fi
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

function _polap_lib_data-example-data {
	local _brg_data="${1}"

	local _data=$(
		cat <<HEREDOC
species,long,short
Anthoceros_agrestis,SRR10190639,SRR10250248
Arabidopsis_thaliana,ERR2173373,ERR2173372
Canavalia_ensiformis,SRR18714551,SRR18714547
Cinchona_pubescens,SRR20784020,SRR20784021
Codonopsis_lanceolata,SRR11585869,SRR11585868
Cucumis_sativus_var_hardwickii,SRR28091980,SRR28091977
Dioscorea_japonica,SRR16108312,SRR16108386
Dunaliella_tertiolecta,SRR22857204,SRR22857205
Eucalyptus_pauciflora,SRR7153095,SRR7161123
Euonymus_alatus,SRR16133411,SRR16122871
Gossypium_herbaceum,SRR17224483,SRR17211914
Juncus_effusus,SRR14298760,SRR14298746
Juncus_inflexus,SRR14298751,SRR14298745
Juncus_roemerianus,SRR21976090,SRR21976092
Juncus_validus,SRR21976089,SRR21976091
Leiosporoceros_dussii,SRR25387688,SRR25387689
Macadamia_jansenii,SRR11191910,SRR11191912
Musa_acuminata_subsp_malaccensis,ERR5455028,ERR3606950
Notothylas_orbicularis,SRR25405055,SRR25405056
Ophrys_lutea,ERR5167480,ERR5303530
Oryza_rufipogon,SRR12104676,SRR12102351
Phaeomegaceros_chiloensis,SRR25430413,SRR25430414
Populus_x_sibirica,SRR15146668,SRR12963707
Prunus_mandshurica,ERR4656977,ERR4656978
Solanum_lycopersicum,SRR11073833,SRR11205474
Spirodela_polyrhiza,SRR11472010,SRR11472009
Vaccinium_vitis-idaea,SRR25468450,SRR25477290
Vitis_vinifera,SRR26163227,SRR26163231
HEREDOC
	)

	echo "${_data}" >"${_brg_data}"
	echo "create ${_brg_data}"
}

help_message_bleeding_edge_polap=$(
	cat <<HEREDOC

  Patch polap with the latest github version.
HEREDOC
)

help_message_patch_polap=$(
	cat <<HEREDOC

  Patch polap with a specific version.
  _polap_version=<VERSION> polap patch-polap
HEREDOC
)

help_message_install_getorganelle=$(
	cat <<HEREDOC

  Install conda environments: getorganelle
  conda create -y --name getorganelle bioconda::getorganelle
  conda activate getorganelle
  get_organelle_config.py --add embplant_pt,embplant_mt
HEREDOC
)

help_message_install_cflye=$(
	cat <<HEREDOC

  Install cflye to polap conda environment
  conda install -y goshng::cflye
HEREDOC
)

help_message_install_fmlrc=$(
	cat <<HEREDOC

  Install conda environments: polap-fmlrc
HEREDOC
)

help_message_install_polap=$(
	cat <<HEREDOC

  Install conda environments: polap
  conda create -y --name polap bioconda::polap
HEREDOC
)

help_message_uninstall=$(
	cat <<HEREDOC

  Uninstall conda environments: polap, polap-fmlrc, getorganelle
  conda remove -n polap --all
  conda remove -n polap-fmlrc --all
  conda remove -n getorganelle --all
  conda env remove -n pmat
HEREDOC
)

function _polap_lib_data-execute-common-subcommand {
	local subcmd1="$1"
	local opt_y_flag="$2"
	local handled=0

	case "$subcmd1" in
	install-conda | setup-conda | uninstall | \
		install-polap | download-polap-github | delete-polap-github | \
		install-fmlrc | patch-polap | bleeding-edge-polap | local-edge-polap | \
		test-polap | \
		install-bandage | install-getorganelle | install-cflye | install-dflye | \
		install-tippo | \
		install-pmat | download-pmat)
		handled=1
		;;
	mkdir-all | rm-empty | rm)
		handled=1
		;;
	esac

	case "$subcmd1" in
	install-conda)
		read -p "Do you want to install miniconda3? (y/N): " confirm
		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			if [[ -d "$HOME/miniconda3" ]]; then
				echo "ERROR: you already have miniconda3: $HOME/miniconda3"
			else
				mkdir -p ~/miniconda3
				wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
					-O ~/miniconda3/miniconda.sh
				bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
				rm ~/miniconda3/miniconda.sh
				echo After installing Miniconda3, close and reopen your terminal application.
			fi
		else
			echo "miniconda3 installation is canceled."
		fi
		;;
	setup-conda)
		read -p "Do you want to install miniconda3? (y/N): " confirm
		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			# Check if conda is available
			if command -v conda &>/dev/null; then
				source "$(conda info --base)/etc/profile.d/conda.sh"
				conda activate base
				conda config --add channels bioconda
				conda config --add channels conda-forge
				conda config --set channel_priority strict
			else
				echo "Error: conda command not found."
				exit 1
			fi
		else
			echo "miniconda3 bioconda setup is canceled."
		fi
		;;
	uninstall)
		echo "Removing the following conda environments: polap, polap-fmlrc, getorganelle"
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to uninstall all conda environments for polap? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			source "$(conda info --base)/etc/profile.d/conda.sh"
			if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
				echo "You're not in the base environment. Chaniging 'base'..."
				conda activate base
			fi

			if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
				if conda info --envs | awk '{print $1}' | grep -Fxq polap; then
					conda env remove -n polap
					# conda remove -y -n polap --all
				fi
				if conda info --envs | awk '{print $1}' | grep -Fxq polap-fmlrc; then
					conda env remove -n polap-fmlrc
				fi
				if conda info --envs | awk '{print $1}' | grep -Fxq getorganelle; then
					conda env remove -n getorganelle
				fi
				if conda info --envs | awk '{print $1}' | grep -Fxq pmat; then
					conda env remove -n pmat
				fi

				echo "conda deactivate if necessary"
			fi
		else
			echo "Uninstallation is canceled."
			echo "${help_message_uninstall}"
		fi
		;;
	install-polap)
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to install polap? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			# Check current conda environment
			# Initialize Conda for non-interactive shells
			source "$(conda info --base)/etc/profile.d/conda.sh"
			if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
				echo "You're not in the base environment. Chaniging 'base'..."
				conda activate base
			fi

			if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
				echo "You're in the base environment. Creating 'polap'..."
				if conda env list | awk '{print $1}' | grep -qx "polap"; then
					echo "ERROR: Conda environment 'polap' already exists."
				else
					conda create -y --name polap bioconda::polap
				fi
			else
				echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
				exit 1
			fi
		else
			echo "polap installation is canceled."
			echo "${help_message_install_polap}"
		fi
		;;
	download-polap-github)
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to download polap from github? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			if [[ -d "polap-${_polap_version}" ]]; then
				echo "ERROR: you already have polap-${_polap_version}"
				echo "  delete it if you want to redownload it."
				exit 0
			fi
			wget -q https://github.com/goshng/polap/archive/refs/tags/${_polap_version}.zip
			unzip -o -q ${_polap_version}.zip
			echo "polap github source is created at polap-${_polap_version}"
		else
			echo "polap download from github is canceled."
		fi
		;;
	delete-polap-github)
		echo "Deleting the following:"
		echo "  ${_polap_version}.zip and duplicates"
		echo "  polap-${_polap_version}"
		echo "  polap"
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to delete all the downloaded polap source? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			rm -f ${_polap_version}.zip*
			rm -rf polap-${_polap_version}
			rm -rf polap
			echo "polap download polap-${_polap_version} and its zipped file have been deleted."
		else
			echo "Deleting the polap download from github is canceled."
		fi
		;;
	install-fmlrc)
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to install polap-fmlrc? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			# Check current conda environment
			# Initialize Conda for non-interactive shells
			source "$(conda info --base)/etc/profile.d/conda.sh"
			if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
				echo "You're not in the base environment. Chaniging 'base'..."
				conda activate base
			fi

			if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
				echo "You're in the base environment. Creating 'polap-fmlrc'..."
				if conda env list | awk '{print $1}' | grep -qx "polap-fmlrc"; then
					echo "ERROR: Conda environment 'polap-fmlrc' already exists."
				else
					wget -q https://github.com/goshng/polap/archive/refs/tags/${_polap_version}.zip
					if [[ -s "${_polap_version}.zip" ]]; then
						unzip -o -q ${_polap_version}.zip
						cd polap-${_polap_version}
						conda env create -f src/polaplib/polap-conda-environment-fmlrc.yaml
					else
						echo "Error: no such file: ${_polap_version}.zip - no such polap version"
						echo "Suggestion: _polap_version=0.4.3.7.4 $0 $subcmd1"
					fi
				fi
			else
				echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
				exit 1
			fi
		else
			echo "polap installation is canceled."
			echo "${help_message_install_fmlrc}"
		fi
		;;
	patch-polap)
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to replace conda env polap with the version ${_polap_version}? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			if conda env list | awk '{print $1}' | grep -qx "polap"; then
				echo "Updating the Polap Conda environment to version ${_polap_version}."
				wget -q https://github.com/goshng/polap/archive/refs/tags/${_polap_version}.zip
				if [[ -s "${_polap_version}.zip" ]]; then
					unzip -o -q ${_polap_version}.zip
					cd polap-${_polap_version}/src
					bash polaplib/polap-build.sh >../build.sh
					cd ..
					PREFIX="$(conda info --base)/envs/polap" bash build.sh
				else
					echo "Error: no such file: ${_polap_version}.zip - no such polap version"
					echo "Suggestion: _polap_version=0.4.3.7.4 $0 $subcmd1"
				fi
			else
				echo "Error: You do not have polap environment. Please activate polap before running this script."
			fi
		else
			echo "polap patch is canceled."
			echo "${help_message_patch_polap}"
		fi
		;;
	bleeding-edge-polap)
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to update conda env polap with the latest polap github? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			if conda env list | awk '{print $1}' | grep -qx "polap"; then
				echo "Updating the Polap Conda environment to its most current state on GitHub ensures access to the latest features and updates."
				if [[ -d "polap" ]]; then
					rm -rf polap
				fi
				git clone --quiet https://github.com/goshng/polap.git
				cd polap/src
				bash polaplib/polap-build.sh >../build.sh
				cd ..
				PREFIX="$(conda info --base)/envs/polap" bash build.sh
			else
				echo "Error: You do not have polap environment. Please activate polap before running this script."
			fi
		else
			echo "polap bleeding-edge version update is canceled."
			echo "${help_message_bleeding_edge_polap}"
		fi
		;;
	local-edge-polap)
		if [[ "${_arg2}" == arg2 ]]; then
			echo "Help: ${subcmd1} <polap github base path>"
			echo "  $0 ${subcmd1} polap/github"
			exit 0
		fi
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to update conda env polap with the local polap source? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			if conda env list | awk '{print $1}' | grep -qx "polap"; then
				echo "Updating the Polap Conda environment to the local source at ${_arg2}."
				if [[ -d "${_arg2}/src" ]]; then
					cd "${_arg2}/src"
					bash polaplib/polap-build.sh >../build.sh
					cd ..
					PREFIX="$(conda info --base)/envs/polap" bash build.sh
				else
					echo "Error: The polap base path does not have src directory: ${_arg2}"
				fi
			else
				echo "Error: You do not have polap environment. Please activate polap before running this script."
			fi
		else
			echo "polap local version update is canceled."
		fi
		;;
	test-polap)
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to test polap? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			if [[ ! -d "polap-${_polap_version}" ]]; then
				echo "ERROR: no such folder: polap-${_polap_version}"
				echo "  fetch the polap github source"
				echo "$0 download-polap-github"
				exit 0
			fi
			cd polap-${_polap_version}/test
			source "$(conda info --base)/etc/profile.d/conda.sh"
			if [[ "$CONDA_DEFAULT_ENV" != "polap" ]]; then
				echo "You're not in the polap environment. Chaniging 'polap'..."
				conda activate polap
			fi

			if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then
				polap assemble --test
			else
				echo "ERROR: no such conda environment: polap"
				echo "$0 install-polap"
			fi
		else
			echo "polap test is canceled."
		fi
		;;
	install-bandage)
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to install bandage? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			wget -q https://github.com/rrwick/Bandage/releases/download/v0.8.1/Bandage_Ubuntu_static_v0_8_1.zip
			unzip Bandage_Ubuntu_static_v0_8_1.zip
			mkdir -p bin
			cp Bandage "bin"
			source <(echo 'export PATH="$PWD/bin:$PATH"')
		else
			echo "bandage installation is canceled."
			echo "${help_message_install_bandage}"
		fi
		;;
	install-getorganelle)
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to install getorganelle? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			# Check current conda environment
			# Initialize Conda for non-interactive shells
			source "$(conda info --base)/etc/profile.d/conda.sh"
			if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
				echo "You're not in the base environment. Chaniging 'base'..."
				conda activate base
			fi

			if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
				echo "You're in the base environment. Creating 'getorganelle'..."
				if conda env list | awk '{print $1}' | grep -qx "getorganelle"; then
					echo "ERROR: Conda environment 'getorganelle' already exists."
				else
					conda create -y --name getorganelle bioconda::getorganelle
					conda activate getorganelle
					get_organelle_config.py --add embplant_pt,embplant_mt
				fi
			else
				echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
				exit 1
			fi
		else
			echo "getorganelle installation is canceled."
			echo "${help_message_install_getorganelle}"
		fi
		;;
	install-cflye)
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to install cflye? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			# Check current conda environment
			# Initialize Conda for non-interactive shells
			source "$(conda info --base)/etc/profile.d/conda.sh"
			if [[ "$CONDA_DEFAULT_ENV" != "polap" ]]; then
				echo "You're not in the polap environment. Chaniging 'polap'..."
				conda activate polap
			fi

			if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then
				echo "You're in the polap environment. Installing 'cflye'..."
				conda install -y goshng::cflye
			else
				echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate polap before running this script."
				exit 1
			fi
		else
			echo "cflye or read-coverage filtering version Flye installation is canceled."
		fi
		;;
	install-dflye)
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to install dflye? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			# Check current conda environment
			# Initialize Conda for non-interactive shells
			source "$(conda info --base)/etc/profile.d/conda.sh"
			if [[ "$CONDA_DEFAULT_ENV" != "polap" ]]; then
				echo "You're not in the polap environment. Chaniging 'polap'..."
				conda activate polap
			fi

			if [[ "$CONDA_DEFAULT_ENV" == "polap" ]]; then
				echo "You're in the polap environment. Installing 'dflye'..."
				conda install -y goshng::dflye
			else
				echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate polap before running this script."
				exit 1
			fi
		else
			echo "dflye or read-coverage filtering version Flye installation is canceled."
		fi
		;;
	install-tippo)
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to install tippo? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			# Check current conda environment
			# Initialize Conda for non-interactive shells
			source "$(conda info --base)/etc/profile.d/conda.sh"
			if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
				echo "You're not in the base environment. Chaniging 'base'..."
				conda activate base
			fi

			if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
				echo "You're in the base environment. Creating 'tippo'..."
				if conda env list | awk '{print $1}' | grep -qx "tippo"; then
					echo "ERROR: Conda environment 'tippo' already exists."
				else
					conda create -y --name tippo bioconda::tipp
					conda activate tippo
				fi
			else
				echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
				exit 1
			fi
		else
			echo "tippo installation is canceled."
			echo "${help_message_install_tippo}"
		fi
		;;
	download-pmat)
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to download pmat? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			wget https://github.com/bichangwei/PMAT/archive/refs/tags/v1.5.3.tar.gz
			tar -zxvf v1.5.3.tar.gz
			source <(echo 'export PATH="$PWD/PMAT-1.5.3/bin:$PATH"')
			cd PMAT-1.5.3/bin
			chmod +x PMAT
		else
			echo "pmat download is canceled."
		fi
		;;
	install-pmat)
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to install pmat? (y/N): " confirm
		else
			confirm="yes"
		fi

		if [[ "${confirm,,}" == "yes" || "${confirm,,}" == "y" ]]; then
			# Check current conda environment
			# Initialize Conda for non-interactive shells
			source "$(conda info --base)/etc/profile.d/conda.sh"
			if [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
				echo "You're not in the base environment. Chaniging 'base'..."
				conda activate base
			fi

			if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
				echo "You're in the base environment. Creating 'pmat'..."
				if conda env list | awk '{print $1}' | grep -qx "pmat"; then
					echo "ERROR: Conda environment 'pmat' already exists."
				else
					conda create -y --name pmat apptainer nextdenovo canu blast
					# conda activate pmat
					# conda install -y conda-forge::apptainer
					# conda install -y bioconda::canu
					# conda install -y bioconda::nextdenovo
					# conda install -y bioconda::blast
					# conda install -y bioconda::aspera-cli
				fi
			else
				echo "Error: You're in the '$CONDA_DEFAULT_ENV' environment. Please activate base before running this script."
				exit 1
			fi
		else
			echo "pmat installation is canceled."
		fi
		;;
	mkdir-all)
		handled=1
		if [[ "${opt_y_flag}" == false ]]; then
			read -p "Do you want to create all species folders? (y/N): " confirm
		else
			confirm="yes"
		fi

		case "$confirm" in
		[yY] | [yY][eE][sS])
			echo "creating all species folder $i ..."
			for i in "${Sall[@]}"; do
				mkdir ${i}
			done
			;;
		*)
			echo "Creating all species folders is canceled."
			;;
		esac
		;;
	rm-empty)
		handled=1
		find . -type d -empty -delete
		echo "Deleting empty folders ... done."
		;;
	rm)
		handled=1
		read -p "Do you really want to delete all? (YES/NO): " confirm
		case "$confirm" in
		YES)
			read -p "Really, type in YES? (YES): " confirm
			if [[ "$confirm" == "YES" ]]; then
				for i in "${Sall[@]}"; do
					echo "deleting folder $i ..."
					rm -rf ${i}
				done
			fi
			;;
		*)
			echo "Deleting all is canceled."
			;;
		esac
		;;
	esac

	[[ $handled -eq 1 ]] && return 0 || return 1
}
