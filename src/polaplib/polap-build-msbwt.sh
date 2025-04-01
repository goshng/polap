_short_read1="${1}"
_msbwt_dir="${2}"

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate polap-fmlrc

rm -rf "${_msbwt_dir}"
cat "${_short_read1}" |
	awk 'NR % 4 == 2' | sort | tr NT TN |
	ropebwt2 -LR |
	tr NT TN |
	msbwt convert "${_msbwt_dir}" \
		>/dev/null

conda deactivate
