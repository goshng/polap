SR1=4C-pse_1.fastq.gz
SR2=4C-pse_2.fastq.gz

if [[ $(file -b --mime-type "$SR1") == "application/gzip" ]]; then
	gunzip -c "$SR1" >"${SR1%.gz}"
	SR1="${SR1%.gz}"
fi

echo $SR1
echo $SR2
