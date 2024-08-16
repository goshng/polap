a2ps run-polap2.sh --landscape --columns=1 --chars-per-line=132 --line-numbers=1 -o 1.ps
ps2pdf 1.ps
mv 1.pdf ~/samba/carex-working/jelee/mtdna/data/4.pdf

a2ps --landscape --columns=1 --chars-per-line=132 -o 1.ps run-polap-pairs.R
ps2pdf 1.ps
mv 1.pdf ~/samba/carex-working/jelee/mtdna/data/5.pdf
