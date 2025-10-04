outdir=$1
rm -rf $outdir-nextdenovo
rm -rf $outdir-pmat
pkill -f execute_nextdenovo.sh 
pkill -f nextDenovo
pkill -f minimap2-nd
