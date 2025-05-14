echo "Peak memory usage: $(awk -F',' 'NR>1{if($2>max) max=$2} END{print max}' memlog.csv) KB"

awk -F',' 'NR>1{if($2>max) max=$2} END{
    mb = max / 1024
    gb = max / 1048576
    printf "Peak memory usage: %.2f MB (%.2f GB)\n", mb, gb
}' memlog.csv



