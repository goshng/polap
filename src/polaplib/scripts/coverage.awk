# coverage.awk
$4=="plus"{cov[$3]+=$5} END{for(k in cov){printf "%.3f\n", cov[k]/k}}