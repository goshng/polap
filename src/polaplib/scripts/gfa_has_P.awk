#!/usr/bin/awk -f
BEGIN{FS="\t"; found=0}
$1=="P"{found=1; exit}
END{ exit found?0:1 }
