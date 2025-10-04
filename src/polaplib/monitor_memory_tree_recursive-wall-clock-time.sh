awk -F',' 'NR==2 {start=$1} END {
    delta = $1 - start
    hrs = int(delta / 3600)
    min = int((delta % 3600) / 60)
    sec = delta % 60
    printf("Wall clock time: %02d:%02d:%02d\n", hrs, min, sec)
}' memlog.csv

