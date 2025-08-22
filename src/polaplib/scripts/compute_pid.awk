# compute_pid.awk
# Input: two lines from `seqkit fx2tab` => "<name>\t<aligned_sequence>"
# Vars: RAW=1 -> print numeric only (e.g., 99.932); RAW=0 -> "99.932%"
# scripts/compute_pid.awk
# Compute global percent identity from 2 aligned sequences (seqkit fx2tab output)
#
# scripts/compute_pid.awk
# Compute global percent identity from 2 aligned sequences (seqkit fx2tab output)

NR==1 {
    s1 = $2
    next
}
NR==2 {
    s2 = $2
}
END {
    ungap = 0
    matches = 0
    len1 = length(s1)
    len2 = length(s2)
    if (len1 != len2) {
        printf("NA\n")
        exit 1
    }
    for (i = 1; i <= len1; i++) {
        a = substr(s1, i, 1)
        b = substr(s2, i, 1)
        A = toupper(a)
        B = toupper(b)
        if (a != "-" && b != "-") {
            ungap++
            if (A == B) {
                matches++
            }
        }
    }
    if (ungap > 0) {
        pid = 100.0 * matches / ungap
    } else {
        pid = 0.0
    }
    if (RAW) {
        printf("%.3f\n", pid)
    } else {
        printf("%.3f%%\n", pid)
    }
}
