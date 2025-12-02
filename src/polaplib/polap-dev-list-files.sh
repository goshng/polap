# List all files under current tree, newest first (reverse sort)
find $1 -type f -printf '%T@ %s %p\n' |
	sort -n |
	awk '{ctime=strftime("%Y-%m-%d %H:%M:%S", $1); printf "%s\t%10d %s\n", ctime, $2, $3}'
