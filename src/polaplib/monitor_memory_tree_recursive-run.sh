_script=$1
# Step 1: Start your command and get its PID
( bash ${_script} & echo $! > pid.txt )

# Step 2: Start the monitor
bash monitor_memory_tree_recursive.sh $(cat pid.txt) memlog.csv

