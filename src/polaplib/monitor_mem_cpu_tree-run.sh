( bash execute_nextdenovo.sh & echo $! > pid.txt )
bash monitor_mem_cpu_tree.sh $(cat pid.txt) memcpu_log.csv

