import sys
import os


# Continuously check the output
dot = "."
while True:
    ret_val = os.popen("""pcluster list-clusters | grep clusterStatus | cut -d '"' -f 4""").readline().rstrip()
    # Read a line of output from the process
    
    if ret_val not in ['CREATE_IN_PROGRESS']:
        print('Cluster Creation: '+ret_val)
        break
    else:
        print(ret_val+": "+dot)
        dot = dot+"."
        os.system('sleep 5')
