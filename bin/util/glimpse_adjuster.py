#!/usr/bin/env python3

import random
import os
import sys
import random

backupt = 41
if random.randint(0, 1) in [1]:
    print("success!", file=sys.stderr)
    pass
else:
    print("next time", file=sys.stderr)
    os._exit(0)

try:
    nthreads = int(sys.argv[1])
    ctrl_file = sys.argv[2]
    glances_file = f"{ctrl_file}.glances"
    new_cores = nthreads
    print("a")
    os.system(
        f"timeout 10 glances --stdout-csv cpu.cpucore,cpu.total >> {glances_file} || sleep .1  "
    )
    cpu_info = os.popen(f"tail -n 1 {ctrl_file}.glances").readline().rstrip().split(",")
    curr_use = float(str(cpu_info[0]).replace("'", ""))
    max_c = float(str(cpu_info[1]).replace("'", ""))
    print("tt")
    delta_used = max_c - curr_use
    pct_remaining = float(curr_use) / float(max_c)
    if pct_remaining < 0.25:
        new_cores = nthreads
    else:
        new_cores = nthreads + int(delta_used * 0.4)
    print("ff")
    os.system(f"echo {new_cores} > {ctrl_file}")
    os._exit(0)

except Exception as e:
    print(e, file=sys.stderr)
    os.system("pkill -n glances || sleep .1")
    try:
        os.system(f"echo {backupt} > {ctrl_file}")
        print("Q")
        os._exit(0)
    except Exception as e:
        print(e, file=sys.stderr)
        print("t")
        os.system(f"echo {backupt}  > {ctrl_file}")
        os._exit(0)

os.system(f"echo {backupt} > {ctrl_file}")

os._exit(0)
