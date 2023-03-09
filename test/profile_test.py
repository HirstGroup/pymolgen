import os
import sys

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("I require a single argument - the python script to profile!")
        quit()

    cmd = f"python -m cProfile -o profile.out {' '.join(sys.argv[1:])}"
    cmd = cmd.replace('(', '\\(').replace(')', '\\)')

    print("Running\n")
    print(cmd)
    print()
    os.system(cmd)

    pwd = os.path.dirname(os.path.realpath(__file__))

    os.system(f"python {pwd}/gprof2dot.py -f pstats profile.out | dot -Tpng -o profile.png")

"""
import os
import sys
import pstats
from pstats import SortKey

to_profile = sys.argv[1]
if len(sys.argv) > 2:
    percentage = int(sys.argv[2])
else:
    percentage = 5

print(f"Profiling {to_profile}")
os.system(f"python -m cProfile -o profile.out {to_profile}")
os.system(f"gprof2dot -f pstats -n {percentage} profile.out | dot -Tpng -o profile.png && eog profile.png")
"""