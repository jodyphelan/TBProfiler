import sys
import re

positions = [int(x) for x in sys.argv[1].split(",")]
maxp = max(positions)
minp = min(positions)

for l in sys.stdin:
    if l[0]=="@":
        sys.stdout.write(l)
        continue
    arr = l.rstrip().split()
    if "S" in arr[5]: continue
    start = int(arr[3])
    end = start+sum([int(x[:-1]) for x in re.findall("\d+M",arr[5])])
    if start<minp-10 and end>maxp+10:
        sys.stdout.write(l)
