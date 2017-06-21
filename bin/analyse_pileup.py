import sys
from collections import defaultdict

positions = [x for x in sys.argv[1].split(",")]

res = {}
temp_refs = {}
for l in sys.stdin:
    arr = l.rstrip().split()
    if arr[1] in positions:
        res[arr[1]] = arr[4].upper()
        temp_refs[arr[1]] = arr[2]
print temp_refs
print res

new_res = {}
for p in positions:
    new_res[p] = []
    for x in res[p]:
        if x=="." or x==",":
            new_res[p].append(temp_refs[p])
        else:
            new_res[p][-1] = new_res[p][-1]+x

print new_res
counts = defaultdict(int)
for i in range(len(res[positions[0]])):
    k = tuple([res[x][i] for x in positions])
    counts[k]+=1

print counts
