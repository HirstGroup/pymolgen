import sys,os

import matplotlib.pyplot as plt
import numpy as np


input = sys.argv[1]
procs = int(sys.argv[2])

weights = []

with open(input) as infile:
    for line in infile:
        a = line.strip('\n').split()
        a = [int(i) for i in a]
        weights.extend(a)

weights = sorted(weights)
weights.reverse()


def assign(weights, procs):

    ass = [0] * len(weights)

    totals = {i: 0.0 for i in range(procs)}
    for i, w in enumerate(weights):

        min_proc = min(totals, key=totals.get)
        ass[i] = min_proc
        totals[min_proc] += w

    return ass

ass = assign(weights, procs)

plt.subplot(221)
plt.plot(weights)

weight_per_proc = []
for i in range(procs):

    total = 0.0
    for j, w in enumerate(weights):
        if ass[j] == i:
            total += w
    weight_per_proc.append(total)
   
plt.subplot(222)
plt.bar(range(procs), weight_per_proc)


plt.show()