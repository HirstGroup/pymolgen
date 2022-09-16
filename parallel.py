from multiprocessing import Pool
from itertools import product
from functools import partial
import time

def test(a, b):
    time.sleep(1)
    return a + b

a = 0

func = partial(test,a)
"""
p = Pool(8)
outputs = p.map(func, range(100) )
"""

for i in range(100):
    test(a,i)

print('DONE')
"""
c = []

for i in range(100000000):
    c.append(a + i)

print('DONE')
"""