from multiprocessing import Pool
from itertools import product
from functools import partial
import time

def test(a, b, c=None, d=None):
    time.sleep(0.1)
    return a + d

a = 0
b=0
c=None
func = partial(test,a=a, b=b, c=c)

p = Pool(8)
outputs = p.map(func, range(10) )

print(outputs)



print('DONE')