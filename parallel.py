from multiprocessing import Pool
from itertools import product
from functools import partial

def test(a, b, c, d=None):
    print(a,b,c,d)
    return a, b


a = 1
b = 2
c = 3

func = partial(test,a,b,c)

p = Pool(8)
outputs = p.map(func, range(5) )

print(outputs)
