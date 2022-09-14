from functools import partial
from multiprocessing import Pool

def f(a, b, c):
    print("{} {} {}".format(a, b, c))

def main():
    iterable = [1, 2, 3, 4, 5]
    pool = Pool()
    a = "hi"
    b = "there"
    func = partial(f, a, b)
    pool.map(func, iterable)
    pool.close()
    pool.join()

if __name__ == "__main__":
    main()