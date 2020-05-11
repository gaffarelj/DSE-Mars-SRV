import multiprocessing as mp
import numpy as np
def test(n):
    return np.array([n,n*2])

if __name__ == '__main__':
    with mp.Pool(5) as p:
        print(np.sum(p.map(test, [1, 2, 3])))
