from contextlib import contextmanager
from time import perf_counter

import numpy as np
from memory_profiler import memory_usage
from scipy.sparse import eye
from scipy.sparse.linalg import LinearOperator, eigs


@contextmanager
def catchtime(title:str = "") -> float:
    start = perf_counter()
    yield lambda: perf_counter() - start
    print(f'Time {title}: {perf_counter() - start:.3f} seconds')

def get_eigenvalues(A:np.ndarray, n:int) -> tuple:
    """Use Krylov to compute n lowest eigenvalues and eigenvectors of AV=VE

    Args:
        A (np.ndarray): Matrix of compute eigenvalues of
        n (int): Number of eigenvalues to return
    """
    w,v = eigs(A, k=n, which="SR")
    return w,v 
    

def main():
    n = 1000
    up = eye(n, n, k=1)
    down = eye(n, n, k=-1)
    get_eigenvalues(0.3*up+0.7*down,20)


if __name__ == "__main__":
    with catchtime():
        mem_usage = memory_usage(main)
    print(f'Memory usage (in chunks of .1 seconds): {mem_usage}')
    print(f'Maximum memory usage: {max(mem_usage)-min(mem_usage)}')
