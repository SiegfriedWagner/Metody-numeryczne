import numpy as np
from numba import jit


def calc_next_am(signal: np.ndarray):
    return (signal[::2] + signal[1::2]) / np.sqrt(2)


def calc_next_dm(signal: np.ndarray):
    return (signal[::2] - signal[1::2]) / np.sqrt(2)


def haar(original_signal: np.ndarray, n: int):
    result = []
    a_n = original_signal
    while n > 0:
        n -= 1
        result.insert(0, calc_next_dm(a_n))
        a_n = calc_next_am(a_n)
    result.insert(0, a_n)
    return result


def gen_v(k: int, n: int) -> np.array:
    k2 = pow(2, k)
    arr = np.zeros((n // k2, n))
    for n in range(n // k2):
        arr[n, k2*n:k2*(n+1)] = 1
    return arr / pow(2, k/2)


def gen_w(k: int, n: int) -> np.array:
    k2 = pow(2, k)
    arr = np.zeros((n // k2, n))
    for n in range(n // k2):
        arr[n, k2 * n : k2 * (n+1) - k] = 1
        arr[n, k2 * (n+1) - k : k2 * (n+1)] = -1
    return arr / pow(2, k/2)


@jit(parallel=True)
def recreate_signal(a, details, orginal_signal_size):
    k = int(np.log2(orginal_signal_size // a.size))
    k2 = int(pow(2, k))
    recreated = np.zeros(orginal_signal_size)
    print("Recreating: A")
    helping = np.ones(k2) / pow(2, k / 2)
    for i in range(a.size):
        recreated[i * k2: (i + 1) * k2] += helping * a[i]
    for detail in details:
        print(f"Recreating: Detail")
        k = int(np.log2(orginal_signal_size // detail.size))
        k2 = int(pow(2, k))
        helping = np.ones(int(k2))
        helping[k2 // 2:] = -1
        helping /= pow(2, k / 2)
        for i in range(detail.size):
            recreated[i * k2: (i + 1) * k2] += helping * detail[i]
    return recreated
