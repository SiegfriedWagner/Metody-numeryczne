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


# @jit(parallel=True)
def haar2d(original_signal: np.ndarray, n: int):
    original_signal = original_signal.copy()
    rows = original_signal.shape[0]
    cols = original_signal.shape[1]
    for _ in range(n):
        new_col = cols // 2
        for row in range(rows):
            a = calc_next_am(original_signal[row, :cols])
            d = calc_next_dm(original_signal[row, :cols])
            original_signal[row, :new_col] = a
            original_signal[row, new_col:cols] = d
        cols = new_col
        new_row = rows // 2
        for col in range(cols):
            a = calc_next_am(original_signal[:rows, col])
            d = calc_next_dm(original_signal[:rows, col])
            original_signal[:new_row, col] = a
            original_signal[new_row:rows, col] = d
        rows = new_row
    return original_signal


def recreate_signal2d(signal: np.ndarray, n: int):
    k = pow(2, n)
    rows = signal.shape[0] // k
    cols = signal.shape[1] // k
    for _ in range(n):
        new_rows = rows * 2
        col_array = np.zeros(new_rows)
        for col in range(cols):
            col_array[::2] = signal[:rows, col]
            col_array[1::2] = signal[:rows, col]
            col_array[::2] += signal[rows:new_rows, col]
            col_array[1::2] -= signal[rows:new_rows, col]
            signal[:new_rows, col] = col_array
        rows = new_rows
        new_cols = cols * 2
        row_array = np.zeros(new_cols)
        for row in range(rows):
            row_array[::2] = signal[row, :cols]
            row_array[1::2] = signal[row, :cols]
            row_array[::2] += signal[row, cols:new_cols]
            row_array[1::2] -= signal[row, cols:new_cols]
            signal[row, :new_cols] = row_array
        cols = new_cols
    return signal


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
