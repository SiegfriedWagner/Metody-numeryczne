import numpy as np
import matplotlib.pyplot as plt

freq = 4
def fun(x):
    return np.abs(np.mod(x, freq) - 2) - 1

if __name__ == '__main__':
    sampling_frequency = 100
    range_end = 8
    print(f"Generowanie sygnału trójkątnego o okresie {freq} zakresie od 0 do {range_end} i częstotliwości {1 / freq}[Hz]")
    samples = range_end * sampling_frequency
    x = np.linspace(0, range_end, samples)
    y = fun(x)

    ft = np.fft.fft(y)
    freqfft_all = np.fft.fftfreq(ft.size, d=1 / sampling_frequency)
    plt.plot(freqfft_all, ft.real)
    power, freqfft = max((element for element in zip(ft, freqfft_all.real) if element[1] > 0), key=lambda pair: pair[0])
    freqfft = abs(freqfft)
    print(f"Most freq [Hz]: {freqfft}")
    periodfft = 1 / (freqfft)
    print(f"Period[s]: {periodfft}")

