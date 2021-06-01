import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from haar import *
import scipy as sc
from numba import jit
from function import fun, period
plt.rcParams['figure.figsize'] = [12, 8]
plt.rcParams['figure.dpi'] = 100 # 200 e.g. is really fine, but slower

# GLOBAL PARAMS
sampling_frequency = 100
range_begin = 0
range_end = 8
samples = (range_end - range_begin) * sampling_frequency

# BASIC SIGNAL
x = np.linspace(range_begin, range_end, samples)
y = fun(x)
plt.plot(x, y)
plt.ylabel("f(t)")
plt.xlabel("t[s]")
plt.title(f"Przebieg sygnału w funkcji czasu (częstotliwość próbkowania {sampling_frequency}Hz)")
plt.show(block=False)

# FOURIER TRANSFORM
ft = np.fft.fft(y)
freqfft_all = np.fft.fftfreq(ft.size, d=1/sampling_frequency)
plt.plot(freqfft_all, ft.real)
plt.ylabel("F[f(x)]")
plt.xlabel("Częstotliwość [Hz]")
plt.title("Transformata fouriera funkcji f(x)")
power, freqfft = max((element for element in zip(ft, freqfft_all.real) if element[1] > 0), key=lambda pair: pair[0])
tfreqfft = abs(freqfft)
print("Fourier transform")
print(f"Most freq [Hz]: {tfreqfft}")
periodfft = 1 / (tfreqfft)
print(f"Period[s]: {periodfft}")
plt.show(block=False)

# RECREATING FUNCTION FROM FOURIER TRANSFORM
samples = sampling_frequency*range_end
terms = 3
locx = np.linspace(range_begin, periodfft, samples, endpoint=False)
a0 = 2. / periodfft * simps(y, x)
an = lambda n: 2.0 / periodfft * simps(y * np.cos(2. * np.pi * n * x / periodfft), locx)
bn = lambda n: 2.0 / periodfft * simps(y * np.sin(2. * np.pi * n * x / periodfft), locx)
# sum of the series
s = a0 / 2. + sum([an(k) * np.cos(2. * np.pi * k * x / periodfft) + bn(k) * np.sin(2. * np.pi *
                                                                             k * x / periodfft) for k in range(1, terms + 1)])
plt.plot(x, y, label="Orginalny sygnał")
plt.plot(x, s, label=f"Syngał odtworzony z {terms} pierwszych wyrazów szeregu fouriera")
plt.xlabel("x")
plt.ylabel("y=f(x)")
plt.legend(loc='best',prop={'size':10})
plt.title("Wykres sygnału i sygnału odtworzonego z częściowego szerego Fouriera")
plt.show()

# SIMPLE HAAR EXAMPLE
print("\n\nSimple haar example:")
a = np.arange(0, 8, 1)
print(f"Signal: {a}")
v1 = gen_v(1, 8)
print(f"V1:")
print(v1)
w1 = gen_w(1, 8)
print(f"W1:")
print(w1)
h = haar(a, 1)
print("Haar result")
print(f"a1: {h[0]}")
print(f"d1: {h[1]}")

A = np.zeros(8)
D = np.zeros(8)
for n in range(4):
    A += v1[n] * h[0][n]
    D += w1[n] * h[1][n]
print(f"A: {A}")
print(f"d: {D}")
print(f"A + D: {A + D}")

# HAAR SIGNAL COMPRESSION
for n in range(1, 6, 2):
    h = haar(y, n)
    recreated = recreate_signal(h[0], h[1:], y.size)
    recreated = np.array(recreated).reshape(x.shape)
    plt.plot(x, recreated, label=f"haar({n})")
plt.legend()
plt.show()
