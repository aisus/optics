from typing import List
import numpy as np
import matplotlib.pyplot as plt


def gauss_beam(x: float):
    return np.exp(-np.pi * (x**2))


def custom_func(x: float):
    if abs(x) > 0.5:
        return 0
    return (x - 1) / 4


def anal_integral(x: float, ksi: float):
    return np.exp(-6 * np.pi * 1j * ksi) - np.exp(-2 * np.pi * 1j * ksi)


def finit_fourier(f, step: float, xs: List[float]) -> List[float]:
    fs = list(map(lambda x: f(x), xs))
    fs = np.fft.fftshift(fs)
    fourier = list(map(lambda x: x * step, np.fft.fft(fs)))
    return np.fft.fftshift(fourier)


def calculus_exp(u: float, x: float) -> float:
    return np.exp(-2 * np.pi * 1j * u * x)


def calculus_fourier(f, step: float, xs: [float],
                     us: List[float]) -> List[float]:
    fourier = []
    for u in us:
        newx = 0
        for x in xs:
            newx += f(x) * calculus_exp(u, x)
        newx *= step
        fourier.append(newx)
    return fourier


def plot_calculus(f, a: float, b: float, step_count: int):
    step_src = (abs(b) + abs(a)) / step_count
    c = step_count / (2 * (abs(a) + abs(b)))
    step = (2 * abs(c)) / step_count
    us = np.arange(-c, c, step)
    xs = np.arange(a, b, step_src)
    print(xs[0], xs[-1])
    print(us[0], us[-1])
    fourier = calculus_fourier(f, step_src, xs, us)
    print(f"Calculus => {fourier[:4]}")
    phase = list(map(np.angle, fourier))
    amplitude = list(map(abs, fourier))
    plt.figure("Phase fft")
    plt.plot(us, phase)
    plt.grid()
    plt.figure("Amplitude fft")
    plt.plot(us, amplitude)
    plt.grid()


def plot_fft(f, a: float, b: float, step_count: int):
    step = (abs(a) + abs(b)) / step_count
    xs = np.arange(a, b, step)
    fourier = finit_fourier(f, step, xs)
    print(f"FFT => {fourier[:4]}")
    phase = list(map(np.angle, fourier))
    amplitude = list(map(abs, fourier))
    plt.figure("Phase fft")
    plt.plot(xs, phase)
    plt.grid()
    plt.figure("Amplitude fft")
    plt.plot(xs, amplitude)
    plt.grid()


if __name__ == "__main__":
    #plot_fft(gauss_beam, -5, 5, 256)
    #plot_calculus(gauss_beam, -5, 5, 256)
    plot_fft(custom_func, -5, 5, 256)
    plot_calculus(custom_func, -5, 5, 256)
    plt.show()
