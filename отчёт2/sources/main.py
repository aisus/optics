from typing import List
import numpy as np
import matplotlib.pyplot as plt


def gauss_beam(x):
    return np.exp(-np.pi * (x**2))


def custom_func(x: float):
    # if abs(x) > 0.5:
    #     return 0
    # return (x - 1) / 4
    if x > 3 or x < -1:
        return 0
    return 1


def custom_func_fourier(x: float, xi: float):
    # return (2 * np.sqrt(2 / np.pi) * np.sin(xi) * np.cos(xi) * (np.cos(xi) + 1j * np.sin(xi))) / xi
    # return np.exp(-2*np.pi*1j*xi) * np.sin(4 * np.pi * xi) / np.pi * xi
    
    flags = np.abs(xi) < 1e-6
    around_zero = 4 * 1j 
    #around_zero = 1j 
    others = -1/(2 * np.pi * xi) * (np.exp(-6 * np.pi * x * 1j * xi) - (np.exp(2 * np.pi * x * 1j * xi)))

    return others * (1-flags) + around_zero * flags


def fourier_builitn(f, step: float, xs: List[float]) -> List[float]:
    fs = list(map(lambda x: f(x), xs))
    fs = np.fft.fftshift(fs)
    fourier = list(map(lambda x: x * step, np.fft.fft(fs)))
    return np.fft.fftshift(fourier)


def fourier_exp(u: float, x: float) -> float:
    return np.exp(-2 * np.pi * 1j * u * x)


def fourier_calc(f, step: float, xs: [float],
                 us: List[float]) -> List[float]:
    fourier = []
    for u in us:
        newx = 0
        for x in xs:
            newx += f(x) * fourier_exp(u, x)
        newx *= step
        fourier.append(newx)
    return fourier


def plot_calculus(f, a: float, b: float, step_count: int, fig_name):
    step_src = abs(b - a) / step_count
    c = step_count / (2 * abs(b - a))
    step = (2 * abs(c)) / step_count
    us = np.arange(-c, c, step)
    xs = np.arange(a, b, step_src)
    fourier = fourier_calc(f, step_src, xs, us)
    phase = np.angle(fourier)
    amplitude = np.abs(fourier)
    #print(f"Calculus => {fourier[:4]}")
    plt.figure(f"Phase {fig_name}")
    plt.plot(us, phase)
    plt.grid()
    plt.figure(f"Amplitude {fig_name}")
    plt.plot(us, amplitude)
    plt.grid()


def plot_fft(f, a: float, b: float, step_count: int, fig_name):
    step_src = (b - a) / step_count
    c = step_count / (2 * (abs(a) + abs(b)))
    step = (2 * abs(c)) / step_count
    #  step = (abs(a) + abs(b)) / step_count
    xs = np.arange(a, b, step_src)
    nxs = np.arange(-c, c, step)
    fourier = fourier_builitn(f, step_src, xs)
    phase = np.angle(fourier)
    amplitude = np.abs(fourier)
    plt.figure(f"Phase {fig_name}")
    plt.plot(nxs, phase)
    plt.grid()
    plt.figure(f"Amplitude {fig_name}")
    plt.plot(nxs, amplitude)
    plt.grid()


def plot_complex_function(f, a, b, step_count, fig_name):
    c = step_count / (2 * (abs(a) + abs(b)))
    step = (2 * abs(c)) / step_count
    nxs = np.arange(-c, c, step)
    us = np.arange(-c, c, step)
    ys = f(nxs, us)
    #ys = f(nxs)
    amplitude = np.abs(ys)
    phase = np.angle(ys)
    plt.figure(f"Phase {fig_name}")
    plt.plot(nxs, phase)
    plt.grid()
    plt.figure(f"Amplitude {fig_name}")
    plt.plot(nxs, amplitude)
    plt.grid()


if __name__ == "__main__":
    #plot_complex_function(gauss_beam, -5, 5, 64, "gauss")
    #plot_fft(gauss_beam, -5, 5, 256, "gauss")
    #plot_calculus(gauss_beam, -5, 5, 256, "gauss")
    plot_fft(custom_func, -2, 2, 512, "custom")
    #plot_calculus(custom_func, -5, 5, 512, "custom")
    #plot_complex_function(custom_func_fourier, -5, 5, 512, "custom")
   
    plt.show()
