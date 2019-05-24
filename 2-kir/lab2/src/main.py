from typing import List
import numpy as np
import matplotlib.pyplot as plt


def gauss_beam(x):
    return np.exp(-np.pi * (x**2))


def bokunofunction(x):
    """bokunofunction
    Function that calculates:
    $f(x) = e^{2 \pi i x} + e^{-5 \pi i x}$
    """
    return np.exp(2 * np.pi * 1j * x) + np.exp(-5 * np.pi * 1j * x)


def finit_fourier(f, step: float, xs: List[float]) -> List[float]:
    """finit_fourier
    Returns finite Fourier transform using fft.
    :param f: function to process.
    :param step: step.
    :type step: float
    :param xs: List of xs.
    :type xs: List[float]
    """
    fs = list(map(lambda x: f(x), xs))
    fs = np.fft.fftshift(fs)
    fourier = list(map(lambda x: x * step, np.fft.fft(fs)))
    return np.fft.fftshift(fourier)


def calculus_exp(u: float, x: float) -> float:
    """calculus_exp
    Returns $e^{-2 \pi i u x}$ from fourier transform.
    :param u: u value
    :type u: float
    :param x: x value
    :type x: float
    :rtype: float
    """
    return np.exp(-2 * np.pi * 1j * u * x)


def analit_func(a, b, u):
    """analit_func
        $\dfrac{e^{-2 \pi i b (u-1)} - e^{-2 \pi i a (u-1)}}{-2 \pi i (u-1)} + \dfrac{((e^{-2 \pi i b
                    (u + 2.5)}) - e^{(-2 \pi i a (u + 2.5))})}{(-2 \pi i (u + 2.5))}$
    """
    return (np.exp(-2 * np.pi * 1j * b * (u - 1)) - np.exp(
        -2 * np.pi * 1j * a * (u - 1))) / (-2 * np.pi * 1j * (u - 1)) + (
            (np.exp(-2 * np.pi * 1j * b *
                    (u + 2.5)) - np.exp(-2 * np.pi * 1j * a * (u + 2.5))) /
            (-2 * np.pi * 1j * (u + 2.5)))


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


def trapezium(f, x, h):
    return (f(x) + f(x + h)) / 2.0


def simpson(f, x, h):
    return (f(x) + 4 * f(x + h / 2) + f(x + h)) / 6.0


def integrate(f, a, b, steps, meth):
    h = (b - a) / steps
    ival = h * sum(meth(f, a + i * h, h) for i in range(steps))
    return ival


def plot_calculus(f,
                  a: float,
                  b: float,
                  step_count: int,
                  fig: str = "Calculus"):
    step_src = abs(b - a) / step_count
    c = step_count / (2 * abs(b - a))
    step = (2 * abs(c)) / step_count
    us = np.arange(-c, c, step)
    xs = np.arange(a, b, step_src)
    fourier = calculus_fourier(f, step_src, xs, us)
    phase = np.angle(fourier)
    amplitude = np.abs(fourier)
    plt.figure(f"Phase {fig}")
    plt.grid()
    plt.plot(us, phase)
    plt.figure(f"Amplitude {fig}")
    plt.plot(us, amplitude)
    plt.grid()


def plot_fft(f, a: float, b: float, step_count: int, fig: str = "fft"):
    step_src = (b - a) / step_count
    c = step_count / (2 * (abs(a) + abs(b)))
    step = (2 * abs(c)) / step_count
    #  step = (abs(a) + abs(b)) / step_count
    xs = np.arange(a, b, step_src)
    nxs = np.arange(-c, c, step)
    fourier = finit_fourier(f, step_src, xs)
    phase = np.angle(fourier)
    amplitude = np.abs(fourier)
    plt.figure(f"Phase {fig}")
    plt.plot(nxs, phase)
    plt.grid()
    plt.figure(f"Amplitude {fig}")
    plt.plot(nxs, amplitude)
    plt.grid()


def plot_function(f, a, b, step_count, fig):
    """plot_function
    Plot function on the interval from a to b
    with specified step count.
    :param f: $f(x)$
    :param a: start of the interval.
    :param b: interval end.
    :param step_count: how many points needs to be plotted.
    :param fig: plot name.
    """
    step = (b - a) / step_count
    xs = np.arange(a, b, step)
    ys = f(xs)
    plt.figure(f"Amplitude {fig}")
    plt.plot(xs, ys)
    plt.grid()


def plot_complex(f, a, b, step_count, fig="Source"):
    c = step_count / (2 * (abs(a) + abs(b)))
    step = (2 * abs(c)) / step_count
    nxs = np.arange(-c, c, step)
    ys = f(nxs)
    amplitude = np.abs(ys)
    phase = np.angle(ys)
    plt.figure(f"Phase {fig}")
    plt.plot(nxs, phase)
    plt.grid()
    plt.figure(f"Amplitude {fig}")
    plt.plot(nxs, amplitude)
    plt.grid()


if __name__ == "__main__":
    #  plot_function(gauss_beam, -5, 5, 512, "frf")
    #  plot_fft(gauss_beam, -5, 5, 512, "frf")
    #  plot_calculus(gauss_beam, -5, 5, 512, "frf")
    #  plot_complex(bokunofunction, -5, 5, 256, "bokunosource")
    plot_complex(lambda x: analit_func(-5, 5, x), -5, 5, 256, "test")
    plot_fft(bokunofunction, -5, 5, 256, "test")
    #  plot_calculus(bokunofunction, -5, 5, 256, "bokunocalculus")
    plt.show()
