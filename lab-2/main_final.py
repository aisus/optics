from typing import List
import numpy as np
import matplotlib.pyplot as plt


def gauss_beam(x):
    return np.exp(-np.pi * (x**2))


def bokunofunction(x):
    #  x = x * 4
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


#  def bokunofunction(u: float, x: float):
#  """bokunofunction
#  Returns analytics function according to my optionem.
#  :param u:
#  :type u: float
#  :param x:
#  :type x: float
#  """
#  return np.exp(-2 * np.pi * 1j )
#  pass


def calculus_exp(u: float, x: float) -> float:
    """calculus_exp
    Returns e^{-2 \\pi i u x} from fourier transform.
    :param u: u value
    :type u: float
    :param x: x value
    :type x: float
    :rtype: float
    """
    return np.exp(-2 * np.pi * 1j * u * x)


def analit_func(a, b, u):
    return (np.exp(-2 * np.pi * 1j * b * (u - 1)) - np.exp(
        -2 * np.pi * 1j * a * (u - 1))) / (-2 * np.pi * 1j * (u - 1)) + (
            (np.exp(-2 * np.pi * 1j * b *
                    (u + 2.5)) - np.exp(-2 * np.pi * 1j * a * (u + 2.5))) /
            (-2 * np.pi * 1j * (u + 2.5)))


def calculus_fourier(f, step: float, xs: [float],
                     us: List[float]) -> List[float]:
    fourier = []
    for u in us:
        #  val = integrate(lambda x: (f(x)
        # * calculus_exp(u, x)), xs[0], xs[-1],
        #  1024, trapezium)
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
    plt.plot(us, phase)
    plt.grid()
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
    step = (b - a) / step_count
    xs = np.arange(a, b, step)
    ys = f(xs)
    plt.figure(f"Amplitude {fig}")
    plt.plot(xs, ys)


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
    plot_fft(bokunofunction, -5, 5, 256, "byaka")
    plot_calculus(bokunofunction, -5, 5, 256, "byaka")
    plot_complex(lambda x: analit_func(-5, 5, x), -5, 5, 256, "byaka")
    plt.show()
