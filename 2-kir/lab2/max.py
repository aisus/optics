import numpy
import math
import cmath
import matplotlib.pyplot as pl


def f(x):
    return 1 / (1 + x**2)


def gaussian_beam(x):
    return numpy.exp(-(x * x))


def amplitude(func):
    a = []
    for i in func:
        a.append(abs(i))
    return a


def phase(func):
    a = []
    for i in func:
        a.append(cmath.phase(i))
    return a


def append_zeros(func,
                 number_of_zeros):  # заполняли вектор необходимым числом нулей
    f_with_zeros = numpy.zeros(2 * number_of_zeros + len(func))
    for i in range(len(func)):
        f_with_zeros[number_of_zeros + i] = func[i]
    return f_with_zeros


def extract_function(func, N, number_of_zeros):
    f_without_zeros = []
    for i in range(number_of_zeros, number_of_zeros + N):
        f_without_zeros.append(func[i])
    return f_without_zeros


# ЗАДАНИЕ №1
def finite_fourier_transform(func, N, M, hx, number_of_zeros):
    f1 = append_zeros(func, int(number_of_zeros))
    tmp1 = f1[0:1023]
    tmp2 = f1[1024:2048]
    f1 = [*tmp2, *tmp1]

    # f1 = numpy.fft.fftshift(f1)  #поменяли местаи ([ 0.,  1.,  2.,  3.,  4., -5., -4., -3., -2., -1.]),
    # стало array([-5., -4., -3., -2., -1.,  0.,  1.,  2.,  3.,  4.])

    F = hx * numpy.fft.fft(
        f1
    )  # Выполнить БПФ от 𝐟 и умножить результат на шаг ℎ𝑥, получив вектор F
    tmp1 = F[0:1023]
    tmp2 = F[1024:2048]
    F = [*tmp2, *tmp1]

    # F = numpy.fft.ifftshift(F) # Разбить вектор 𝐅 на две половины и поменять их местами.
    F = extract_function(
        F, N, int(number_of_zeros)
    )  # «Вырезать» центральную часть вектора 𝐅, оставив 𝑁 элементов.

    return F


def my_ftt(func, a, b, N, hx, hu):
    sumf = [0j] * M
    for i in range(M):
        u = -b + i * hu
        for k in range(N):
            r = cmath.exp(-2 * math.pi * u * (-a + k * hx) * 1j)
            f_mas = func(-a + k * hx)
            sumf[i] += f_mas * r
        sumf[i] *= hx
    return sumf


def jordans(M, hu):
    sumf = [0j] * M
    for i in range(M):
        u = -b + i * hu
        sumf[i] = math.pi * numpy.exp(2 * math.pi * (-1) * abs(u))
    return sumf


if __name__ == '__main__':
    a = 5
    N = 256
    M = 2048
    b = N**2 / (4 * a * M)
    number_of_zeros = (M - N) / 2
    hx = 2 * a / N  # шаг по x
    hu = 2 * b / M  # шаг по u

    x = numpy.linspace(-a, a, N)
    x1 = numpy.linspace(-b, b, N)
    x2 = numpy.linspace(-b, b, M)
    x3 = numpy.linspace(-a, a, M)

    gaus_beam = gaussian_beam(x)
    gaus_beam_f = finite_fourier_transform(gaus_beam, N, M, hx,
                                           number_of_zeros)
    gaus_beam_f2 = my_ftt(gaussian_beam, a, b, N, hx, hu)

    # my_beam = f(x)
    # my_beam_f = finite_fourier_transform(my_beam, N, M, hx, number_of_zeros)
    # my_beam_f2 = my_ftt(f, a, b, N, hx, hu)

    # ЗАДАНИЕ №2

    pl.subplot(131)
    pl.plot(x, gaussian_beam(x))
    pl.title("Пучок Гаусса")
    pl.grid()

    # ЗАДАНИЕ №3 , 4

    pl.subplot(132)
    pl.title('Амплитуда')
    pl.plot(x, amplitude(gaus_beam), label='Исходный')
    pl.plot(x1, amplitude(gaus_beam_f), label='Фурье встр.')
    pl.plot(x2, amplitude(gaus_beam_f2), label='Фурье числ.')
    pl.grid()
    pl.legend()

    pl.subplot(133)
    pl.title('Фаза')
    pl.plot(x, phase(gaus_beam), label='Исходный')
    pl.plot(x1, phase(gaus_beam_f), label='Фурье встр.')
    pl.plot(x2, phase(gaus_beam_f2), label='Фурье числ.')
    pl.gca().set_ylim(-math.pi, math.pi)
    pl.grid()
    pl.legend()
    pl.show()

    # ЗАДАНИЕ № 6
    my_beam_f = finite_fourier_transform(f(x), N, M, hx, number_of_zeros)
    pl.subplot(131)
    pl.plot(x, f(x))
    pl.title("Исходный мой пучок")
    pl.grid()
    pl.subplot(132)
    pl.title('Амплитуда')
    pl.plot(x, amplitude(f(x)), label='Исходный')
    pl.plot(x1, amplitude(my_beam_f), label='Фурье встр.')
    # pl.plot(x2, amplitude(f(x)), label='Моё')

    pl.plot(x2, amplitude(jordans(M, hu)), label='Аналитика', color='c')
    pl.grid()
    pl.legend()
    pl.subplot(133)
    pl.title('Фаза')
    pl.plot(x, phase(f(x)), label='Исходный')
    pl.plot(x1, phase(my_beam_f), label='Фурье встр.')
    # pl.plot(x2, phase(f(x)), label='Моё')
    pl.plot(x2, phase(jordans(M, hu)), label='Аналитика', color='c')
    # pl.gca().set_ylim(-math.pi, math.pi)
    pl.grid()
    pl.legend()

    pl.show()
