import numpy as np


def normalize(v):
    norm = np.linalg.norm(v, ord=1)
    if norm == 0:
        norm = np.finfo(v.dtype).eps
    return v / norm


def reflection(e, n):
    reflection_point = e - 2 * np.dot(e, n) * n
    print(reflection_point)
    return reflection_point


def refraction(e, n, n1, n2):
    under_root = 1 - ((n1 ** 2) / (n2 ** 2)) * (1 - np.dot(e, n) ** 2)
    e = np.array(e)
    refraction_point = (e * n1 - n * (
            n1 * np.dot(e, n) - n2 * np.sqrt(under_root))) / n2
    return refraction_point
