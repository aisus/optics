import numpy as np

from utils import normalize


class Ray(object):
    def __init__(self, origin: np.array, direction: np.array):
        self.origin = np.array(origin)
        self.direction = normalize(np.array(direction))
