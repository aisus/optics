import logging

import matplotlib.pyplot as plt
import numpy as np

from ray import Ray
from utils import refraction, reflection

logger = logging.getLogger("CIRCLE")


class Circle(object):
    def __init__(self, origin, radius):
        self.origin = np.array(origin)
        self.radius = radius

    def draw(self, figure="Fig"):
        plt.figure(figure)
        ax = plt.gca()

        circle = plt.Circle(xy=(self.origin[0], self.origin[1]), radius=self.radius, fill=False, label="Circle")
        ax.add_patch(circle)
        plt.gca().set_aspect("equal")
        return circle

    def sphere_intersection(self, ray: Ray):
        t1 = (np.dot(ray.origin - self.origin, ray.direction) + np.sqrt(
            (np.dot(ray.origin - self.origin, ray.direction)) ** 2 - np.dot(ray.direction, ray.direction) * (
                    np.dot(ray.origin - self.origin, ray.origin - self.origin) - self.radius ** 2))) / np.dot(
            ray.direction,
            ray.direction)
        t2 = (np.dot(ray.origin - self.origin, ray.direction) - np.sqrt(
            (np.dot(ray.origin - self.origin, ray.direction)) ** 2 - np.dot(ray.direction, ray.direction) * (
                    np.dot(ray.origin - self.origin, ray.origin - self.origin) - self.radius ** 2))) / np.dot(
            ray.direction,
            ray.direction)
        return t1, t2

    def intersect_ray(self, ray: Ray, en1, en2):
        t1, t2 = self.sphere_intersection(ray)
        r1 = ray.origin + ray.direction * -t1
        r2 = ray.origin + ray.direction * -t2
        n1 = -(r1 - self.origin) / np.sqrt(np.dot(r1 - self.origin, r1 - self.origin))
        n2 = (r2 - self.origin) / np.sqrt(np.dot(r2 - self.origin, r2 - self.origin))
        rfl1 = -reflection(ray.direction, n1)
        rfl2 = -reflection(ray.direction, n2)
        rfr1 = -refraction(ray.direction, n1, en1, en2)
        rfr2 = -refraction(ray.direction, n2, en2, en1)
        print(f"Length: {t1} ; {t2}")
        return r1, r2, rfl1, rfl2, rfr1, rfr2, n1, n2

    def draw_intersection(self, ray, on1, on2, figure="Fig"):
        plt.figure(figure)
        r1, r2, rfl1, rfl2, rfr1, rfr2, n1, n2 = self.intersect_ray(ray, on1, on2)
        ray_plt = plt.plot([ray.origin[0], r1[0]], [ray.origin[1], r1[1]], label="Ray")
        plt.plot([r1[0], r1[0] - n1[0]],
                 [r1[1], r1[1] - n1[1]], label="Ray")
        plt.plot([ray.origin[0], r2[0]],
                 [ray.origin[1], r2[1]], label="Ray")
        plt.plot([r2[0], r2[0] - n2[0]],
                 [r2[1], r2[1] - n2[1]], label="Ray")
        n1_plt, = plt.plot([r1[0], r1[0] - rfl1[0]],
                         [r1[1], r1[1] - rfl1[1]], label="Normal")
        plt.plot([r2[0], r2[0] - rfl2[0]],
                 [r2[1], r2[1] - rfl2[1]], label="Ray")
        n2_plt, = plt.plot([r1[0], r1[0] - rfr1[0]], [r1[1], r1[1] - rfr1[1]], label="Normal")
        plt.plot([r2[0], r2[0] - rfr2[0]],
                 [r2[1], r2[1] - rfr2[1]], label="Ray")
        plt.plot()
        print(f"Sphere int. point: {r1} ; {r2}")
        return ray_plt, n1_plt, n2_plt
