import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse as Elp

import utils
from ray import Ray


class Ellipse(object):
    def __init__(self, origin, a, b):
        self.origin = np.array(origin)
        self.a = a
        self.b = b

    def norm(self, x0, y0):
        return -(y0 * (self.a ** 2)) / (x0 * (self.b ** 2))

    # def norm2(self, x0, y0):
    #     return -(y0 * (self.a ** 2)) / (x0 * (self.b ** 2))

    def intersection_ellipse_points(self, ray: Ray):
        m = np.array([[self.b, 0], [0, self.a]])
        # mr = np.dot(m, ray.origin - self.origin)
        # print(mr)
        # me = np.dot(m, ray.direction)
        a = np.dot(np.dot(m, ray.direction), np.dot(m, ray.direction))
        b = 2 * np.dot(np.dot(m, ray.direction), np.dot(m, ray.origin - self.origin))
        c = np.dot(np.dot(m, ray.origin - self.origin), np.dot(m, ray.origin - self.origin)) - (m[0][0] * m[1][1]) ** 2

        # print(a, b, c)
        return np.roots([a, b, c])

    # m = np.array([[self.b, 0], [0, self.a]])
    # mr = np.dot(m, ray.origin - self.origin)
    # a = np.dot(np.dot(m, ray.direction), np.dot(m, ray.direction))
    # b = 2 * np.dot(np.dot(m, ray.direction), mr)
    # c = np.dot(mr, mr) - (m[0][0] * m[1][1]) ** 2
    # print(a, b, c)
    # return np.roots([a, b, c])

    def draw(self, figure="Fig"):
        plt.figure(figure)
        ax = plt.gca()
        ellipse = Elp((self.origin[0], self.origin[1]), 2 * self.a, 2 * self.b,
                      fill=False,
                      label="Ellipse")
        ax.add_patch(ellipse)
        # plt.axis('scaled')
        plt.gca().set_aspect("equal")
        plt.grid(True)

    def intersect_ray(self, ray: Ray, nn1, nn2):
        t = self.intersection_ellipse_points(ray)
        r = ray.origin + np.array(ray.direction) * t[1]
        r2 = ray.origin + np.array(ray.direction) * t[0]
        a = self.norm(r[0], r[1])
        a = np.arctan(a)
        a = np.array([np.cos(a), np.sin(a)])
        a2 = - self.norm(r2[0], r2[1])
        a2 = np.arctan(a2)
        a2 = np.array([np.cos(a2), np.sin(a2)])
        refl = - utils.reflection(ray.direction, a)
        refl2 = - utils.reflection(ray.direction, a2)
        refr = - utils.refraction(ray.direction, a, nn1, nn2)
        refr2 = - utils.refraction(ray.direction, a2, nn1, nn2)
        ray_lgnd, = plt.plot([ray.origin[0], r[0], r2[0]], [ray.origin[1], r[1], r2[1]])
        normal1, = plt.plot([r[0], r[0] - a[0]], [r[1], r[1] - a[1]])
        normal2, = plt.plot([r2[0], r2[0] - a2[0]], [r2[1], r2[1] - a2[1]])
        reflection, = plt.plot([r[0], r[0] - refl[0]], [r[1], r[1] - refl[1]])
        reflection2, = plt.plot([r2[0], r2[0] - refl2[0]], [r2[1], r2[1] - refl2[1]])
        refraction, = plt.plot([r[0], r[0] - refr[0]], [r[1], r[1] - refr[1]])
        refraction2, = plt.plot([r2[0], r2[0] - refr2[0]], [r2[1], r2[1] - refr2[1]])
        plt.legend([ray_lgnd, normal1, normal2],
                   ["Ray", "Normal1", "Normal2"])
        # a = np.array([np.cos(a), np.sin(a)])

    def draw_intersection(self, ray: Ray, nn1, nn2, figure="Fig"):
        plt.figure(figure)
        self.intersect_ray(ray, nn1, nn2)
