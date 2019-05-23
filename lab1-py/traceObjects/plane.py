import matplotlib.pyplot as plt
import numpy as np

from ray import Ray
from utils import normalize


class Plane(object):
    def __init__(self, normal, radius_vec):
        self.normal = normalize(np.array(normal))
        self.radius_vec = np.array(radius_vec)
        self.normal = np.array(self.norm_v(normal))
        print(self.get_len(normal))

        #  модуль вектора

    def get_len(self, v):
        result = 0
        for i in range(len(v)):
            result += v[i] ** 2
        return np.sqrt(result)

    # нормировка
    def norm_v(self, v):
        l = self.get_len(v)
        for i in range(len(v)):
            v[i] = v[i] / l
        return v


    def draw(self, figure="Fig"):
        plt.figure(figure)
        n_point = np.array([self.normal[1], -self.normal[0], 0.0])
        first = self.radius_vec - n_point * 10
        second = self.radius_vec + n_point * 10
        xs = np.array([first[0], second[0]])
        ys = np.array([first[1], second[1]])
        plane_lgnd, = plt.plot(xs, ys, label="Plane")
        return plane_lgnd

    def intersect_ray(self, ray: Ray, n1, n2):
        if abs(self.normal.dot(ray.direction)) <= 0:
            return None
        t = self.normal.dot(self.radius_vec - ray.origin) / self.normal.dot(ray.direction)
        intersection_point = ray.origin + ray.direction * t
        reflected_direction = ray.direction - 2 * ray.direction.dot(self.normal) * self.normal
        reflected_point = intersection_point + t * reflected_direction
        refracted_direction = n1 * ray.direction - self.normal * (
                n1 * (ray.direction.dot(self.normal)) - n2 *
                np.sqrt(1 - ((n1 ** 2) / (n2 ** 2)) * (
                        1 - (ray.direction.dot(self.normal) ** 2))))
        refracted_direction = normalize(refracted_direction)
        refracted_point = intersection_point + refracted_direction * t
        return intersection_point, reflected_point, refracted_point

    def draw_intersection(self, ray, n1, n2, figure="Fig"):
        intersect = self.intersect_ray(ray, n1, n2)
        if intersect is not None:
            (intersection, reflection, refraction) = intersect
            plt.figure(figure)
            ray_fig, = plt.plot([ray.origin[0], intersection[0]], [ray.origin[1], intersection[1]], label="Ray")
            reflection_fig, = plt.plot([intersection[0], reflection[0]], [intersection[1], reflection[1]],
                                       label="Reflected ray")
            refraction_fig, = plt.plot([intersection[0], refraction[0]], [intersection[1], refraction[1]],
                                       label="Refracted ray")
            print("Plane int. point: {}".format(intersection))
            print(f"len: {self.get_len(intersection)}")
            return ray_fig, reflection_fig, refraction_fig
        else:
            print("No intersection found")
