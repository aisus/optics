import matplotlib.pyplot as plt

from ray import Ray
from traceObjects.circle import Circle
from traceObjects.ellipse import Ellipse
from traceObjects.plane import Plane


def plane_example(ray, n1, n2):
    plane = Plane(
        normal=[1.0, 0.0, 0.0],
        radius_vec=[3.0, 2.0, 0.0],
    )
    plane_lgnd = plane.draw("plane")
    ray_lgnd, reflection_lgnd, refraction_lgnd = plane.draw_intersection(ray, n1, n2, "plane")
    plt.legend([plane_lgnd, ray_lgnd, reflection_lgnd, refraction_lgnd])
    plt.grid(True)


def circle_example(ray, n1, n2):
    circle = Circle(
        origin=[1, 1, 0],
        radius=2
    )
    circle_lgnd = circle.draw("circle")
    ray_res, norm_1, norm_2 = circle.draw_intersection(ray, n1, n2, "circle")
    plt.legend([circle_lgnd])
    plt.grid(True)


def ellipse_example(ray, n1, n2):
    ellipse = Ellipse(
        origin=[2, 2],
        a=4,
        b=2
    )
    ray.origin = [ray.origin[0], ray.origin[1]]
    ray.direction = [ray.direction[0], ray.direction[1]]
    ellipse_lgnd = ellipse.draw("ellipse")
    lgnd = ellipse.draw_intersection(ray, n1, n2, "ellipse")
    plt.grid(True)


if __name__ == "__main__":
    ray = Ray(
        origin=[-1, -1, 0],
        direction=[1, 2, 0]
    )
    # plane_example(ray, 0.5, 0.5)
    # circle_example(ray, 1.5, 1)
    ellipse_example(ray, 1, 1.5)
    plt.show()
