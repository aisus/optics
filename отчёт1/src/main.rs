extern crate gnuplot;
extern crate nalgebra as na;

use gnuplot::*;
use na::*;

use std::f64::consts::*;
use std::io;
use std::io::{BufRead, StdinLock};

fn read_f64(handler: &mut StdinLock) -> f64 {
    let mut buffer = String::new();

    handler
        .read_line(&mut buffer)
        .expect("Error while reading line");

    buffer
        .trim()
        .parse::<f64>()
        .expect("Error while parsing float")
}

fn read_vec2(handler: &mut StdinLock) -> Vector2<f64> {
    let mut buffer = String::new();

    handler
        .read_line(&mut buffer)
        .expect("Error while reading line");

    let splitted = buffer.split(',').collect::<Vec<_>>();

    if splitted.len() != 2 {
        panic!("Origin should have 2 coordinates");
    }

    let x = splitted[0]
        .trim()
        .parse::<f64>()
        .expect("Error while parsing x");
    let y = splitted[1]
        .trim()
        .parse::<f64>()
        .expect("Error while parsing y");

    Vector2::new(x, y)
}

#[derive(Debug, PartialEq, Copy, Clone)]
struct Ray {
    origin: Vector2<f64>,
    direction: Vector2<f64>,
}

#[derive(Debug, PartialEq, Copy, Clone)]
struct Plane {
    rad: Vector2<f64>,
    norm: Vector2<f64>,
}

#[derive(Debug, PartialEq, Copy, Clone)]
struct Circle {
    origin: Vector2<f64>,
    radius: f64,
}

#[derive(Debug, PartialEq, Copy, Clone)]
struct Ellipse {
    origin: Vector2<f64>,
    radiuses: Vector2<f64>,
}

type IntersectPoint = Vector2<f64>;
type ReflectRay = Ray;
type RefractRay = Option<Ray>;

struct RefractionData {
    n1: f64,
    n2: f64,
}

type IntersectData = (IntersectPoint, ReflectRay, RefractRay);

trait Intersect {
    fn intersect(&self, ray: &Ray, refract_data: &RefractionData) -> Option<IntersectData>;
}

impl Intersect for Plane {
    fn intersect(&self, ray: &Ray, refract_data: &RefractionData) -> Option<IntersectData> {
        let n1 = refract_data.n1;
        let n2 = refract_data.n2;

        let dir = self.rad - ray.origin;
        let proj = self.norm.dot(&ray.direction);

        if proj.abs() <= 1e-8 {
            return None;
        }

        let t = self.norm.dot(&dir) / proj;

        let intersection = ray.origin + ray.direction * t;

        let reflected_ray = Ray {
            origin: intersection,
            direction: (ray.direction - 2.0 * ray.direction.dot(&self.norm) * self.norm)
                .normalize(),
        };

        let root = 1.0 - (n1 / n2).powi(2) * (1.0 - ray.direction.dot(&self.norm).powi(2));

        let refracted_ray = if root <= 0.0 {
            None
        } else {
            let g = n1 * ray.direction.dot(&self.norm) - n2 * root.sqrt();
            let refracted_ray = Ray {
                origin: intersection,
                direction: ((n1 * ray.direction - self.norm * g) / n2).normalize(),
            };

            let alpha1 = ray.direction.dot(&self.norm).acos();
            let alpha2 = refracted_ray.direction.dot(&self.norm).acos();

            dbg!(alpha1.sin() / alpha2.sin());
            dbg!(n2 / n1);

            Some(refracted_ray)
        };

        Some((intersection, reflected_ray, refracted_ray))
    }
}

impl Intersect for Circle {
    fn intersect(&self, ray: &Ray, refract_data: &RefractionData) -> Option<IntersectData> {
        let mut n1 = refract_data.n1;
        let mut n2 = refract_data.n2;

        let is_inside = (ray.origin - self.origin).magnitude() < self.radius;

        let dir = ray.origin - self.origin;

        let root = dir.dot(&ray.direction).powi(2) - dir.dot(&dir) + self.radius.powi(2);

        if root + 1e-9 <= 0. {
            return None;
        }

        let has_intersection =
            (dir.dot(&ray.direction).powi(2) - (dir.dot(&dir) - self.radius.powi(2))) >= 0.;

        if !has_intersection {
            return None;
        }

        let t1 = -dir.dot(&ray.direction) - root.sqrt();
        let t2 = -dir.dot(&ray.direction) + root.sqrt();

        let intersections = vec![
            ray.origin + t1 * ray.direction,
            ray.origin + t2 * ray.direction,
        ];

        let intersection = intersections
            .iter()
            .filter(|&&intersection| (intersection - ray.origin).dot(&ray.direction) >= 0.0)
            .min_by(|&&intersection1, &&intersection2| {
                let first_magnitude = (intersection1 - ray.origin).magnitude();
                let second_magnitude = (intersection2 - ray.origin).magnitude();

                first_magnitude.partial_cmp(&second_magnitude).unwrap()
            });

        dbg!(&intersection);

        let intersection = *(intersection?);

        let mut norm_vec = -(intersection - self.origin).normalize();

        if is_inside {
            norm_vec = -norm_vec;
            std::mem::swap(&mut n1, &mut n2);
        }

        let reflected_ray = Ray {
            origin: intersection,
            direction: (ray.direction - 2.0 * norm_vec * ray.direction.dot(&norm_vec)).normalize(),
        };

        let root = 1.0 - (n1 / n2).powi(2) * (1.0 - ray.direction.dot(&norm_vec).powi(2));

        let refracted_ray = if root <= 0.0 {
            None
        } else {
            let g = n1 * ray.direction.dot(&norm_vec) - n2 * root.sqrt();
            let refracted_ray = Ray {
                origin: intersection,
                direction: ((n1 * ray.direction - norm_vec * g) / n2).normalize(),
            };

            let alpha1 = ray.direction.dot(&norm_vec).acos();
            let alpha2 = refracted_ray.direction.dot(&norm_vec).acos();

            dbg!(alpha1.sin() / alpha2.sin());
            dbg!(n2 / n1);

            Some(refracted_ray)
        };

        Some((intersection, reflected_ray, refracted_ray))
    }
}

fn power(a: f64, b: i32) -> f64 {
    a.powi(b)
}

fn sqrt(a: f64) -> f64 {
    a.sqrt()
}

impl Intersect for Ellipse {
    fn intersect(&self, ray: &Ray, refract_data: &RefractionData) -> Option<IntersectData> {
        let mut n1 = refract_data.n1;
        let mut n2 = refract_data.n2;

        let is_inside = ((ray.origin.x - self.origin.x).powi(2) / self.radiuses.x.powi(2)
            + (ray.origin.y - self.origin.y).powi(2) / self.radiuses.y.powi(2))
            <= 1.0;

        let ex = ray.direction.x;
        let ey = ray.direction.y;

        let px = self.origin.x;
        let py = self.origin.y;

        let x0 = ray.origin.x;
        let y0 = ray.origin.y;

        let a = self.radiuses.x;
        let b = self.radiuses.y;

        let root = power(b, 2) * power(ex, 2) + power(a, 2) * power(ey, 2)
            - power(ey * (px - x0) + ex * (-py + y0), 2);

        if root + 1e-9 <= 0.0 {
            return None;
        }

        let t1: f64 = (power(b, 2) * ex * px + power(a, 2) * ey * py
            - power(b, 2) * ex * x0
            - power(a, 2) * ey * y0
            - a * b * sqrt(root))
            / (power(b, 2) * power(ex, 2) + power(a, 2) * power(ey, 2));
        let t2: f64 = (power(b, 2) * ex * px + power(a, 2) * ey * py
            - power(b, 2) * ex * x0
            - power(a, 2) * ey * y0
            + a * b * sqrt(root))
            / (power(b, 2) * power(ex, 2) + power(a, 2) * power(ey, 2));

        let intersections: Vec<Vector2<f64>> = vec![
            ray.origin + t1 * ray.direction,
            ray.origin + t2 * ray.direction,
        ];

        let intersection = intersections
            .iter()
            .filter(|&&intersection| (intersection - ray.origin).dot(&ray.direction) >= 0.0)
            .min_by(|&&intersection1, &&intersection2| {
                let first_magnitude = (intersection1 - ray.origin).magnitude();
                let second_magnitude = (intersection2 - ray.origin).magnitude();

                first_magnitude.partial_cmp(&second_magnitude).unwrap()
            });

        dbg!(&intersection);

        let intersection = *(intersection?);

        let x = intersection.x;
        let y = intersection.y;

        let mut norm_vec: Vector2<f64> = Vector2::new(
            (2.0 * (-px + x))
                / (power(a, 2)
                    * sqrt(
                        (4.0 * power(px - x, 2)) / power(a, 4) + (4.0 * power(py - y, 2)) / power(b, 4),
                    )),
            (2.0 * (-py + y))
                / (power(b, 2)
                    * sqrt(
                        (4.0 * power(px - x, 2)) / power(a, 4) + (4.0 * power(py - y, 2)) / power(b, 4),
                    )),
        );

        norm_vec = -norm_vec.normalize();

        if is_inside {
            norm_vec = -norm_vec;
            std::mem::swap(&mut n1, &mut n2);
        }

        let reflected_ray = Ray {
            origin: intersection,
            direction: (ray.direction - 2.0 * norm_vec * ray.direction.dot(&norm_vec)).normalize(),
        };

        let root = 1.0 - (n1 / n2).powi(2) * (1.0 - ray.direction.dot(&norm_vec).powi(2));

        let refracted_ray = if root <= 0.0 {
            None
        } else {
            let g = n1 * ray.direction.dot(&norm_vec) - n2 * root.sqrt();
            let refracted_ray = Ray {
                origin: intersection,
                direction: ((n1 * ray.direction - norm_vec * g) / n2).normalize(),
            };

            let alpha1 = ray.direction.dot(&norm_vec).acos();
            let alpha2 = refracted_ray.direction.dot(&norm_vec).acos();

            dbg!(alpha1.sin() / alpha2.sin());
            dbg!(n2 / n1);

            Some(refracted_ray)
        };

        Some((intersection, reflected_ray, refracted_ray))
    }
}

fn plane_solution(mut handler: &mut std::io::StdinLock) {
    println!("Enter ray origin point (x, y):");
    let origin = read_vec2(&mut handler);
    println!("Enter ray direction (x, y):");
    let direction = read_vec2(&mut handler);

    if direction.magnitude() <= 1e-5 {
        panic!("Bad direction vector");
    }

    let direction = direction.normalize();

    let ray = Ray { origin, direction };

    println!("{:#?}", ray);

    println!("Enter plane rad. vector (x, y):");
    let rad = read_vec2(&mut handler);
    println!("Enter plane norm. vector (x, y):");
    let norm = read_vec2(&mut handler);

    if norm.magnitude() <= 1e-5 {
        panic!("Bad norm vector");
    }

    let norm = norm.normalize();

    let plane = Plane { rad, norm };

    println!("{:#?}", plane);

    println!("Enter refraction coeff. before plane:");
    let before = read_f64(&mut handler);
    println!("Enter refraction coeff. after plane:");
    let after = read_f64(&mut handler);

    let refract_data = RefractionData {
        n1: before,
        n2: after,
    };

    let intersect_data = plane.intersect(&ray, &refract_data);

    if let Some((intersection, reflected_ray, refracted_ray)) = intersect_data {
        let plane_ray = Ray {
            origin: plane.rad,
            direction: Vector2::new(plane.norm.y, -plane.norm.x).normalize(),
        };

        dbg!(&intersection);
        dbg!(&reflected_ray);

        let mut fg = Figure::new();
        let axes = fg
            .axes2d()
            .lines(
                &[
                    plane_ray.origin.x - plane_ray.direction.x,
                    plane_ray.origin.x + plane_ray.direction.x,
                ],
                &[
                    plane_ray.origin.y - plane_ray.direction.y,
                    plane_ray.origin.y + plane_ray.direction.y,
                ],
                &[Caption("Plane"), Color("red")],
            )
            .lines(
                &[ray.origin.x, intersection.x],
                &[ray.origin.y, intersection.y],
                &[Caption("Ray"), Color("black")],
            )
            .lines(
                &[
                    reflected_ray.origin.x,
                    reflected_ray.origin.x + reflected_ray.direction.x,
                ],
                &[
                    reflected_ray.origin.y,
                    reflected_ray.origin.y + reflected_ray.direction.y,
                ],
                &[Caption("Reflected ray"), Color("green")],
            );

        if let Some(refracted_ray) = refracted_ray {
            dbg!(&refracted_ray);

            axes.lines(
                &[
                    refracted_ray.origin.x,
                    refracted_ray.origin.x + refracted_ray.direction.x,
                ],
                &[
                    refracted_ray.origin.y,
                    refracted_ray.origin.y + refracted_ray.direction.y,
                ],
                &[Caption("Refracted ray"), Color("blue")],
            );
        }

        fg.show();
    }
}

fn circle_solution(mut handler: &mut std::io::StdinLock) {
    println!("Enter ray origin point (x, y):");
    let origin = read_vec2(&mut handler);
    println!("Enter ray direction (x, y):");
    let direction = read_vec2(&mut handler);

    if direction.magnitude() <= 1e-5 {
        panic!("Bad direction vector");
    }

    let direction = direction.normalize();

    let ray = Ray { origin, direction };

    println!("{:#?}", ray);

    println!("Enter cirlce origin point (x, y):");
    let circle_origin = read_vec2(&mut handler);
    println!("Enter circle radius (r):");
    let circle_radius = read_f64(&mut handler);

    if circle_radius - 1e-9 <= 0. {
        panic!("Circle radius should be positive.");
    }

    let circle = Circle {
        origin: circle_origin,
        radius: circle_radius,
    };

    println!("{:#?}", circle);

    println!("Enter refraction coeff. outside circle:");
    let outside = read_f64(&mut handler);
    println!("Enter refraction coeff. inside circle:");
    let inside = read_f64(&mut handler);

    let refract_data = RefractionData {
        n1: outside,
        n2: inside,
    };

    let intersect_data = circle.intersect(&ray, &refract_data);

    if let Some((intersection, reflected_ray, refracted_ray)) = intersect_data {
        let mut fg = Figure::new();
        let mut points = vec![];

        let steps = 1000;

        let step = PI / f64::from(steps);

        for i in 0..steps {
            let angle = PI * f64::from(i) * step;
            let x = circle.origin.x + angle.cos() * circle.radius;
            let y = circle.origin.y + angle.sin() * circle.radius;
            points.push((x, y));
        }

        let xs: Vec<_> = points.iter().cloned().map(|(x, _)| x).collect();
        let ys: Vec<_> = points.iter().cloned().map(|(_, y)| y).collect();

        dbg!(&intersection);
        dbg!(&reflected_ray);

        let axes = fg
            .axes2d()
            .lines(&xs, &ys, &[Caption("Circle"), Color("red")])
            .lines(
                &[ray.origin.x, intersection.x],
                &[ray.origin.y, intersection.y],
                &[Caption("Ray"), Color("black")],
            )
            .lines(
                &[
                    reflected_ray.origin.x,
                    reflected_ray.origin.x + reflected_ray.direction.x,
                ],
                &[
                    reflected_ray.origin.y,
                    reflected_ray.origin.y + reflected_ray.direction.y,
                ],
                &[Caption("Reflected ray"), Color("green")],
            );

        if let Some(refracted_ray) = refracted_ray {
            dbg!(&refracted_ray);

            axes.lines(
                &[
                    refracted_ray.origin.x,
                    refracted_ray.origin.x + refracted_ray.direction.x,
                ],
                &[
                    refracted_ray.origin.y,
                    refracted_ray.origin.y + refracted_ray.direction.y,
                ],
                &[Caption("Refracted ray"), Color("blue")],
            );
        }

        fg.show();
    }
}

fn ellipse_solution(mut handler: &mut std::io::StdinLock) {
    println!("Enter ray origin point (x, y):");
    let origin = read_vec2(&mut handler);
    println!("Enter ray direction (x, y):");
    let direction = read_vec2(&mut handler);

    if direction.magnitude() <= 1e-5 {
        panic!("Bad direction vector");
    }

    let direction = direction.normalize();

    let ray = Ray { origin, direction };

    println!("{:#?}", ray);

    println!("Enter ellipse origin point (x, y):");
    let ellipse_origin = read_vec2(&mut handler);
    println!("Enter ellipse radiuses (a, b):");
    let ellipse_radiuses = read_vec2(&mut handler);

    if ellipse_radiuses.min() - 1e-9 <= 0. {
        panic!("Ellipse radius should be positive.");
    }

    let ellipse = Ellipse {
        origin: ellipse_origin,
        radiuses: ellipse_radiuses,
    };

    println!("{:#?}", ellipse);

    println!("Enter refraction coeff. outside ellipse:");
    let outside = read_f64(&mut handler);
    println!("Enter refraction coeff. inside ellipse:");
    let inside = read_f64(&mut handler);

    let refract_data = RefractionData {
        n1: outside,
        n2: inside,
    };

    let intersect_data = ellipse.intersect(&ray, &refract_data);

    if let Some((intersection, reflected_ray, refracted_ray)) = intersect_data {
        let mut fg = Figure::new();
        let mut points = vec![];

        let steps = 1000;

        let step = 2.0 * PI / f64::from(steps);

        for i in 0..steps {
            let angle = 2.0 * PI * f64::from(i) * step;
            let x = ellipse.origin.x + angle.cos() * ellipse.radiuses.x;
            let y = ellipse.origin.y + angle.sin() * ellipse.radiuses.y;
            points.push((x, y));
        }

        let xs: Vec<_> = points.iter().cloned().map(|(x, _)| x).collect();
        let ys: Vec<_> = points.iter().cloned().map(|(_, y)| y).collect();

        dbg!(&intersection);
        dbg!(&reflected_ray);

        let axes = fg
            .axes2d()
            .lines(&xs, &ys, &[Caption("Ellipse"), Color("red")])
            .lines(
                &[ray.origin.x, intersection.x],
                &[ray.origin.y, intersection.y],
                &[Caption("Ray"), Color("black")],
            )
            .lines(
                &[
                    reflected_ray.origin.x,
                    reflected_ray.origin.x + reflected_ray.direction.x,
                ],
                &[
                    reflected_ray.origin.y,
                    reflected_ray.origin.y + reflected_ray.direction.y,
                ],
                &[Caption("Reflected ray"), Color("green")],
            );

        if let Some(refracted_ray) = refracted_ray {
            dbg!(&refracted_ray);

            axes.lines(
                &[
                    refracted_ray.origin.x,
                    refracted_ray.origin.x + refracted_ray.direction.x,
                ],
                &[
                    refracted_ray.origin.y,
                    refracted_ray.origin.y + refracted_ray.direction.y,
                ],
                &[Caption("Refracted ray"), Color("blue")],
            );
        }

        fg.show();
    }
}

fn main() {
    let stdin = io::stdin();
    let mut handler = stdin.lock();
    let mut buffer = String::new();

    println!("Enter type of surface (plane, circle, ellipse):");

    handler.read_line(&mut buffer).unwrap();

    let inpt = buffer.trim();
    let inpt = inpt.to_lowercase();
    buffer.clear();

    match inpt.as_ref() {
        "plane" => plane_solution(&mut handler),
        "circle" => circle_solution(&mut handler),
        "ellipse" => ellipse_solution(&mut handler),
        _ => unimplemented!(),
    }
}
