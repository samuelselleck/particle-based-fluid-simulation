extern crate glutin_window;
extern crate graphics;
extern crate opengl_graphics;
extern crate piston;

use glutin_window::GlutinWindow as Window;
use opengl_graphics::{GlGraphics, OpenGL};
use piston::event_loop::{EventSettings, Events};
use piston::input::{RenderArgs, RenderEvent, UpdateArgs, UpdateEvent};
use piston::window::WindowSettings;

use rand::Rng;
use std::collections::HashMap;

const KERNEL_RADIUS: f64 = 16.0;
const KR_SQUARED: f64 = KERNEL_RADIUS*KERNEL_RADIUS;
const BUCKET_SIZE: f64 = KERNEL_RADIUS;

const WINDOW_WIDTH: f64 = 1000.0;
const WINDOW_HEIGHT: f64 = 800.0;

const HORIZONTAL_BUCKET_COUNT: i32 = 
    (WINDOW_WIDTH as f64/BUCKET_SIZE) as i32 + 1;

pub struct App {
    gl: GlGraphics,
    particles: Vec<Particle>,
    buckets: HashMap<i32, Vec<usize>>,
}

#[derive(Debug)]
pub struct Particle {
    position: Vec2d,
    velocity: Vec2d,
    force: Vec2d,
    density: f64,
    pressure: f64,
}

#[derive(Debug)]
struct Vec2d {
    x: f64,
    y: f64,
}

const OTHER_BUCKET_OFFSETS: [i32; 9] = [
    -HORIZONTAL_BUCKET_COUNT - 1,
    -HORIZONTAL_BUCKET_COUNT,
    -HORIZONTAL_BUCKET_COUNT + 1,
     -1, 0, 1,
     HORIZONTAL_BUCKET_COUNT - 1,
     HORIZONTAL_BUCKET_COUNT,
     HORIZONTAL_BUCKET_COUNT + 1];

struct BucketLocality<> {
    all: Vec<usize>,
    main: Vec<usize>,
    index: usize,
}

impl BucketLocality<> {
    fn new(bucket_index: i32, buckets: &HashMap<i32, Vec<usize>>) -> BucketLocality {

        let mut main = vec![];
        match buckets.get(&bucket_index) {
            Some(points) => for p in points { main.push(*p) },
            None => (),
        }

        let mut all = vec![];
        for offset in OTHER_BUCKET_OFFSETS.iter() {
            match buckets.get(&(bucket_index + offset)) {
                Some(points) => for p in points { all.push(*p) },
                None => (),
            }
        };

        BucketLocality { all, main, index: 0 }
    }
}

impl Iterator for BucketLocality {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        let main_index = self.index / self.all.len();
        if main_index >= self.main.len() { return None }
        let other_index = self.index % self.all.len();
        self.index += 1;
        Some((self.main[main_index], self.all[other_index]))
    }
}

impl Vec2d {

    fn dist_squared(&self, other: &Vec2d) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx*dx + dy*dy
    }

    fn dist(&self, other: &Vec2d) -> f64 {
        self.dist_squared(other).sqrt()
    }

    fn add(&mut self, other: &Vec2d) -> &mut Vec2d {
        self.x += other.x;
        self.y += other.y;
        self
    }

    fn sub(&mut self, other: &Vec2d) -> &mut Vec2d {
        self.x -= other.x;
        self.y -= other.y;
        self
    }

    fn mult(&mut self, factor: f64) -> &mut Vec2d {
        self.x *= factor;
        self.y *= factor;
        self
    }

    fn div(&mut self, factor: f64) -> &mut Vec2d {
        self.x /= factor;
        self.y /= factor;
        self
    }

    fn copy(&self) -> Vec2d {
        Vec2d {x: self.x, y: self.y}
    }
}

impl App {
    fn render(&mut self, args: &RenderArgs) {
        use graphics::*;

        const BACKGROUND: [f32; 4] = [0.9, 0.85, 1.0, 1.0];
        const PARTICLE: [f32; 4] = [0.3, 0.5, 0.9, 1.0];

        let circle = ellipse::circle(0.0, 0.0, 5.0);

        self.gl.draw(args.viewport(), |_c, gl| {
            clear(BACKGROUND, gl);
        });
        for particle in self.particles.iter() {
            let pos = &particle.position;
            self.gl.draw(args.viewport(), |c, gl| {
                let transform = c.transform.trans(pos.x, pos.y);
                ellipse(PARTICLE, circle, transform, gl);
            });
        }
    }

    fn update(&mut self, _args: &UpdateArgs) {

        const PI: f64 = 3.14159265358979323846264338327950288f64;

        // simulation parameters

        const MASS: f64 = 65.0;
        const REST_DENS: f64 = 1000.0;
        const GAS_CONST: f64 = 2000.0;
        
        const VISC: f64 = 250.0;
        
        const G: Vec2d = Vec2d {x: 0.0, y: 12000.0*9.8};

        const DT: f64 = 0.0008;
        const EPS: f64 = KERNEL_RADIUS;
        const BOUND_DAMPING: f64 = -0.5;

        let poly_6: f64 = 315.0/(65.0*PI*KERNEL_RADIUS.powf(9.0));
        let spiky_grad: f64 = -45.0/(PI*KERNEL_RADIUS.powf(6.0));
        let visc_lap: f64 = 45.0/(PI*KERNEL_RADIUS.powf(6.0));

        let p: &mut [Particle] = &mut self.particles;

        //Calculate pressure and density:
        for i in 0..p.len() {
            p[i].density = 0.0;
        }
        
        for bucket_index in self.buckets.keys() {
            let bucket_locality = BucketLocality::new(*bucket_index, &self.buckets);

            for (i, j) in bucket_locality {
                let dist_sq = p[i].position.dist_squared(&p[j].position);
                if dist_sq < KR_SQUARED {
                    p[i].density += MASS*poly_6*(KR_SQUARED - dist_sq).powf(3.0);
                }
            }
        }
        for i in 0..p.len() {
            p[i].pressure = GAS_CONST*(p[i].density - REST_DENS);
            p[i].force = G.copy();
            p[i].force.mult(p[i].density);
        }

        //Calculate forces:
        for bucket_index in self.buckets.keys() {

            let bucket_locality = BucketLocality::new(*bucket_index, &self.buckets);
            for (i, j) in bucket_locality {

                if i == j {
                    continue;
                }

                let distance = p[i].position.dist(&p[j].position);
                let mut dir = p[i].position.copy();
                dir.sub(&p[j].position).div(distance);

                if distance < KERNEL_RADIUS {

                    let factor = 
                        MASS*
                        (p[i].pressure + p[j].pressure)/
                        (2.0*p[j].density)*
                        spiky_grad*(KERNEL_RADIUS - distance).powf(2.0);
                    p[i].force.add(dir.mult(factor));

                    let factor = VISC*MASS/p[j].density*visc_lap*(KERNEL_RADIUS-distance);
                    let mut j_vel_copy = p[j].velocity.copy();
                    p[i].force.add((&mut j_vel_copy).sub(&p[i].velocity).mult(factor));
                }
            }
        }

        //Update positions:
        for i in 0..p.len() {

            let old_bucket = get_pos_hash(&p[i].position);

            p[i].velocity.add(p[i].force.mult(DT/p[i].density));
            p[i].position.add(p[i].velocity.copy().mult(DT));

            if p[i].position.y > WINDOW_HEIGHT - EPS {
                p[i].velocity.y *= BOUND_DAMPING;
                p[i].position.y = WINDOW_HEIGHT - EPS;
            }

            if p[i].position.y < EPS {
                p[i].velocity.y *= BOUND_DAMPING;
                p[i].position.y = EPS;
            }

            if p[i].position.x > WINDOW_WIDTH - EPS {
                p[i].velocity.x *= BOUND_DAMPING;
                p[i].position.x = WINDOW_WIDTH - EPS;
            }

            if p[i].position.x < EPS {
                p[i].velocity.x *= BOUND_DAMPING;
                p[i].position.x = EPS;
            }

            let new_bucket = get_pos_hash(&p[i].position);

            if !(old_bucket == new_bucket) {
                match self.buckets.get_mut(&old_bucket) {
                    Some(points) => {
                        points.retain(|&x| x != i);
                        if points.len() == 0 {
                            self.buckets.remove(&old_bucket);
                        }
                        ()},
                    None => (),
                };
            
                match self.buckets.get_mut(&new_bucket) {
                    Some(points) => {points.push(i); ()},
                    None => {self.buckets.insert(new_bucket, vec![i]); ()},
                };
            }
        }
    }

}

fn get_pos_hash(p: &Vec2d) -> i32 {
    (p.x /BUCKET_SIZE).floor() as i32+
    HORIZONTAL_BUCKET_COUNT*(p.y/BUCKET_SIZE).floor() as i32
}

fn main() {

    let opengl = OpenGL::V3_2;

    let mut window: Window = WindowSettings::new("SPH", [WINDOW_WIDTH, WINDOW_HEIGHT])
        .graphics_api(opengl)
        .exit_on_esc(true)
        .build()
        .unwrap();

    let mut app = App {
        gl: GlGraphics::new(opengl),
        particles: {
            let mut particles = vec![];
            for x in 0..40 {
                for y in 0..40 {

                    let x_jitter = rand::thread_rng()
                        .gen_range(1, 1000) as f64/1000.0;
                    let y_jitter = rand::thread_rng()
                        .gen_range(1, 1000) as f64/1000.0;

                    particles.push(Particle {
                        position: Vec2d {x: x as f64*KERNEL_RADIUS + 100.0 + x_jitter, y: y as f64*KERNEL_RADIUS + 100.0 + y_jitter},
                        velocity: Vec2d {x: 0.0, y: 0.0},
                        force: Vec2d {x: 0.0, y: 0.0},
                        density: 0.0,
                        pressure: 0.0,
                    });
                }
            }
            particles
        },
        buckets: HashMap::new()
    };

    //Init buckets:
    for i in 0..app.particles.len() {
        let bucket = get_pos_hash(&app.particles[i].position);
        match app.buckets.get_mut(&bucket) {
            Some(points) => {points.push(i); ()},
            None => {app.buckets.insert(bucket, vec![i]); ()},
        }
    }

    let mut events = Events::new(EventSettings::new());

    while let Some(e) = events.next(&mut window) {
        if let Some(args) = e.render_args() {
            app.render(&args);
        }

        if let Some(args) = e.update_args() {
            app.update(&args);
        }
    }
}