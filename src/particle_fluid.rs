
use std::collections::HashMap;
use super::vec_math::Vec2d;
use std::f32::consts::PI;
use rand::Rng;

pub struct Parameters {

    kernel_radius: f32,
    container_width: f32,
    container_height: f32,
    particle_mass: f32,
    particle_rest_density: f32,
    gas_const: f32,
    viscosity: f32,
    gravity: Vec2d,
    bound_damping: f32,

    //Derivative values:
    kernel_radius_squared: f32,
    bucket_size: f32,
    horizontal_bucket_count: i32,
    surrounding_bucket_index_offsets: [i32; 8],
    epsilon: f32,
    jitter: f32,

    //Internals:
    poly_6: f32,
    spiky_grad: f32,
    visc_lap: f32,
}

impl Parameters {
    pub fn new(kernel_radius: f32, container_width: f32, container_height: f32) -> Parameters{

        let bucket_size = kernel_radius;
        let kernel_radius_squared = kernel_radius*kernel_radius;
        let horizontal_bucket_count = (container_width/bucket_size) as i32 + 1;

        let surrounding_bucket_index_offsets = [
            -horizontal_bucket_count - 1,
            -horizontal_bucket_count,
            -horizontal_bucket_count + 1,
            -1, 1,
            horizontal_bucket_count - 1,
            horizontal_bucket_count,
            horizontal_bucket_count + 1
        ];

        let epsilon = kernel_radius;
        let jitter = kernel_radius/1e5;

        Parameters {
            kernel_radius,
            container_width,
            container_height,
            kernel_radius_squared,
            bucket_size,
            horizontal_bucket_count,
            surrounding_bucket_index_offsets,
            particle_mass: 65.0,
            particle_rest_density: 1000.0,
            gas_const: 2000.0,
            viscosity: 250.0,
            gravity: Vec2d {x: 0.0, y: 12000.0*9.8},
            epsilon,
            bound_damping: -0.5,
            poly_6: 315.0/(65.0*PI*kernel_radius.powf(9.0)),
            spiky_grad: -45.0/(PI*kernel_radius.powf(6.0)),
            visc_lap: 45.0/(PI*kernel_radius.powf(6.0)),
            jitter,
        }
    }

    fn get_pos_hash(&self, p: &Vec2d) -> i32 {
        (p.x /self.bucket_size).floor() as i32+
        self.horizontal_bucket_count*(p.y/self.bucket_size).floor() as i32
    }
}

pub struct ParticleFluid {
    params: Parameters,
    pub particles: Vec<Particle>,
    buckets: HashMap<i32, Vec<usize>>,
}

impl ParticleFluid {
    pub fn new(params: Parameters, particle_positions: Vec<Vec2d>) -> ParticleFluid {


        let particles: Vec<Particle> = particle_positions.into_iter()
            .map(|p| {

                let jitter = || { rand::thread_rng()
                    .gen_range(0, 1000) as f32*params.jitter };

                let position = Vec2d{x: jitter(), y: jitter()}.add(&p).copy();

                Particle { 
                    position,
                    velocity: Vec2d {x: 0.0, y: 0.0},
                    force: Vec2d {x: 0.0, y: 0.0},
                    density: 0.0,
                    pressure: 0.0,
                }
            }).collect();

        let mut buckets: HashMap<i32, Vec<usize>> = HashMap::new();

        for i in 0..particles.len() {
            let bucket = params.get_pos_hash(&particles[i].position);
            match buckets.get_mut(&bucket) {
                Some(points) => {points.push(i); ()},
                None => {buckets.insert(bucket, vec![i]); ()},
            }
        }

        ParticleFluid {
            params,
            particles,
            buckets,
        }
    }

    pub fn step(&mut self, dt: f32) {

        self.calculate_pressure_and_density();
        self.calculate_forces();
        self.update_positions(dt);
    }

    fn calculate_pressure_and_density(&mut self) {

        let p: &mut [Particle] = &mut self.particles;

        for i in 0..p.len() {
            p[i].density = 0.0;
        }
        
        for bucket_index in self.buckets.keys() {
            let bucket_locality = BucketLocality::new(
                *bucket_index,
                 &self.buckets,
                  self.params.surrounding_bucket_index_offsets
            );

            for (i, j) in bucket_locality {
                let dist_sq = p[i].position.dist_squared(&p[j].position);
                if dist_sq < self.params.kernel_radius_squared {
                    p[i].density += self.params.particle_mass*self.params.poly_6*
                        (self.params.kernel_radius_squared - dist_sq).powf(3.0);
                }
            }
        }

        for i in 0..p.len() {
            p[i].pressure = self.params.gas_const*
                (p[i].density - self.params.particle_rest_density); 
        }
    }

    fn calculate_forces(&mut self) {

        let p: &mut [Particle] = &mut self.particles;

        for i in 0..p.len() {
            p[i].force = self.params.gravity.copy();
            p[i].force.mult(p[i].density);
        }

        for bucket_index in self.buckets.keys() {

            let bucket_locality = BucketLocality::new(
                *bucket_index,
                 &self.buckets,
                  self.params.surrounding_bucket_index_offsets
            );

            for (i, j) in bucket_locality {

                if i == j {
                    continue;
                }

                let distance = p[i].position.dist(&p[j].position);
                let mut dir = p[i].position.copy();
                dir.sub(&p[j].position).div(distance);

                if distance < self.params.kernel_radius {

                    let factor = 
                        self.params.particle_mass*
                        (p[i].pressure + p[j].pressure)/
                        (2.0*p[j].density)*
                        self.params.spiky_grad*(self.params.kernel_radius - distance).powf(2.0);
                    p[i].force.add(dir.mult(factor));

                    let factor = self.params.viscosity*self.params.particle_mass/
                        p[j].density*self.params.visc_lap*(self.params.kernel_radius-distance);
                    let mut j_vel_copy = p[j].velocity.copy();
                    p[i].force.add((&mut j_vel_copy).sub(&p[i].velocity).mult(factor));
                }
            }
        }
    }

    fn update_positions(&mut self, dt: f32) {

        let p: &mut [Particle] = &mut self.particles;

        for i in 0..p.len() {

            let old_bucket = self.params.get_pos_hash(&p[i].position);

            p[i].velocity.add(p[i].force.mult(dt/p[i].density));
            p[i].position.add(p[i].velocity.copy().mult(dt));

            let height = self.params.container_height;
            let width = self.params.container_width;
            let bound_damping = self.params.bound_damping;
            let eps = self.params.epsilon;

            if p[i].position.y > height - eps {
                p[i].velocity.y *= bound_damping;
                p[i].position.y = height - eps;
            }

            if p[i].position.y < eps {
                p[i].velocity.y *= bound_damping;
                p[i].position.y = eps;
            }

            if p[i].position.x > width - eps {
                p[i].velocity.x *= bound_damping;
                p[i].position.x = width - eps;
            }

            if p[i].position.x < eps {
                p[i].velocity.x *= bound_damping;
                p[i].position.x = eps;
            }

            let new_bucket = self.params.get_pos_hash(&p[i].position);

            ParticleFluid::update_bucket(&mut self.buckets, old_bucket, new_bucket, i);
        }
    }

    fn update_bucket(buckets: &mut HashMap<i32, Vec<usize>>, old_bucket: i32, new_bucket: i32, particle_index: usize) {
        if !(old_bucket == new_bucket) {
            match buckets.get_mut(&old_bucket) {
                Some(points) => {
                    points.retain(|&x| x != particle_index);
                    if points.len() == 0 {
                        buckets.remove(&old_bucket);
                    }
                    ()},
                None => (),
            };
        
            match buckets.get_mut(&new_bucket) {
                Some(points) => {points.push(particle_index); ()},
                None => {buckets.insert(new_bucket, vec![particle_index]); ()},
            };
        }
    }
}

#[derive(Debug)]
pub struct Particle {
    pub position: Vec2d,
    velocity: Vec2d,
    force: Vec2d,
    density: f32,
    pub pressure: f32,
}

struct BucketLocality<> {
    all: Vec<usize>,
    main: Vec<usize>,
    index: usize,
}

impl BucketLocality {
    fn new(bucket_index: i32, buckets: &HashMap<i32, Vec<usize>>, surrounding_bucket_index_offsets: [i32; 8]) -> BucketLocality {

        let mut main = vec![];
        match buckets.get(&bucket_index) {
            Some(points) => for p in points { main.push(*p) },
            None => (),
        }

        let mut all = main.clone();
        for offset in surrounding_bucket_index_offsets.iter() {
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