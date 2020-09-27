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

mod vec_math;
mod particle_fluid;

const WINDOW_WIDTH: f32 = 1280.0;
const WINDOW_HEIGHT: f32 = 720.0;

pub struct App {
    gl: GlGraphics,
    particle_fluid: particle_fluid::ParticleFluid,
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
        for particle in self.particle_fluid.particles.iter() {
            let pos = &particle.position;
            self.gl.draw(args.viewport(), |c, gl| {
                let transform = c.transform.trans(pos.x as f64, pos.y as f64);
                ellipse(PARTICLE, circle, transform, gl);
            });
        }
    }

    fn update(&mut self, _args: &UpdateArgs) {
        const DT: f32 = 0.0008;
        self.particle_fluid.step(DT);
    }

}

fn main() {

    let opengl = OpenGL::V3_2;

    let mut window: Window = WindowSettings::new("SPH", [WINDOW_WIDTH as f64, WINDOW_HEIGHT as f64])
        .graphics_api(opengl)
        .exit_on_esc(true)
        .build()
        .unwrap();

    let mut app = App {
        gl: GlGraphics::new(opengl),
        particle_fluid: {

            let kernel_size: f32 = 16.0;
            let mut particles = vec![];
            for x in 0..70 {
                for y in 0..20 {

                    let x_jitter = rand::thread_rng()
                        .gen_range(1, 1000) as f32/1000.0;
                    let y_jitter = rand::thread_rng()
                        .gen_range(1, 1000) as f32/1000.0;

                    particles.push(vec_math::Vec2d {x: x as f32*kernel_size + 100.0 + x_jitter, y: y as f32*kernel_size + 300.0 + y_jitter});
                }
            }
            let params = particle_fluid::Parameters::new(kernel_size, WINDOW_WIDTH, WINDOW_HEIGHT);

            particle_fluid::ParticleFluid::new(params, particles)
        },
    };
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