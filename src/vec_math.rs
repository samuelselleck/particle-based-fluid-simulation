
#[derive(Debug)]
pub struct Vec2d {
    pub x: f32,
    pub y: f32,
}

impl Vec2d {

    pub fn dist_squared(&self, other: &Vec2d) -> f32 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx*dx + dy*dy
    }

    pub fn dist(&self, other: &Vec2d) -> f32 {
        self.dist_squared(other).sqrt()
    }

    pub fn add(&mut self, other: &Vec2d) -> &mut Vec2d {
        self.x += other.x;
        self.y += other.y;
        self
    }

    pub fn sub(&mut self, other: &Vec2d) -> &mut Vec2d {
        self.x -= other.x;
        self.y -= other.y;
        self
    }

    pub fn mult(&mut self, factor: f32) -> &mut Vec2d {
        self.x *= factor;
        self.y *= factor;
        self
    }

    pub fn div(&mut self, factor: f32) -> &mut Vec2d {
        self.x /= factor;
        self.y /= factor;
        self
    }

    pub fn copy(&self) -> Vec2d {
        Vec2d {x: self.x, y: self.y}
    }
}