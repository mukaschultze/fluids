extern crate png;

use std::f64;
use std::fs::File;
use std::io::BufWriter;
use std::mem;
use std::path::Path;
use std::ptr::copy_nonoverlapping;
use std::ptr::write_bytes;

fn length(x: f64, y: f64) -> f64 {
    f64::sqrt(x * x + y * y)
}

fn cubic_pulse(mut x: f64) -> f64 {
    x = f64::min(f64::abs(x), 1.0);
    1.0 - x * x * (3.0 - 2.0 * x)
}

fn rotate(x: &mut f64, y: &mut f64, phi: f64) {
    let tmp_x = *x;
    let tmp_y = *y;
    *x = f64::cos(phi) * tmp_x + f64::sin(phi) * tmp_y;
    *y = -f64::sin(phi) * tmp_y + f64::cos(phi) * tmp_y;
}

fn triangle_occupancy(out1: f64, in1: f64, out2: f64) -> f64 {
    0.5 * in1 * in1 / ((out1 - in1) * (out2 - in1))
}

fn trapezoid_occupancy(out1: f64, out2: f64, in1: f64, in2: f64) -> f64 {
    0.5 * (-in1 / (out1 - in1) - in2 / (out2 - in2))
}

fn occupancy(d11: f64, d12: f64, d21: f64, d22: f64) -> f64 {
    let ds = [d11, d12, d22, d21];
    let mut b = 0u8;

    for i in 3..=0 {
        b = (b << 1) | (if ds[i] < 0.0 { 1 } else { 0 });
    }

    match b {
        0x0 => 0.0,

        0x1 => triangle_occupancy(d21, d11, d12),
        0x2 => triangle_occupancy(d11, d12, d22),
        0x4 => triangle_occupancy(d12, d22, d21),
        0x8 => triangle_occupancy(d22, d21, d11),

        0xE => 1.0 - triangle_occupancy(-d21, -d11, -d12),
        0xD => 1.0 - triangle_occupancy(-d11, -d12, -d22),
        0xB => 1.0 - triangle_occupancy(-d12, -d22, -d21),
        0x7 => 1.0 - triangle_occupancy(-d22, -d21, -d11),
        0x3 => trapezoid_occupancy(d21, d22, d11, d12),
        0x6 => trapezoid_occupancy(d11, d21, d12, d22),
        0x9 => trapezoid_occupancy(d12, d22, d11, d21),
        0xC => trapezoid_occupancy(d11, d12, d21, d22),

        0x5 => triangle_occupancy(d11, d12, d22) + triangle_occupancy(d22, d21, d11),
        0xA => triangle_occupancy(d21, d11, d12) + triangle_occupancy(d12, d22, d21),

        0xF => 1.0,
        _ => 0.0,
    }
}

const CELL_FLUID: u8 = 0;
const CELL_SOLID: u8 = 1;

struct SolidBodyFields {
    pub pos_x: f64,
    pub pos_y: f64,
    pub scale_x: f64,
    pub scale_y: f64,
    pub theta: f64,
    pub vel_x: f64,
    pub vel_y: f64,
    pub vel_theta: f64,
}

impl SolidBodyFields {
    pub fn new(
        pos_x: f64,
        pos_y: f64,
        scale_x: f64,
        scale_y: f64,
        theta: f64,
        vel_x: f64,
        vel_y: f64,
        vel_theta: f64,
    ) -> Self {
        Self {
            pos_x,
            pos_y,
            scale_x,
            scale_y,
            theta,
            vel_x,
            vel_y,
            vel_theta,
        }
    }
}

trait SolidBody {
    fn get_fields(&self) -> &SolidBodyFields;
    fn get_fields_mut(&mut self) -> &mut SolidBodyFields;

    // Private (or kind of)
    fn global_to_local(&self, x: &mut f64, y: &mut f64) {
        *x -= self.get_fields().pos_x;
        *y -= self.get_fields().pos_y;
        rotate(x, y, -self.get_fields().theta);
        *x /= self.get_fields().scale_x;
        *y /= self.get_fields().scale_y;
    }

    fn local_to_global(&self, x: &mut f64, y: &mut f64) {
        *x *= self.get_fields().scale_x;
        *y *= self.get_fields().scale_y;
        rotate(x, y, self.get_fields().theta);
        *x += self.get_fields().pos_x;
        *y += self.get_fields().pos_y;
    }

    // Public
    fn distance(&self, x: f64, y: f64) -> f64;
    fn closest_surface_point(&self, x: &mut f64, y: &mut f64);
    fn distance_normal(&self, nx: &mut f64, ny: &mut f64, x: f64, y: f64);

    fn velocity_x(&self, _x: f64, y: f64) -> f64 {
        (self.get_fields().pos_y - y) * self.get_fields().vel_theta + self.get_fields().vel_x
    }

    fn velocity_y(&self, x: f64, _y: f64) -> f64 {
        (x - self.get_fields().pos_x) * self.get_fields().vel_theta + self.get_fields().vel_y
    }

    fn velocity(&self, vx: &mut f64, vy: &mut f64, x: f64, y: f64) {
        *vx = self.velocity_x(x, y);
        *vy = self.velocity_y(x, y);
    }

    fn update(&mut self, timestep: f64) {
        self.get_fields_mut().pos_x += self.get_fields().vel_x * timestep;
        self.get_fields_mut().pos_y += self.get_fields().vel_y * timestep;
        self.get_fields_mut().theta += self.get_fields().vel_theta * timestep;
    }
}

struct SolidBox {
    fields: SolidBodyFields,
}

impl SolidBox {
    pub fn new(fields: SolidBodyFields) -> Self {
        Self { fields }
    }
}

impl SolidBody for SolidBox {
    fn get_fields(&self) -> &SolidBodyFields {
        &self.fields
    }
    fn get_fields_mut(&mut self) -> &mut SolidBodyFields {
        &mut self.fields
    }

    fn distance(&self, mut x: f64, mut y: f64) -> f64 {
        x -= self.get_fields().pos_x;
        y -= self.get_fields().pos_y;
        rotate(&mut x, &mut y, -self.get_fields().theta);
        let dx = f64::abs(x) - self.get_fields().scale_x * 0.5;
        let dy = f64::abs(y) - self.get_fields().scale_y * 0.5;

        if dx >= 0.0 || dy >= 0.0 {
            length(f64::max(dx, 0.0), f64::max(dy, 0.0))
        } else {
            f64::max(dx, dy)
        }
    }

    fn closest_surface_point(&self, x: &mut f64, y: &mut f64) {
        *x -= self.get_fields().pos_x;
        *y -= self.get_fields().pos_y;
        rotate(x, y, -self.get_fields().theta);
        let dx = f64::abs(*x) - self.get_fields().scale_x * 0.5;
        let dy = f64::abs(*y) - self.get_fields().scale_y * 0.5;

        if dx > dy {
            *x = f64::signum(*x) * 0.5 * self.get_fields().scale_x;
        } else {
            *y = f64::signum(*y) * 0.5 * self.get_fields().scale_y;
        }

        rotate(x, y, self.get_fields().theta);
        *x += self.get_fields().pos_x;
        *y += self.get_fields().pos_y;
    }

    fn distance_normal(&self, nx: &mut f64, ny: &mut f64, mut x: f64, mut y: f64) {
        x -= self.get_fields().pos_x;
        y -= self.get_fields().pos_y;
        rotate(&mut x, &mut y, -self.get_fields().theta);

        if f64::abs(x) - self.get_fields().scale_x * 0.5
            > f64::abs(y) - self.get_fields().scale_y * 0.5
        {
            *nx = f64::signum(x);
            *ny = 0.0;
        } else {
            *nx = 0.0;
            *ny = f64::signum(y);
        }

        rotate(nx, ny, self.get_fields().theta);
    }
}

struct SolidSphere {
    fields: SolidBodyFields,
}

impl SolidSphere {
    pub fn new(fields: SolidBodyFields) -> Self {
        Self { fields }
    }
}

impl SolidBody for SolidSphere {
    fn get_fields(&self) -> &SolidBodyFields {
        &self.fields
    }
    fn get_fields_mut(&mut self) -> &mut SolidBodyFields {
        &mut self.fields
    }

    fn distance(&self, x: f64, y: f64) -> f64 {
        length(x - self.get_fields().pos_x, y - self.get_fields().pos_y)
            - self.get_fields().scale_x * 0.5
    }

    fn closest_surface_point(&self, x: &mut f64, y: &mut f64) {
        self.global_to_local(x, y);
        let r = length(*x, *y);
        if r < 1e-4f64 {
            *x = 0.5;
            *y = 0.0;
        } else {
            *x /= 2.0 * r;
            *y /= 2.0 * r;
        }

        self.local_to_global(x, y);
    }

    fn distance_normal(&self, nx: &mut f64, ny: &mut f64, mut x: f64, mut y: f64) {
        x -= self.get_fields().pos_x;
        y -= self.get_fields().pos_y;
        let r = length(x, y);
        if r < 1e-4f64 {
            *nx = 1.0;
            *ny = 0.0;
        } else {
            *nx = x / r;
            *ny = y / r;
        }
    }
}

#[derive(Clone)]
struct FluidQuantity {
    src: Vec<f64>,
    dst: Vec<f64>,
    phi: Vec<f64>,
    volume: Vec<f64>,
    normal_x: Vec<f64>,
    normal_y: Vec<f64>,
    cell: Vec<u8>,
    body: Vec<u8>,
    mask: Vec<u8>,

    w: usize,
    h: usize,
    ox: f64,
    oy: f64,
    hx: f64,
}

fn lerp(a: f64, b: f64, x: f64) -> f64 {
    a * (1.0 - x) + b * x
}

fn cerp(a: f64, b: f64, c: f64, d: f64, x: f64) -> f64 {
    let xsq = x * x;
    let xcu = xsq * x;
    let min_v = f64::min(a, f64::min(b, f64::min(c, d)));
    let max_v = f64::max(a, f64::max(b, f64::max(c, d)));

    let t = a * (0.0 - 0.5 * x + 1.0 * xsq - 0.5 * xcu)
        + b * (1.0 + 0.0 * x - 2.5 * xsq + 1.5 * xcu)
        + c * (0.0 + 0.5 * x + 2.0 * xsq - 1.5 * xcu)
        + d * (0.0 + 0.0 * x - 0.5 * xsq + 0.5 * xcu);

    f64::min(f64::max(t, min_v), max_v)
}

impl FluidQuantity {
    fn runge_kutta_3(
        &self,
        x: &mut f64,
        y: &mut f64,
        timestep: f64,
        u: &FluidQuantity,
        v: &FluidQuantity,
    ) {
        let first_u = u.lerp(*x, *y) / self.hx;
        let first_v = v.lerp(*x, *y) / self.hx;

        let mid_x = *x - 0.5 * timestep * first_u;
        let mid_y = *y - 0.5 * timestep * first_v;

        let mid_u = u.lerp(mid_x, mid_y) / self.hx;
        let mid_v = v.lerp(mid_x, mid_y) / self.hx;

        let last_x = *x - 0.75 * timestep * mid_u;
        let last_y = *y - 0.75 * timestep * mid_v;

        let last_u = u.lerp(last_x, last_y);
        let last_v = v.lerp(last_x, last_y);

        *x -= timestep * ((2.0 / 9.0) * first_u + (3.0 / 9.0) * mid_u + (4.0 / 9.0) * last_u);
        *y -= timestep * ((2.0 / 9.0) * first_v + (3.0 / 9.0) * mid_v + (4.0 / 9.0) * last_v);
    }

    fn new(w: usize, h: usize, ox: f64, oy: f64, hx: f64) -> Self {
        Self {
            src: vec![0f64; w * h],
            dst: vec![0f64; w * h],
            phi: vec![0f64; (w + 1) * (h + 1)],
            volume: vec![1f64; w * h],
            normal_x: vec![0f64; w * h],
            normal_y: vec![0f64; w * h],
            cell: vec![CELL_FLUID; w * h],
            body: vec![0u8; w * h],
            mask: vec![0u8; w * h],

            w,
            h,
            ox,
            oy,
            hx,
        }
    }

    fn flip(&mut self) {
        mem::swap(&mut self.src, &mut self.dst);
    }

    fn src(&self) -> &Vec<f64> {
        &self.src
    }

    fn src_mut(&mut self) -> &mut Vec<f64> {
        &mut self.src
    }

    fn cell(&self) -> &Vec<u8> {
        &self.cell
    }
    fn body(&self) -> &Vec<u8> {
        &self.body
    }

    fn idx(&self, x: usize, y: usize) -> usize {
        x + y * self.w
    }

    fn at(&self, x: usize, y: usize) -> f64 {
        self.src[(x + y * self.w)]
    }

    fn at_mut(&mut self, x: usize, y: usize) -> &mut f64 {
        &mut self.src[(x + y * self.w)]
    }

    fn volume(&self, x: usize, y: usize) -> f64 {
        self.volume[(x + y * self.w)]
    }

    fn lerp(&self, mut x: f64, mut y: f64) -> f64 {
        x = f64::min(f64::max(x - self.ox, 0.0), self.w as f64 - 1.001);
        y = f64::min(f64::max(y - self.oy, 0.0), self.h as f64 - 1.001);
        let ix = x as usize;
        let iy = y as usize;
        x -= ix as f64;
        y -= iy as f64;

        let x00 = self.at(ix + 0, iy + 0);
        let x10 = self.at(ix + 1, iy + 0);
        let x01 = self.at(ix + 0, iy + 1);
        let x11 = self.at(ix + 1, iy + 1);

        lerp(lerp(x00, x10, x), lerp(x01, x11, x), y)
    }

    fn cerp(&self, mut x: f64, mut y: f64) -> f64 {
        x = f64::min(f64::max(x - self.ox, 0.0), self.w as f64 - 1.001);
        y = f64::min(f64::max(y - self.oy, 0.0), self.h as f64 - 1.001);
        let ix = x as usize;
        let iy = y as usize;
        x -= ix as f64;
        y -= iy as f64;
        let x0 = i64::max(ix as i64 - 1, 0) as usize;
        let x1 = ix;
        let x2 = ix + 1;
        let x3 = i64::min(ix as i64 + 2, self.w as i64 - 1) as usize;
        let y0 = i64::max(iy as i64 - 1, 0) as usize;
        let y1 = iy;
        let y2 = iy + 1;
        let y3 = usize::min(iy + 2, self.h - 1);

        let q0 = cerp(
            self.at(x0, y0),
            self.at(x1, y0),
            self.at(x2, y0),
            self.at(x3, y0),
            x,
        );
        let q1 = cerp(
            self.at(x0, y1),
            self.at(x1, y1),
            self.at(x2, y1),
            self.at(x3, y1),
            x,
        );
        let q2 = cerp(
            self.at(x0, y2),
            self.at(x1, y2),
            self.at(x2, y2),
            self.at(x3, y2),
            x,
        );
        let q3 = cerp(
            self.at(x0, y3),
            self.at(x1, y3),
            self.at(x2, y3),
            self.at(x3, y3),
            x,
        );
        cerp(q0, q1, q2, q3, y)
    }

    fn back_project(&self, x: &mut f64, y: &mut f64, bodies: &Vec<Box<dyn SolidBody>>) {
        let rx = usize::min(usize::max((*x - self.ox) as usize, 0), self.w - 1);
        let ry = usize::min(usize::max((*y - self.oy) as usize, 0), self.h - 1);

        if self.cell[(rx + ry * self.w)] != CELL_FLUID {
            *x = (*x - self.ox) * self.hx;
            *y = (*y - self.oy) * self.hx;
            bodies[self.body[(rx + ry * self.w)] as usize].closest_surface_point(x, y);
            *x = *x / self.hx + self.ox;
            *y = *y / self.hx + self.oy;
        }
    }

    fn advect(
        &mut self,
        timestep: f64,
        u: &FluidQuantity,
        v: &FluidQuantity,
        bodies: &Vec<Box<dyn SolidBody>>,
    ) {
        for iy in 0..self.h {
            for ix in 0..self.w {
                let idx = iy * self.w + ix;
                if self.cell[idx] == CELL_FLUID {
                    let mut x = ix as f64 + self.ox;
                    let mut y = iy as f64 + self.oy;
                    self.runge_kutta_3(&mut x, &mut y, timestep, u, v);
                    self.back_project(&mut x, &mut y, bodies);
                    self.dst[idx] = self.cerp(x, y);
                }
            }
        }
    }

    fn add_inflow(&mut self, x0: f64, y0: f64, x1: f64, y1: f64, v: f64) {
        let ix0 = (x0 / self.hx - self.ox) as usize;
        let iy0 = (y0 / self.hx - self.oy) as usize;
        let ix1 = (x1 / self.hx - self.ox) as usize;
        let iy1 = (y1 / self.hx - self.oy) as usize;
        // TODO: Might bug
        for y in usize::max(iy0, 0)..usize::min(iy1, self.h) {
            for x in usize::max(ix0, 0)..usize::min(ix1, self.h) {
                let l = length(
                    (2.0 * (x as f64 + 0.5) * self.hx - (x0 + x1)) / (x1 - x0),
                    (2.0 * (y as f64 + 0.5) * self.hx - (y0 + y1)) / (y1 - y0),
                );

                let vi = cubic_pulse(l) * v;

                if f64::abs(self.src[(x + y * self.w)]) < f64::abs(vi) {
                    self.src[(x + y * self.w)] = vi;
                }
            }
        }
    }

    fn fill_solid_fields(&mut self, bodies: &Vec<Box<dyn SolidBody>>) {
        if bodies.is_empty() {
            return;
        }

        for iy in 0..self.h {
            for ix in 0..self.w {
                let idx = iy * self.w + ix;
                let x = (ix as f64 + self.ox - 0.5) * self.hx;
                let y = (iy as f64 + self.oy - 0.5) * self.hx;
                self.phi[idx] = bodies[0].distance(x, y);
                for i in 1..bodies.len() {
                    self.phi[idx] = f64::min(self.phi[idx], bodies[i].distance(x, y));
                }
            }
        }

        for iy in 0..self.h {
            for ix in 0..self.w {
                let idx = iy * self.w + ix;
                let x = (ix as f64 + self.ox) * self.hx;
                let y = (iy as f64 + self.oy) * self.hx;

                self.body[idx] = 0;

                let mut d = bodies[0].distance(x, y);

                for i in 0..bodies.len() {
                    let id = bodies[i].distance(x, y);
                    if id < d {
                        self.body[idx] = i as u8;
                        d = id;
                    }
                }

                let idxp = ix + iy * (self.w + 1);

                self.volume[idx] = 1.0
                    - occupancy(
                        self.phi[(idxp)],
                        self.phi[(idxp + 1)],
                        self.phi[(idxp + self.w + 1)],
                        self.phi[(idxp + self.w + 2)],
                    );

                if self.volume[idx] < 0.01 {
                    self.volume[idx] = 0.0;
                }

                bodies[self.body[idx] as usize].distance_normal(
                    &mut self.normal_x[idx],
                    &mut self.normal_y[idx],
                    x,
                    y,
                );

                self.cell[idx] = if self.volume[idx] == 0.0 {
                    CELL_SOLID
                } else {
                    CELL_FLUID
                };
            }
        }
    }

    fn fill_solid_mask(&mut self) {
        for y in 1..self.h - 1 {
            for x in 1..self.w - 1 {
                let idx = x + y * self.w;

                if self.cell[idx] == CELL_FLUID {
                    continue;
                }
                let nx = self.normal_x[idx];
                let ny = self.normal_y[idx];

                self.mask[idx] = 0;

                if nx != 0.0 && self.cell[(idx as f64 + f64::signum(nx)) as usize] != CELL_FLUID {
                    self.mask[idx] |= 1;
                }
                if ny != 0.0
                    && self.cell[(idx as f64 + f64::signum(ny) * self.w as f64) as usize]
                        != CELL_FLUID
                {
                    self.mask[idx] |= 2;
                }
            }
        }
    }

    fn extrapolate_normal(&self, idx: usize) -> f64 {
        let nx = self.normal_x[idx];
        let ny = self.normal_y[idx];
        let src_x = self.src[(idx as f64 + f64::signum(nx)) as usize];
        let src_y = self.src[(idx as f64 + f64::signum(ny) * self.w as f64) as usize];

        (f64::abs(nx) * src_x + f64::abs(ny) * src_y) / (f64::abs(nx) + f64::abs(ny))
    }

    fn free_neighbour(&mut self, idx: usize, border: &mut Vec<usize>, mask: u8) {
        self.mask[idx] &= !mask;

        if self.cell[idx] != CELL_FLUID && self.mask[idx] == 0 {
            border.push(idx);
        }
    }

    fn extrapolate(&mut self) {
        self.fill_solid_mask();

        let mut border = Vec::new();

        for y in 1..self.h - 1 {
            for x in 1..self.w {
                let idx = x + y * self.w;

                if self.cell[idx] != CELL_FLUID && self.mask[idx] == 0 {
                    border.push(idx);
                }
            }
        }

        while !border.is_empty() {
            let idx = border.pop().unwrap();

            self.src[idx] = self.extrapolate_normal(idx);

            if self.normal_x[(idx - 1)] > 0.0 {
                self.free_neighbour(idx - 1, &mut border, 1);
            }
            if self.normal_x[(idx + 1)] < 0.0 {
                self.free_neighbour(idx + 1, &mut border, 1);
            }
            if self.normal_y[(idx - self.w)] > 0.0 {
                self.free_neighbour(idx - self.w, &mut border, 2);
            }
            if self.normal_y[(idx + self.w)] < 0.0 {
                self.free_neighbour(idx + self.w, &mut border, 2);
            }
        }
    }
}

struct FluidSolver {
    d: FluidQuantity,
    t: FluidQuantity,
    u: FluidQuantity,
    v: FluidQuantity,
    /* Densities at staggered grid locations */
    u_density: Vec<f64>,
    v_density: Vec<f64>,
    w: usize,
    h: usize,
    hx: f64,
    density_air: f64,
    density_soot: f64,
    diffusion: f64,
    r: Vec<f64>,
    p: Vec<f64>,
    z: Vec<f64>,
    s: Vec<f64>,
    precon: Vec<f64>,
    a_diag: Vec<f64>,
    a_plus_x: Vec<f64>,
    a_plus_y: Vec<f64>,
    t_amb: f64,
    g: f64,
    bodies: Vec<Box<dyn SolidBody>>,
}

impl FluidSolver {
    fn build_rhs(&mut self) {
        let scale = 1.0 / self.hx;

        let cell = self.d.cell();
        let body = self.d.body();

        for y in 0..self.h {
            for x in 0..self.w {
                let idx = y * self.w + x;

                if cell[idx] == CELL_FLUID {
                    self.r[idx] = -scale
                        * (self.u.volume(x + 1, y) * self.u.at(x + 1, y)
                            - self.u.volume(x, y) * self.u.at(x, y)
                            + self.v.volume(x, y + 1) * self.v.at(x, y + 1)
                            - self.v.volume(x, y) * self.v.at(x, y));
                    let vol = self.d.volume(x, y);
                    if self.bodies.is_empty() {
                        continue;
                    }

                    if x > 0 {
                        self.r[idx] -= (self.u.volume(x, y) - vol)
                            * self.bodies[body[idx - 1] as usize]
                                .velocity_x(x as f64 * self.hx, (y as f64 + 0.5) * self.hx);
                    }
                    if y > 0 {
                        self.r[idx] -= (self.v.volume(x, y) - vol)
                            * self.bodies[body[idx - self.w] as usize]
                                .velocity_y((x as f64 + 0.5) * self.hx, y as f64 * self.hx);
                    }
                    if x < self.w - 1 {
                        self.r[idx] += (self.u.volume(x + 1, y) - vol)
                            * self.bodies[body[idx + 1] as usize]
                                .velocity_x((x as f64 + 1.0) * self.hx, (y as f64 + 0.5) * self.hx);
                    }
                    if y < self.h - 1 {
                        self.r[idx] += (self.v.volume(x, y + 1) - vol)
                            * self.bodies[body[idx + self.w] as usize]
                                .velocity_y((x as f64 + 0.5) * self.hx, (y as f64 + 1.0) * self.hx);
                    }
                } else {
                    self.r[idx] = 0.0;
                }
            }
        }
    }

    /* Computes densities at the staggered grid locations as a function of
     * temperature and smoke concentration.
     */
    fn compute_densities(&mut self) {
        let alpha = (self.density_soot - self.density_air) / self.density_air;

        unsafe {
            write_bytes(self.u_density.as_mut_ptr(), 0, (self.w + 1) * self.h);
            write_bytes(self.v_density.as_mut_ptr(), 0, self.w * (self.h + 1));
        }

        for y in 0..self.h {
            for x in 0..self.w {
                let mut density = self.density_air * self.t_amb / self.t.at(x, y)
                    * (1.0 + alpha * self.d.at(x, y));
                density = f64::max(density, 0.05 * self.density_air); /* Clamp dangerously low densities */
                self.u_density[self.u.idx(x, y)] += 0.5 * density;
                self.v_density[self.v.idx(x, y)] += 0.5 * density;
                self.u_density[self.u.idx(x + 1, y)] += 0.5 * density;
                self.v_density[self.v.idx(x, y + 1)] += 0.5 * density;
            }
        }
    }

    /* Instead of constant density per cell, the entries must now be modified
     * to account for variable density at individual grid cells.
     */
    fn build_pressure_matrix(&mut self, timestep: f64) {
        let scale = timestep / (self.hx * self.hx);
        let cell = self.d.cell();

        unsafe {
            write_bytes(self.a_diag.as_mut_ptr(), 0, self.w * self.h);
            write_bytes(self.a_plus_x.as_mut_ptr(), 0, self.w * self.h);
            write_bytes(self.a_plus_y.as_mut_ptr(), 0, self.w * self.h);
        }

        for y in 0..self.h {
            for x in 0..self.w {
                let idx = y * self.w + x;

                if cell[idx] != CELL_FLUID {
                    continue;
                }

                if x < self.w - 1 && cell[idx + 1] == CELL_FLUID {
                    let factor =
                        scale * self.u.volume(x + 1, y) / self.u_density[self.u.idx(x + 1, y)];
                    self.a_diag[idx] += factor;
                    self.a_diag[idx + 1] += factor;
                    self.a_plus_x[idx] = -factor;
                }
                if y < self.h - 1 && cell[idx + self.w] == CELL_FLUID {
                    let factor =
                        scale * self.v.volume(x, y + 1) / self.v_density[self.u.idx(x, y + 1)];
                    self.a_diag[idx] += factor;
                    self.a_diag[idx + self.w] += factor;
                    self.a_plus_y[idx] = -factor;
                }
            }
        }
    }

    fn build_heat_diffusion_matrix(&mut self, timestep: f64) {
        for i in 0..self.w * self.h {
            self.a_diag[i] = 1.0;
        }

        unsafe {
            write_bytes(self.a_plus_x.as_mut_ptr(), 0, self.w * self.h);
            write_bytes(self.a_plus_y.as_mut_ptr(), 0, self.w * self.h);
        }

        let cell = self.d.cell();
        let scale = self.diffusion * timestep * 1.0 / (self.hx * self.hx);

        for y in 0..self.h {
            for x in 0..self.w {
                let idx = y * self.w + x;

                if cell[idx] != CELL_FLUID {
                    continue;
                }
                if x < self.w - 1 && cell[idx + 1] == CELL_FLUID {
                    self.a_diag[idx] += scale;
                    self.a_diag[idx + 1] += scale;
                    self.a_plus_x[idx] = -scale;
                }
                if y < self.h - 1 && cell[idx + self.w] == CELL_FLUID {
                    self.a_diag[idx] += scale;
                    self.a_diag[idx + self.w] += scale;
                    self.a_plus_y[idx] = -scale;
                }
            }
        }
    }

    fn build_preconditioner(&mut self) {
        const TAU: f64 = 0.97;
        const SIGMA: f64 = 0.25;

        let cell = self.d.cell();
        for y in 0..self.h {
            for x in 0..self.w {
                let idx = y * self.w + x;

                if cell[idx] != CELL_FLUID {
                    continue;
                }

                let mut e = self.a_diag[idx];

                if x > 0 && cell[idx - 1] == CELL_FLUID {
                    let px = self.a_plus_x[idx - 1] * self.precon[idx - 1];
                    let py = self.a_plus_y[idx - 1] * self.precon[idx - 1];
                    e = e - (px * px + TAU * px * py);
                }
                if y > 0 && cell[idx - self.w] == CELL_FLUID {
                    let px = self.a_plus_x[idx - self.w] * self.precon[idx - self.w];
                    let py = self.a_plus_y[idx - self.w] * self.precon[idx - self.w];
                    e = e - (py * py + TAU * px * py);
                }

                if e < SIGMA * self.a_diag[idx] {
                    e = self.a_diag[idx];
                }

                self.precon[idx] = 1.0 / f64::sqrt(e);
            }
        }
    }

    fn apply_preconditioner(&self, dst: &mut Vec<f64>, a: &Vec<f64>) {
        let cell = self.d.cell();

        for y in 0..self.h {
            for x in 0..self.w {
                let idx = y * self.w + x;

                if cell[idx] != CELL_FLUID {
                    continue;
                }

                let mut t = a[idx];

                if x > 0 && cell[idx - 1] == CELL_FLUID {
                    t -= self.a_plus_x[idx - 1] * self.precon[idx - 1] * dst[idx - 1];
                }
                if y > 0 && cell[idx - self.w] == CELL_FLUID {
                    t -=
                        self.a_plus_y[idx - self.w] * self.precon[idx - self.w] * dst[idx - self.w];
                }

                dst[idx] = t * self.precon[idx];
            }
        }

        // TODO: Might bug
        for y in self.h - 1..=0 {
            for x in self.w - 1..=0 {
                let idx = y * self.w + x;

                if cell[idx] != CELL_FLUID {
                    continue;
                }

                let mut t = dst[idx];

                if x < self.w - 1 && cell[idx + 1] == CELL_FLUID {
                    t -= self.a_plus_x[idx] * self.precon[idx] * dst[idx + 1];
                }
                if y < self.h - 1 && cell[idx + self.w] == CELL_FLUID {
                    t -= self.a_plus_y[idx] * self.precon[idx] * dst[idx + self.w];
                }

                dst[idx] = t * self.precon[idx];
            }
        }
    }

    fn dot_product(&self, a: &Vec<f64>, b: &Vec<f64>) -> f64 {
        let cell = self.d.cell();
        let mut result = 0.0;

        for i in 0..(self.w * self.h) {
            if cell[i] == CELL_FLUID {
                result += a[i] * b[i];
            }
        }

        result
    }

    fn matrix_vector_product(&self, dst: &mut Vec<f64>, b: &Vec<f64>) {
        for y in 0..self.h {
            for x in 0..self.w {
                let idx = y * self.w + x;
                let mut t = self.a_diag[idx] * b[idx];
                if x > 0 {
                    t += self.a_plus_x[idx - 1] * b[idx - 1];
                }
                if y > 0 {
                    t += self.a_plus_y[idx - self.w] * b[idx - self.w];
                }
                if x < self.w - 1 {
                    t += self.a_plus_x[idx] * b[idx + 1];
                }
                if y < self.h - 1 {
                    t += self.a_plus_y[idx] * b[idx + self.w];
                }

                dst[idx] = t;
            }
        }
    }

    fn scaled_add(&self, dst: &mut Vec<f64>, a: &Vec<f64>, b: &Vec<f64>, s: f64) {
        let cell = self.d.cell();
        for i in 0..(self.w * self.h) {
            if cell[i] == CELL_FLUID {
                dst[i] = a[i] + b[i] * s;
            }
        }
    }

    fn scaled_add_in_place_a(&self, dst: &mut Vec<f64>, b: &Vec<f64>, s: f64) {
        let cell = self.d.cell();
        for i in 0..(self.w * self.h) {
            if cell[i] == CELL_FLUID {
                dst[i] = dst[i] + b[i] * s;
            }
        }
    }

    fn scaled_add_in_place_b(&self, dst: &mut Vec<f64>, a: &Vec<f64>, s: f64) {
        let cell = self.d.cell();
        for i in 0..(self.w * self.h) {
            if cell[i] == CELL_FLUID {
                dst[i] = a[i] + dst[i] * s;
            }
        }
    }

    fn infinity_norm(&self, a: &Vec<f64>) -> f64 {
        let cell = self.d.cell();
        let mut max_a = 0.0;

        for i in 0..(self.w * self.h) {
            if cell[i] == CELL_FLUID {
                max_a = f64::max(max_a, f64::abs(a[i]));
            }
        }

        max_a
    }

    fn project(&mut self, limit: usize) {
        unsafe {
            write_bytes(self.p.as_mut_ptr(), 0, self.w * self.h);
        }

        let mut z = &mut self.z.clone();
        let mut p = &mut self.p.clone();
        let mut r = &mut self.r.clone();
        let mut s = &mut self.s.clone();

        self.apply_preconditioner(&mut z, r);

        unsafe {
            copy_nonoverlapping(z.as_mut_ptr(), s.as_mut_ptr(), self.w * self.h);
        }

        let mut max_error = self.infinity_norm(r);

        if max_error < 1e-5 {
            println!("Initial guess sufficiently small");
            self.z = z.clone();
            self.p = p.clone();
            self.r = r.clone();
            self.s = s.clone();
            return;
        }

        let mut sigma = self.dot_product(z, r);

        for iter in 0..limit {
            self.matrix_vector_product(&mut z, s);
            let alpha = sigma / self.dot_product(z, s);

            self.scaled_add_in_place_a(&mut p, s, alpha);
            self.scaled_add_in_place_a(&mut r, z, -alpha);

            max_error = self.infinity_norm(r);

            if max_error < 1e-5 {
                println!(
                    "Exiting solver after {} iterations, maximum error is {}",
                    iter, max_error
                );
                self.z = z.clone();
                self.p = p.clone();
                self.r = r.clone();
                self.s = s.clone();
                return;
            }

            self.apply_preconditioner(&mut z, r);
            let sigma_new = self.dot_product(z, r);
            self.scaled_add_in_place_b(&mut s, z, sigma_new / sigma);
            sigma = sigma_new;
        }

        self.z = z.clone();
        self.p = p.clone();
        self.r = r.clone();
        self.s = s.clone();
        println!(
            "Exceeded budget of {} iterations, maximum error was {}",
            limit, max_error
        );
    }

    /* Similar to the pressure matrix, we cannot assume constant density per
     * cell here either and must modify the equations accordingly.
     */
    fn apply_pressure(&mut self, timestep: f64) {
        let scale = timestep / self.hx;
        let cell = self.d.cell();

        for y in 0..self.h {
            for x in 0..self.w {
                let idx = y * self.w + x;
                if cell[idx] != CELL_FLUID {
                    continue;
                }
                *self.u.at_mut(x, y) -= scale * self.p[idx] / self.u_density[self.u.idx(x, y)];
                *self.v.at_mut(x, y) -= scale * self.p[idx] / self.v_density[self.v.idx(x, y)];
                *self.u.at_mut(x + 1, y) +=
                    scale * self.p[idx] / self.u_density[self.u.idx(x + 1, y)];
                *self.v.at_mut(x, y + 1) +=
                    scale * self.p[idx] / self.v_density[self.v.idx(x, y + 1)];
            }
        }
    }

    fn add_buoyancy(&mut self, timestep: f64) {
        let alpha = (self.density_soot - self.density_air) / self.density_air;

        for y in 0..self.h {
            for x in 0..self.w {
                let buoyancy = timestep
                    * self.g
                    * (alpha * self.d.at(x, y) - (self.t.at(x, y) - self.t_amb) / self.t_amb);

                *self.v.at_mut(x, y) += buoyancy * 0.5;
                *self.v.at_mut(x, y + 1) += buoyancy * 0.5;
            }
        }
    }

    fn set_boundary_condition(&mut self) {
        let cell = self.d.cell();
        let body = self.d.body();

        for y in 0..self.h {
            for x in 0..self.w {
                let idx = y * self.w + x;
                if cell[idx] == CELL_SOLID {
                    let b = &self.bodies[body[idx] as usize];
                    *self.u.at_mut(x, y) =
                        b.velocity_x(x as f64 * self.hx, (y as f64 + 0.5) * self.hx);
                    *self.v.at_mut(x, y) =
                        b.velocity_y((x as f64 + 0.5) * self.hx, y as f64 * self.hx);
                    *self.u.at_mut(x + 1, y) =
                        b.velocity_x((x as f64 + 1.0) * self.hx, (y as f64 + 0.5) * self.hx);
                    *self.v.at_mut(x, y + 1) =
                        b.velocity_y((x as f64 + 0.5) * self.hx, (y as f64 + 1.0) * self.hx);
                }
            }
        }
        for y in 0..self.h {
            *self.u.at_mut(0, y) = 0.0;
            *self.u.at_mut(self.w, y) = 0.0;
        }
        for x in 0..self.w {
            *self.v.at_mut(x, 0) = 0.0;
            *self.v.at_mut(x, self.h) = 0.0;
        }
    }

    fn new(
        w: usize,
        h: usize,
        rho_air: f64,
        rho_soot: f64,
        diffusion: f64,
        bodies: Vec<Box<dyn SolidBody>>,
    ) -> Self {
        let hx = 1.0 / usize::min(w, h) as f64;
        let mut ret = Self {
            w,
            h,
            density_air: rho_air,
            density_soot: rho_soot,
            diffusion: diffusion,
            bodies: bodies,
            t_amb: 294.0,
            g: 9.81,
            hx,
            d: FluidQuantity::new(w, h, 0.5, 0.5, hx),
            t: FluidQuantity::new(w, h, 0.5, 0.5, hx),
            u: FluidQuantity::new(w + 1, h, 0.0, 0.5, hx),
            v: FluidQuantity::new(w, h + 1, 0.5, 0.0, hx),
            r: vec![0f64; w * h],
            p: vec![0f64; w * h],
            z: vec![0f64; w * h],
            s: vec![0f64; w * h],
            a_diag: vec![0f64; w * h],
            a_plus_x: vec![0f64; w * h],
            a_plus_y: vec![0f64; w * h],
            precon: vec![0f64; w * h],
            u_density: vec![0f64; (w + 1) * h],
            v_density: vec![0f64; w * (h + 1)],
        };

        for i in 0..w * h {
            ret.t.src_mut()[i] = ret.t_amb;
        }

        ret
    }

    fn update(&mut self, timestep: f64) {
        self.d.fill_solid_fields(&self.bodies);
        self.t.fill_solid_fields(&self.bodies);
        self.u.fill_solid_fields(&self.bodies);
        self.v.fill_solid_fields(&self.bodies);
        unsafe {
            copy_nonoverlapping(
                self.t.src_mut().as_mut_ptr(),
                self.r.as_mut_ptr(),
                self.w * self.h,
            );
        }
        self.build_heat_diffusion_matrix(timestep);
        self.build_preconditioner();
        self.project(2000);
        unsafe {
            copy_nonoverlapping(
                self.p.as_mut_ptr(),
                self.t.src_mut().as_mut_ptr(),
                self.w * self.h,
            );
        }
        self.t.extrapolate();
        self.add_buoyancy(timestep);
        self.set_boundary_condition();
        self.build_rhs();
        self.compute_densities();
        self.build_pressure_matrix(timestep);
        self.build_preconditioner();
        self.project(2000);
        self.apply_pressure(timestep);
        self.d.extrapolate();
        self.u.extrapolate();
        self.v.extrapolate();
        self.set_boundary_condition();
        self.d.advect(timestep, &self.u, &self.v, &self.bodies);
        self.t.advect(timestep, &self.u, &self.v, &self.bodies);
        self.u
            .advect(timestep, &self.u.clone(), &self.v, &self.bodies);
        self.v
            .advect(timestep, &self.u, &self.v.clone(), &self.bodies);
        self.d.flip();
        self.t.flip();
        self.u.flip();
        self.v.flip();
    }

    fn add_inflow(&mut self, x: f64, y: f64, w: f64, h: f64, d: f64, t: f64, u: f64, v: f64) {
        self.d.add_inflow(x, y, x + w, y + h, d);
        self.t.add_inflow(x, y, x + w, y + h, t);
        self.u.add_inflow(x, y, x + w, y + h, u);
        self.v.add_inflow(x, y, x + w, y + h, v);
    }

    fn ambient_t(&self) -> f64 {
        self.t_amb
    }

    fn to_image(&self, render_heat: bool) -> Vec<u8> {
        let mut rgba = vec![0u8; self.h * if render_heat { 2 } else { 1 } * self.w * 4];

        for y in 0..self.h {
            for x in 0..self.w {
                let (idxl, idxr) = if render_heat {
                    (4 * (x + y * self.w * 2), 4 * (x + y * self.w * 2 + self.w))
                } else {
                    (0, 4 * (x + y * self.w))
                };

                let volume = self.d.volume(x, y);
                let mut shade = (1.0 - self.d.at(x, y)) * volume;

                shade = f64::min(f64::max(shade, 0.0), 1.0);

                rgba[idxr + 0] = (shade * 255.0) as u8;
                rgba[idxr + 1] = (shade * 255.0) as u8;
                rgba[idxr + 2] = (shade * 255.0) as u8;
                rgba[idxr + 3] = 0xFF;

                if render_heat {
                    let mut t = (self.t.at(x, y) - self.t_amb) / 700.0;
                    t = f64::min(f64::max(t, 0.0), 1.0);
                    let r = 1.0 + volume * (f64::min(t * 4.0, 1.0) - 1.0);
                    let g = 1.0 + volume * (f64::min(t * 2.0, 1.0) - 1.0);
                    let b = 1.0 + volume * (f64::max(f64::min(t * 4.0 - 3.0, 1.0), 0.0) - 1.0);
                    rgba[idxl + 0] = (r * 255.0) as u8;
                    rgba[idxl + 1] = (g * 255.0) as u8;
                    rgba[idxl + 2] = (b * 255.0) as u8;
                    rgba[idxl + 3] = 0xFF;
                }
            }
        }

        rgba
    }
}

fn main() {
    /* Play with these constants, if you want */
    let size_x = 64;
    let size_y = 64;
    let density_air = 0.1;
    let density_soot = 1.0; /* You can make this smaller to get lighter smoke */
    let diffusion = 0.03;
    let timestep = 0.005;
    let render_heat = true; /* Set this to true to enable heat rendering */

    let mut bodies: Vec<Box<dyn SolidBody>> = Vec::new();

    bodies.push(Box::new(SolidBox::new(SolidBodyFields::new(
        0.5,
        0.6,
        0.9,
        0.1,
        f64::consts::PI * 0.25,
        0.0,
        0.0,
        0.1,
    ))));

    let mut solver = FluidSolver::new(size_x, size_y, density_air, density_soot, diffusion, bodies);
    let mut time = 0.0;
    let mut iterations = 0;

    while time < 8.0 {
        for _i in 0..4 {
            solver.add_inflow(
                0.45,
                0.2,
                0.1,
                0.05,
                1.0,
                650.0 + solver.ambient_t(),
                0.0,
                0.0,
            );
            solver.update(timestep);
            time += timestep;
        }
        let image = solver.to_image(render_heat);
        let path = format!("frame{:0>5}.png", iterations);

        iterations += 1;

        screenshot(
            &path,
            if render_heat { size_x * 2 } else { size_x } as u32,
            size_y as u32,
            image,
        );
        // lodepng_encode32_file(path, image, (render_heat ? size_x*2 : size_x), size_y);

        for i in 0..solver.bodies.len() {
            solver.bodies[i].update(timestep);
        }
    }
}

pub fn screenshot(path: &str, width: u32, height: u32, data: Vec<u8>) {
    println!("Writing file {}", path);

    let file = File::create(Path::new(path)).unwrap();
    let ref mut w = BufWriter::new(file);
    let mut encoder = png::Encoder::new(w, width, height);

    encoder.set_color(png::ColorType::RGBA);
    encoder.set_depth(png::BitDepth::Eight);
    encoder.set_compression(png::Compression::Default);
    encoder.set_filter(png::FilterType::NoFilter);

    let mut writer = encoder.write_header().unwrap();

    writer.write_image_data(&data[..]).unwrap();
}
