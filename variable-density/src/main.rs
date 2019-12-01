extern crate png;

use std::f64;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::mem;
use std::path::Path;
use std::ptr::copy_nonoverlapping;
use std::ptr::write_bytes;

fn length(x: f64, y: f64) -> f64 {
    f64::sqrt(x * x + y * y)
}

fn cubicPulse(mut x: f64) -> f64 {
    x = f64::min(f64::abs(x), 1.0);
    1.0 - x * x * (3.0 - 2.0 * x)
}

fn rotate(x: &mut f64, y: &mut f64, phi: f64) {
    let tmpX = *x;
    let tmpY = *y;
    *x = f64::cos(phi) * tmpX + f64::sin(phi) * tmpY;
    *y = -f64::sin(phi) * tmpY + f64::cos(phi) * tmpY;
}

fn triangleOccupancy(out1: f64, in1: f64, out2: f64) -> f64 {
    0.5 * in1 * in1 / ((out1 - in1) * (out2 - in1))
}

fn trapezoidOccupancy(out1: f64, out2: f64, in1: f64, in2: f64) -> f64 {
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

        0x1 => triangleOccupancy(d21, d11, d12),
        0x2 => triangleOccupancy(d11, d12, d22),
        0x4 => triangleOccupancy(d12, d22, d21),
        0x8 => triangleOccupancy(d22, d21, d11),

        0xE => 1.0 - triangleOccupancy(-d21, -d11, -d12),
        0xD => 1.0 - triangleOccupancy(-d11, -d12, -d22),
        0xB => 1.0 - triangleOccupancy(-d12, -d22, -d21),
        0x7 => 1.0 - triangleOccupancy(-d22, -d21, -d11),
        0x3 => trapezoidOccupancy(d21, d22, d11, d12),
        0x6 => trapezoidOccupancy(d11, d21, d12, d22),
        0x9 => trapezoidOccupancy(d12, d22, d11, d21),
        0xC => trapezoidOccupancy(d11, d12, d21, d22),

        0x5 => triangleOccupancy(d11, d12, d22) + triangleOccupancy(d22, d21, d11),
        0xA => triangleOccupancy(d21, d11, d12) + triangleOccupancy(d12, d22, d21),

        0xF => 1.0,
        _ => 0.0,
    }
}

const CELL_FLUID: u8 = 0;
const CELL_SOLID: u8 = 1;

struct SolidBodyFields {
    pub _posX: f64,
    pub _posY: f64,
    pub _scaleX: f64,
    pub _scaleY: f64,
    pub _theta: f64,
    pub _velX: f64,
    pub _velY: f64,
    pub _velTheta: f64,
}

impl SolidBodyFields {
    pub fn new(
        posX: f64,
        posY: f64,
        scaleX: f64,
        scaleY: f64,
        theta: f64,
        velX: f64,
        velY: f64,
        velTheta: f64,
    ) -> Self {
        Self {
            _posX: posX,
            _posY: posY,
            _scaleX: scaleX,
            _scaleY: scaleY,
            _theta: theta,
            _velX: velX,
            _velY: velY,
            _velTheta: velTheta,
        }
    }
}

trait SolidBody {
    fn getFields(&self) -> &SolidBodyFields;
    fn getFieldsMut(&mut self) -> &mut SolidBodyFields;

    // Private (or kind of)
    fn globalToLocal(&self, x: &mut f64, y: &mut f64) {
        *x -= self.getFields()._posX;
        *y -= self.getFields()._posY;
        rotate(x, y, -self.getFields()._theta);
        *x /= self.getFields()._scaleX;
        *y /= self.getFields()._scaleY;
    }

    fn localToGlobal(&self, x: &mut f64, y: &mut f64) {
        *x *= self.getFields()._scaleX;
        *y *= self.getFields()._scaleY;
        rotate(x, y, self.getFields()._theta);
        *x += self.getFields()._posX;
        *y += self.getFields()._posY;
    }

    // Public
    fn distance(&self, x: f64, y: f64) -> f64;
    fn closestSurfacePoint(&self, x: &mut f64, y: &mut f64);
    fn distanceNormal(&self, nx: &mut f64, ny: &mut f64, x: f64, y: f64);

    fn velocityX(&self, x: f64, y: f64) -> f64 {
        (self.getFields()._posY - y) * self.getFields()._velTheta + self.getFields()._velX
    }

    fn velocityY(&self, x: f64, y: f64) -> f64 {
        (x - self.getFields()._posX) * self.getFields()._velTheta + self.getFields()._velY
    }

    fn velocity(&self, vx: &mut f64, vy: &mut f64, x: f64, y: f64) {
        *vx = self.velocityX(x, y);
        *vy = self.velocityY(x, y);
    }

    fn update(&mut self, timestep: f64) {
        self.getFieldsMut()._posX += self.getFields()._velX * timestep;
        self.getFieldsMut()._posY += self.getFields()._velY * timestep;
        self.getFieldsMut()._theta += self.getFields()._velTheta * timestep;
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
    fn getFields(&self) -> &SolidBodyFields {
        &self.fields
    }
    fn getFieldsMut(&mut self) -> &mut SolidBodyFields {
        &mut self.fields
    }

    fn distance(&self, mut x: f64, mut y: f64) -> f64 {
        x -= self.getFields()._posX;
        y -= self.getFields()._posY;
        rotate(&mut x, &mut y, -self.getFields()._theta);
        let dx = f64::abs(x) - self.getFields()._scaleX * 0.5;
        let dy = f64::abs(y) - self.getFields()._scaleY * 0.5;

        if dx >= 0.0 || dy >= 0.0 {
            length(f64::max(dx, 0.0), f64::max(dy, 0.0))
        } else {
            f64::max(dx, dy)
        }
    }

    fn closestSurfacePoint(&self, x: &mut f64, y: &mut f64) {
        *x -= self.getFields()._posX;
        *y -= self.getFields()._posY;
        rotate(x, y, -self.getFields()._theta);
        let dx = f64::abs(*x) - self.getFields()._scaleX * 0.5;
        let dy = f64::abs(*y) - self.getFields()._scaleY * 0.5;

        if dx > dy {
            *x = f64::signum(*x) * 0.5 * self.getFields()._scaleX;
        } else {
            *y = f64::signum(*y) * 0.5 * self.getFields()._scaleY;
        }

        rotate(x, y, self.getFields()._theta);
        *x += self.getFields()._posX;
        *y += self.getFields()._posY;
    }

    fn distanceNormal(&self, nx: &mut f64, ny: &mut f64, mut x: f64, mut y: f64) {
        x -= self.getFields()._posX;
        y -= self.getFields()._posY;
        rotate(&mut x, &mut y, -self.getFields()._theta);

        if f64::abs(x) - self.getFields()._scaleX * 0.5
            > f64::abs(y) - self.getFields()._scaleY * 0.5
        {
            *nx = f64::signum(x);
            *ny = 0.0;
        } else {
            *nx = 0.0;
            *ny = f64::signum(y);
        }

        rotate(nx, ny, self.getFields()._theta);
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
    fn getFields(&self) -> &SolidBodyFields {
        &self.fields
    }
    fn getFieldsMut(&mut self) -> &mut SolidBodyFields {
        &mut self.fields
    }

    fn distance(&self, x: f64, y: f64) -> f64 {
        length(x - self.getFields()._posX, y - self.getFields()._posY)
            - self.getFields()._scaleX * 0.5
    }

    fn closestSurfacePoint(&self, x: &mut f64, y: &mut f64) {
        self.globalToLocal(x, y);
        let r = length(*x, *y);
        if r < 1e-4f64 {
            *x = 0.5;
            *y = 0.0;
        } else {
            *x /= 2.0 * r;
            *y /= 2.0 * r;
        }

        self.localToGlobal(x, y);
    }

    fn distanceNormal(&self, nx: &mut f64, ny: &mut f64, mut x: f64, mut y: f64) {
        x -= self.getFields()._posX;
        y -= self.getFields()._posY;
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
    _src: Vec<f64>,
    _dst: Vec<f64>,
    _phi: Vec<f64>,
    _volume: Vec<f64>,
    _normalX: Vec<f64>,
    _normalY: Vec<f64>,
    _cell: Vec<u8>,
    _body: Vec<u8>,
    _mask: Vec<u8>,

    _w: u32,
    _h: u32,
    _ox: f64,
    _oy: f64,
    _hx: f64,
}

fn lerp(a: f64, b: f64, x: f64) -> f64 {
    a * (1.0 - x) + b * x
}

fn cerp(a: f64, b: f64, c: f64, d: f64, x: f64) -> f64 {
    let xsq = x * x;
    let xcu = xsq * x;
    let minV = f64::min(a, f64::min(b, f64::min(c, d)));
    let maxV = f64::max(a, f64::max(b, f64::max(c, d)));

    let t = a * (0.0 - 0.5 * x + 1.0 * xsq - 0.5 * xcu)
        + b * (1.0 + 0.0 * x - 2.5 * xsq + 1.5 * xcu)
        + c * (0.0 + 0.5 * x + 2.0 * xsq - 1.5 * xcu)
        + d * (0.0 + 0.0 * x - 0.5 * xsq + 0.5 * xcu);

    f64::min(f64::max(t, minV), maxV)
}

impl FluidQuantity {
    fn rungeKutta3(
        &self,
        x: &mut f64,
        y: &mut f64,
        timestep: f64,
        u: &FluidQuantity,
        v: &FluidQuantity,
    ) {
        let firstU = u.lerp(*x, *y) / self._hx;
        let firstV = v.lerp(*x, *y) / self._hx;

        let midX = *x - 0.5 * timestep * firstU;
        let midY = *y - 0.5 * timestep * firstV;

        let midU = u.lerp(midX, midY) / self._hx;
        let midV = v.lerp(midX, midY) / self._hx;

        let lastX = *x - 0.75 * timestep * midU;
        let lastY = *y - 0.75 * timestep * midV;

        let lastU = u.lerp(lastX, lastY);
        let lastV = v.lerp(lastX, lastY);

        *x -= timestep * ((2.0 / 9.0) * firstU + (3.0 / 9.0) * midU + (4.0 / 9.0) * lastU);
        *y -= timestep * ((2.0 / 9.0) * firstV + (3.0 / 9.0) * midV + (4.0 / 9.0) * lastV);
    }

    fn new(w: u32, h: u32, ox: f64, oy: f64, hx: f64) -> Self {
        Self {
            _src: vec![0f64; (w * h) as usize],
            _dst: vec![0f64; (w * h) as usize],
            _phi: vec![0f64; ((w + 1) * (h + 1)) as usize],
            _volume: vec![1f64; (w * h) as usize],
            _normalX: vec![0f64; (w * h) as usize],
            _normalY: vec![0f64; (w * h) as usize],
            _cell: vec![CELL_FLUID; (w * h) as usize],
            _body: vec![0u8; (w * h) as usize],
            _mask: vec![0u8; (w * h) as usize],

            _w: w,
            _h: h,
            _ox: ox,
            _oy: oy,
            _hx: hx,
        }
    }

    fn flip(&mut self) {
        mem::swap(&mut self._src, &mut self._dst);
    }

    fn src(&self) -> &Vec<f64> {
        &self._src
    }

    fn src_mut(&mut self) -> &mut Vec<f64> {
        &mut self._src
    }

    fn cell(&self) -> &Vec<u8> {
        &self._cell
    }
    fn body(&self) -> &Vec<u8> {
        &self._body
    }
    fn idx(&self, x: u32, y: u32) -> u32 {
        x + y * self._w
    }

    fn at(&self, x: u32, y: u32) -> f64 {
        self._src[(x + y * self._w) as usize]
    }

    fn at_mut(&mut self, x: u32, y: u32) -> &mut f64 {
        &mut self._src[(x + y * self._w) as usize]
    }

    fn volume(&self, x: u32, y: u32) -> f64 {
        self._volume[(x + y * self._w) as usize]
    }

    fn lerp(&self, mut x: f64, mut y: f64) -> f64 {
        x = f64::min(f64::max(x - self._ox, 0.0), self._w as f64 - 1.001);
        y = f64::min(f64::max(y - self._oy, 0.0), self._h as f64 - 1.001);
        let ix = x as u32;
        let iy = y as u32;
        x -= ix as f64;
        y -= iy as f64;

        let x00 = self.at(ix + 0, iy + 0);
        let x10 = self.at(ix + 1, iy + 0);
        let x01 = self.at(ix + 0, iy + 1);
        let x11 = self.at(ix + 1, iy + 1);

        lerp(lerp(x00, x10, x), lerp(x01, x11, x), y)
    }

    fn cerp(&self, mut x: f64, mut y: f64) -> f64 {
        x = f64::min(f64::max(x - self._ox, 0.0), self._w as f64 - 1.001);
        y = f64::min(f64::max(y - self._oy, 0.0), self._h as f64 - 1.001);
        let ix = x as u32;
        let iy = y as u32;
        x -= ix as f64;
        y -= iy as f64;
        let x0 = i64::max(ix as i64 - 1, 0) as u32;
        let x1 = ix;
        let x2 = ix + 1;
        let x3 = i64::min(ix as i64 + 2, self._w as i64 - 1) as u32;
        let y0 = i64::max(iy as i64 - 1, 0) as u32;
        let y1 = iy;
        let y2 = iy + 1;
        let y3 = u32::min(iy + 2, self._h - 1);

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

    fn backProject(&self, x: &mut f64, y: &mut f64, bodies: &Vec<Box<dyn SolidBody>>) {
        let rx = u32::min(u32::max((*x - self._ox) as u32, 0), self._w - 1);
        let ry = u32::min(u32::max((*y - self._oy) as u32, 0), self._h - 1);

        if self._cell[(rx + ry * self._w) as usize] != CELL_FLUID {
            *x = (*x - self._ox) * self._hx;
            *y = (*y - self._oy) * self._hx;
            bodies[self._body[(rx + ry * self._w) as usize] as usize].closestSurfacePoint(x, y);
            *x = *x / self._hx + self._ox;
            *y = *y / self._hx + self._oy;
        }
    }

    fn advect(
        &mut self,
        timestep: f64,
        u: &FluidQuantity,
        v: &FluidQuantity,
        bodies: &Vec<Box<SolidBody>>,
    ) {
        for iy in 0..self._h {
            for ix in 0..self._w {
                let idx = (iy * self._w + ix) as usize;
                if (self._cell[idx] == CELL_FLUID) {
                    let mut x = ix as f64 + self._ox;
                    let mut y = iy as f64 + self._oy;
                    self.rungeKutta3(&mut x, &mut y, timestep, u, v);
                    self.backProject(&mut x, &mut y, bodies);
                    self._dst[idx] = self.cerp(x, y);
                }
            }
        }
    }

    fn addInflow(&mut self, x0: f64, y0: f64, x1: f64, y1: f64, v: f64) {
        let ix0 = (x0 / self._hx - self._ox) as u32;
        let iy0 = (y0 / self._hx - self._oy) as u32;
        let ix1 = (x1 / self._hx - self._ox) as u32;
        let iy1 = (y1 / self._hx - self._oy) as u32;
        // TODO: Might bug
        for y in u32::max(iy0, 0)..u32::min(iy1, self._h) {
            for x in u32::max(ix0, 0)..u32::min(ix1, self._h) {
                let l = length(
                    (2.0 * (x as f64 + 0.5) * self._hx - (x0 + x1)) / (x1 - x0),
                    (2.0 * (y as f64 + 0.5) * self._hx - (y0 + y1)) / (y1 - y0),
                );

                let vi = cubicPulse(l) * v;

                if f64::abs(self._src[(x + y * self._w) as usize]) < f64::abs(vi) {
                    self._src[(x + y * self._w) as usize] = vi;
                }
            }
        }
    }

    fn fillSolidFields(&mut self, bodies: &Vec<Box<dyn SolidBody>>) {
        if bodies.is_empty() {
            return;
        }

        for iy in 0..self._h {
            for ix in 0..self._w {
                let idx = (iy * self._w + ix) as usize;
                let x = (ix as f64 + self._ox - 0.5) * self._hx;
                let y = (iy as f64 + self._oy - 0.5) * self._hx;
                self._phi[idx] = bodies[0].distance(x, y);
                for i in 1..bodies.len() {
                    self._phi[idx] = f64::min(self._phi[idx], bodies[i].distance(x, y));
                }
            }
        }

        for iy in 0..self._h {
            for ix in 0..self._w {
                let idx = (iy * self._w + ix) as usize;
                let x = (ix as f64 + self._ox) * self._hx;
                let y = (iy as f64 + self._oy) * self._hx;

                self._body[idx] = 0;

                let mut d = bodies[0].distance(x, y);

                for i in 0..bodies.len() {
                    let id = bodies[i].distance(x, y);
                    if id < d {
                        self._body[idx] = i as u8;
                        d = id;
                    }
                }

                let idxp = ix + iy * (self._w + 1);

                self._volume[idx] = 1.0
                    - occupancy(
                        self._phi[(idxp) as usize],
                        self._phi[(idxp + 1) as usize],
                        self._phi[(idxp + self._w + 1) as usize],
                        self._phi[(idxp + self._w + 2) as usize],
                    );

                if self._volume[idx] < 0.01 {
                    self._volume[idx] = 0.0;
                }

                bodies[self._body[idx as usize] as usize].distanceNormal(
                    &mut self._normalX[idx as usize],
                    &mut self._normalY[idx as usize],
                    x,
                    y,
                );

                self._cell[idx] = if self._volume[idx as usize] == 0.0 {
                    CELL_SOLID
                } else {
                    CELL_FLUID
                };
            }
        }
    }

    fn fillSolidMask(&mut self) {
        for y in 1..self._h - 1 {
            for x in 1..self._w - 1 {
                let idx = (x + y * self._w) as usize;

                if self._cell[idx] == CELL_FLUID {
                    continue;
                }
                let nx = self._normalX[idx];
                let ny = self._normalY[idx];

                self._mask[idx] = 0;

                if nx != 0.0 && self._cell[idx + f64::signum(nx) as usize] != CELL_FLUID {
                    self._mask[idx] |= 1;
                }
                if ny != 0.0
                    && self._cell[idx + (f64::signum(ny) as u32 * self._w) as usize] != CELL_FLUID
                {
                    self._mask[idx] |= 2;
                }
            }
        }
    }

    fn extrapolateNormal(&self, idx: u32) -> f64 {
        let nx = self._normalX[idx as usize];
        let ny = self._normalY[idx as usize];
        let srcX = self._src[(idx + f64::signum(nx) as u32) as usize];
        let srcY = self._src[(idx + f64::signum(ny) as u32 * self._w) as usize];

        (f64::abs(nx) * srcX + f64::abs(ny) * srcY) / (f64::abs(nx) + f64::abs(ny))
    }

    fn freeNeighbour(&mut self, idx: u32, border: &mut Vec<u32>, mask: u8) {
        self._mask[idx as usize] &= !mask;

        if self._cell[idx as usize] != CELL_FLUID && self._mask[idx as usize] == 0 {
            border.push(idx);
        }
    }

    fn extrapolate(&mut self) {
        self.fillSolidMask();

        let mut border = Vec::new();

        for y in 1..self._h - 1 {
            for x in 1..self._w {
                let idx = x + y * self._w;

                if self._cell[idx as usize] != CELL_FLUID && self._mask[idx as usize] == 0 {
                    border.push(idx);
                }
            }
        }

        while !border.is_empty() {
            let idx = border.pop().unwrap();

            self._src[idx as usize] = self.extrapolateNormal(idx);

            if self._normalX[(idx - 1) as usize] > 0.0 {
                self.freeNeighbour(idx - 1, &mut border, 1);
            }
            if self._normalX[(idx + 1) as usize] < 0.0 {
                self.freeNeighbour(idx + 1, &mut border, 1);
            }
            if self._normalY[(idx - self._w) as usize] > 0.0 {
                self.freeNeighbour(idx - self._w, &mut border, 2);
            }
            if self._normalY[(idx + self._w) as usize] < 0.0 {
                self.freeNeighbour(idx + self._w, &mut border, 2);
            }
        }
    }
}

struct FluidSolver {
    _d: FluidQuantity,
    _t: FluidQuantity,
    _u: FluidQuantity,
    _v: FluidQuantity,
    /* Densities at staggered grid locations */
    _uDensity: Vec<f64>,
    _vDensity: Vec<f64>,
    _w: u32,
    _h: u32,
    _hx: f64,
    _densityAir: f64,
    _densitySoot: f64,
    _diffusion: f64,
    _r: Vec<f64>,
    _p: Vec<f64>,
    _z: Vec<f64>,
    _s: Vec<f64>,
    _precon: Vec<f64>,
    _aDiag: Vec<f64>,
    _aPlusX: Vec<f64>,
    _aPlusY: Vec<f64>,
    _tAmb: f64,
    _g: f64,
    _bodies: Vec<Box<SolidBody>>,
}

impl FluidSolver {
    fn buildRhs(&mut self) {
        let scale = 1.0 / self._hx;

        let cell = self._d.cell();
        let body = self._d.body();

        for y in 0..self._h {
            for x in 0..self._w {
                let idx = (y * self._w + x) as usize;

                if (cell[idx] == CELL_FLUID) {
                    self._r[idx] = -scale
                        * (self._u.volume(x + 1, y) * self._u.at(x + 1, y)
                            - self._u.volume(x, y) * self._u.at(x, y)
                            + self._v.volume(x, y + 1) * self._v.at(x, y + 1)
                            - self._v.volume(x, y) * self._v.at(x, y));
                    let vol = self._d.volume(x, y);
                    if self._bodies.is_empty() {
                        continue;
                    }

                    if (x > 0) {
                        self._r[idx] -= (self._u.volume(x, y) - vol)
                            * self._bodies[body[idx - 1] as usize]
                                .velocityX(x as f64 * self._hx, (y as f64 + 0.5) * self._hx);
                    }
                    if (y > 0) {
                        self._r[idx] -= (self._v.volume(x, y) - vol)
                            * self._bodies[body[idx - self._w as usize] as usize]
                                .velocityY((x as f64 + 0.5) * self._hx, y as f64 * self._hx);
                    }
                    if (x < self._w - 1) {
                        self._r[idx] += (self._u.volume(x + 1, y) - vol)
                            * self._bodies[body[idx + 1] as usize].velocityX(
                                (x as f64 + 1.0) * self._hx,
                                (y as f64 + 0.5) * self._hx,
                            );
                    }
                    if (y < self._h - 1) {
                        self._r[idx] += (self._v.volume(x, y + 1) - vol)
                            * self._bodies[body[idx + self._w as usize] as usize].velocityY(
                                (x as f64 + 0.5) * self._hx,
                                (y as f64 + 1.0) * self._hx,
                            );
                    }
                } else {
                    self._r[idx] = 0.0;
                }
            }
        }
    }

    /* Computes densities at the staggered grid locations as a function of
     * temperature and smoke concentration.
     */
    fn computeDensities(&mut self) {
        let alpha = (self._densitySoot - self._densityAir) / self._densityAir;

        unsafe {
            write_bytes(
                self._uDensity.as_mut_ptr(),
                0,
                ((self._w + 1) * self._h) as usize,
            );
            write_bytes(
                self._vDensity.as_mut_ptr(),
                0,
                (self._w * (self._h + 1)) as usize,
            );
        }

        for y in 0..self._h {
            for x in 0..self._w {
                let mut density = self._densityAir * self._tAmb / self._t.at(x, y)
                    * (1.0 + alpha * self._d.at(x, y));
                density = f64::max(density, 0.05 * self._densityAir); /* Clamp dangerously low densities */
                self._uDensity[self._u.idx(x, y) as usize] += 0.5 * density;
                self._vDensity[self._v.idx(x, y) as usize] += 0.5 * density;
                self._uDensity[self._u.idx(x + 1, y) as usize] += 0.5 * density;
                self._vDensity[self._v.idx(x, y + 1) as usize] += 0.5 * density;
            }
        }
    }

    /* Instead of constant density per cell, the entries must now be modified
     * to account for variable density at individual grid cells.
     */
    fn buildPressureMatrix(&mut self, timestep: f64) {
        let scale = timestep / (self._hx * self._hx);
        let cell = self._d.cell();

        unsafe {
            write_bytes(self._aDiag.as_mut_ptr(), 0, (self._w * self._h) as usize);
            write_bytes(self._aPlusX.as_mut_ptr(), 0, (self._w * self._h) as usize);
            write_bytes(self._aPlusY.as_mut_ptr(), 0, (self._w * self._h) as usize);
        }

        for y in 0..self._h {
            for x in 0..self._w {
                let idx = (y * self._w + x) as usize;

                if (cell[idx] != CELL_FLUID) {
                    continue;
                }

                if (x < self._w - 1 && cell[idx + 1] == CELL_FLUID) {
                    let factor = scale * self._u.volume(x + 1, y)
                        / self._uDensity[self._u.idx(x + 1, y) as usize];
                    self._aDiag[idx] += factor;
                    self._aDiag[idx + 1] += factor;
                    self._aPlusX[idx] = -factor;
                }
                if (y < self._h - 1 && cell[idx + self._w as usize] == CELL_FLUID) {
                    let factor = scale * self._v.volume(x, y + 1)
                        / self._vDensity[self._u.idx(x, y + 1) as usize];
                    self._aDiag[idx] += factor;
                    self._aDiag[idx + self._w as usize] += factor;
                    self._aPlusY[idx] = -factor;
                }
            }
        }
    }

    fn buildHeatDiffusionMatrix(&mut self, timestep: f64) {
        for i in 0..self._w * self._h {
            self._aDiag[i as usize] = 1.0;
        }

        unsafe {
            write_bytes(self._aPlusX.as_mut_ptr(), 0, (self._w * self._h) as usize);
            write_bytes(self._aPlusY.as_mut_ptr(), 0, (self._w * self._h) as usize);
        }

        let cell = self._d.cell();
        let scale = self._diffusion * timestep * 1.0 / (self._hx * self._hx);

        for y in 0..self._h {
            for x in 0..self._w {
                let idx = (y * self._w + x) as usize;

                if (cell[idx] != CELL_FLUID) {
                    continue;
                }
                if (x < self._w - 1 && cell[idx + 1] == CELL_FLUID) {
                    self._aDiag[idx] += scale;
                    self._aDiag[idx + 1] += scale;
                    self._aPlusX[idx] = -scale;
                }
                if (y < self._h - 1 && cell[idx + self._w as usize] == CELL_FLUID) {
                    self._aDiag[idx] += scale;
                    self._aDiag[idx + self._w as usize] += scale;
                    self._aPlusY[idx] = -scale;
                }
            }
        }
    }

    fn buildPreconditioner(&mut self) {
        const tau: f64 = 0.97;
        const sigma: f64 = 0.25;

        let cell = self._d.cell();
        for y in 0..self._h {
            for x in 0..self._w {
                let idx = (y * self._w + x) as usize;

                if (cell[idx] != CELL_FLUID) {
                    continue;
                }

                let mut e = self._aDiag[idx];

                if (x > 0 && cell[idx - 1] == CELL_FLUID) {
                    let px = self._aPlusX[idx - 1] * self._precon[idx - 1];
                    let py = self._aPlusY[idx - 1] * self._precon[idx - 1];
                    e = e - (px * px + tau * px * py);
                }
                if (y > 0 && cell[idx - self._w as usize] == CELL_FLUID) {
                    let px =
                        self._aPlusX[idx - self._w as usize] * self._precon[idx - self._w as usize];
                    let py =
                        self._aPlusY[idx - self._w as usize] * self._precon[idx - self._w as usize];
                    e = e - (py * py + tau * px * py);
                }

                if (e < sigma * self._aDiag[idx]) {
                    e = self._aDiag[idx];
                }

                self._precon[idx] = 1.0 / f64::sqrt(e);
            }
        }
    }

    fn applyPreconditioner(&self, dst: &mut Vec<f64>, a: &Vec<f64>) {
        let cell = self._d.cell();

        for y in 0..self._h {
            for x in 0..self._w {
                let idx = (y * self._w + x) as usize;

                if (cell[idx] != CELL_FLUID) {
                    continue;
                }

                let mut t = a[idx];

                if (x > 0 && cell[idx - 1] == CELL_FLUID) {
                    t -= self._aPlusX[idx - 1] * self._precon[idx - 1] * dst[idx - 1];
                }
                if (y > 0 && cell[idx - self._w as usize] == CELL_FLUID) {
                    t -= self._aPlusY[idx - self._w as usize]
                        * self._precon[idx - self._w as usize]
                        * dst[idx - self._w as usize];
                }

                dst[idx] = t * self._precon[idx];
            }
        }

        // TODO: Might bug
        for y in self._h - 1..=0 {
            for x in self._w - 1..=0 {
                let idx = (y * self._w + x) as usize;

                if (cell[idx] != CELL_FLUID) {
                    continue;
                }

                let mut t = dst[idx];

                if (x < self._w - 1 && cell[idx + 1] == CELL_FLUID) {
                    t -= self._aPlusX[idx] * self._precon[idx] * dst[idx + 1];
                }
                if (y < self._h - 1 && cell[idx + self._w as usize] == CELL_FLUID) {
                    t -= self._aPlusY[idx] * self._precon[idx] * dst[idx + self._w as usize];
                }

                dst[idx] = t * self._precon[idx];
            }
        }
    }

    fn dotProduct(&self, a: &Vec<f64>, b: &Vec<f64>) -> f64 {
        let cell = self._d.cell();
        let mut result = 0.0;

        for i in 0..(self._w * self._h) as usize {
            if (cell[i] == CELL_FLUID) {
                result += a[i] * b[i];
            }
        }

        result
    }

    fn matrixVectorProduct(&self, dst: &mut Vec<f64>, b: &Vec<f64>) {
        for y in 0..self._h {
            for x in 0..self._w {
                let idx = (y * self._w + x) as usize;
                let mut t = self._aDiag[idx] * b[idx];
                if (x > 0) {
                    t += self._aPlusX[idx - 1] * b[idx - 1];
                }
                if (y > 0) {
                    t += self._aPlusY[idx - self._w as usize] * b[idx - self._w as usize];
                }
                if (x < self._w - 1) {
                    t += self._aPlusX[idx] * b[idx + 1];
                }
                if (y < self._h - 1) {
                    t += self._aPlusY[idx] * b[idx + self._w as usize];
                }

                dst[idx] = t;
            }
        }
    }

    fn scaledAdd(&self, dst: &mut Vec<f64>, a: &Vec<f64>, b: &Vec<f64>, s: f64) {
        let cell = self._d.cell();
        for i in 0..(self._w * self._h) as usize {
            if (cell[i] == CELL_FLUID) {
                dst[i] = a[i] + b[i] * s;
            }
        }
    }

    fn scaledAddInPlaceA(&self, dst: &mut Vec<f64>, b: &Vec<f64>, s: f64) {
        let cell = self._d.cell();
        for i in 0..(self._w * self._h) as usize {
            if (cell[i] == CELL_FLUID) {
                dst[i] = dst[i] + b[i] * s;
            }
        }
    }

    fn scaledAddInPlaceB(&self, dst: &mut Vec<f64>, a: &Vec<f64>, s: f64) {
        let cell = self._d.cell();
        for i in 0..(self._w * self._h) as usize {
            if (cell[i] == CELL_FLUID) {
                dst[i] = a[i] + dst[i] * s;
            }
        }
    }

    fn infinityNorm(&self, a: &Vec<f64>) -> f64 {
        let cell = self._d.cell();
        let mut maxA = 0.0;

        for i in 0..(self._w * self._h) as usize {
            if (cell[i] == CELL_FLUID) {
                maxA = f64::max(maxA, f64::abs(a[i]));
            }
        }

        maxA
    }

    fn project(&mut self, limit: usize) {
        unsafe {
            write_bytes(self._p.as_mut_ptr(), 0, (self._w * self._h) as usize);
        }

        let mut z = &mut self._z.clone();
        let mut p = &mut self._p.clone();
        let mut r = &mut self._r.clone();
        let mut s = &mut self._s.clone();

        self.applyPreconditioner(&mut z, r);

        unsafe {
            copy_nonoverlapping(z.as_mut_ptr(), s.as_mut_ptr(), (self._w * self._h) as usize);
        }

        let mut maxError = self.infinityNorm(r);

        if (maxError < 1e-5) {
            println!("Initial guess sufficiently small");
            self._z = z.clone();
            self._p = p.clone();
            self._r = r.clone();
            self._s = s.clone();
            return;
        }

        let mut sigma = self.dotProduct(z, r);

        for iter in 0..limit {
            self.matrixVectorProduct(&mut z, s);
            let alpha = sigma / self.dotProduct(z, s);

            self.scaledAddInPlaceA(&mut p, s, alpha);
            self.scaledAddInPlaceA(&mut r, z, -alpha);

            maxError = self.infinityNorm(r);

            if (maxError < 1e-5) {
                println!(
                    "Exiting solver after {} iterations, maximum error is {}",
                    iter, maxError
                );
                self._z = z.clone();
                self._p = p.clone();
                self._r = r.clone();
                self._s = s.clone();
                return;
            }

            self.applyPreconditioner(&mut z, r);
            let sigmaNew = self.dotProduct(z, r);
            self.scaledAddInPlaceB(&mut s, z, sigmaNew / sigma);
            sigma = sigmaNew;
        }

        self._z = z.clone();
        self._p = p.clone();
        self._r = r.clone();
        self._s = s.clone();
        println!(
            "Exceeded budget of {} iterations, maximum error was {}",
            limit, maxError
        );
    }

    /* Similar to the pressure matrix, we cannot assume constant density per
     * cell here either and must modify the equations accordingly.
     */
    fn applyPressure(&mut self, timestep: f64) {
        let scale = timestep / self._hx;
        let cell = self._d.cell();

        for y in 0..self._h {
            for x in 0..self._w {
                let idx = (y * self._w + x) as usize;
                if (cell[idx] != CELL_FLUID) {
                    continue;
                }
                *self._u.at_mut(x, y) -=
                    scale * self._p[idx] / self._uDensity[self._u.idx(x, y) as usize];
                *self._v.at_mut(x, y) -=
                    scale * self._p[idx] / self._vDensity[self._v.idx(x, y) as usize];
                *self._u.at_mut(x + 1, y) +=
                    scale * self._p[idx] / self._uDensity[self._u.idx(x + 1, y) as usize];
                *self._v.at_mut(x, y + 1) +=
                    scale * self._p[idx] / self._vDensity[self._v.idx(x, y + 1) as usize];
            }
        }
    }

    fn addBuoyancy(&mut self, timestep: f64) {
        let alpha = (self._densitySoot - self._densityAir) / self._densityAir;

        for y in 0..self._h {
            for x in 0..self._w {
                let buoyancy = timestep
                    * self._g
                    * (alpha * self._d.at(x, y) - (self._t.at(x, y) - self._tAmb) / self._tAmb);

                *self._v.at_mut(x, y) += buoyancy * 0.5;
                *self._v.at_mut(x, y + 1) += buoyancy * 0.5;
            }
        }
    }

    fn setBoundaryCondition(&mut self) {
        let cell = self._d.cell();
        let body = self._d.body();

        for y in 0..self._h {
            for x in 0..self._w {
                let idx = (y * self._w + x) as usize;
                if (cell[idx] == CELL_SOLID) {
                    let b = &self._bodies[body[idx] as usize];
                    *self._u.at_mut(x, y) =
                        b.velocityX(x as f64 * self._hx, (y as f64 + 0.5) * self._hx);
                    *self._v.at_mut(x, y) =
                        b.velocityY((x as f64 + 0.5) * self._hx, y as f64 * self._hx);
                    *self._u.at_mut(x + 1, y) =
                        b.velocityX((x as f64 + 1.0) * self._hx, (y as f64 + 0.5) * self._hx);
                    *self._v.at_mut(x, y + 1) =
                        b.velocityY((x as f64 + 0.5) * self._hx, (y as f64 + 1.0) * self._hx);
                }
            }
        }
        for y in 0..self._h {
            *self._u.at_mut(0, y) = 0.0;
            *self._u.at_mut(self._w, y) = 0.0;
        }
        for x in 0..self._w {
            *self._v.at_mut(x, 0) = 0.0;
            *self._v.at_mut(x, self._h) = 0.0;
        }
    }

    fn new(
        w: usize,
        h: usize,
        rhoAir: f64,
        rhoSoot: f64,
        diffusion: f64,
        bodies: Vec<Box<SolidBody>>,
    ) -> Self {
        let _hx = 1.0 / usize::min(w, h) as f64;
        let mut ret = Self {
            _w: w as u32,
            _h: h as u32,
            _densityAir: rhoAir,
            _densitySoot: rhoSoot,
            _diffusion: diffusion,
            _bodies: bodies,
            _tAmb: 294.0,
            _g: 9.81,
            _hx,
            _d: FluidQuantity::new(w as u32, h as u32, 0.5, 0.5, _hx),
            _t: FluidQuantity::new(w as u32, h as u32, 0.5, 0.5, _hx),
            _u: FluidQuantity::new(w as u32 + 1, h as u32, 0.0, 0.5, _hx),
            _v: FluidQuantity::new(w as u32, h as u32 + 1, 0.5, 0.0, _hx),
            _r: vec![0f64; w * h],
            _p: vec![0f64; w * h],
            _z: vec![0f64; w * h],
            _s: vec![0f64; w * h],
            _aDiag: vec![0f64; w * h],
            _aPlusX: vec![0f64; w * h],
            _aPlusY: vec![0f64; w * h],
            _precon: vec![0f64; w * h],
            _uDensity: vec![0f64; (w + 1) * h],
            _vDensity: vec![0f64; w * (h + 1)],
        };

        for i in 0..w * h {
            ret._t.src_mut()[i] = ret._tAmb;
        }

        ret
    }

    fn update(&mut self, timestep: f64) {
        self._d.fillSolidFields(&self._bodies);
        self._t.fillSolidFields(&self._bodies);
        self._u.fillSolidFields(&self._bodies);
        self._v.fillSolidFields(&self._bodies);
        unsafe {
            copy_nonoverlapping(
                self._t.src_mut().as_mut_ptr(),
                self._r.as_mut_ptr(),
                (self._w * self._h) as usize,
            );
        }
        self.buildHeatDiffusionMatrix(timestep);
        self.buildPreconditioner();
        self.project(2000);
        unsafe {
            copy_nonoverlapping(
                self._p.as_mut_ptr(),
                self._t.src_mut().as_mut_ptr(),
                (self._w * self._h) as usize,
            );
        }
        self._t.extrapolate();
        self.addBuoyancy(timestep);
        self.setBoundaryCondition();
        self.buildRhs();
        self.computeDensities();
        self.buildPressureMatrix(timestep);
        self.buildPreconditioner();
        self.project(2000);
        self.applyPressure(timestep);
        self._d.extrapolate();
        self._u.extrapolate();
        self._v.extrapolate();
        self.setBoundaryCondition();
        self._d.advect(timestep, &self._u, &self._v, &self._bodies);
        self._t.advect(timestep, &self._u, &self._v, &self._bodies);
        self._u
            .advect(timestep, &self._u.clone(), &self._v, &self._bodies);
        self._v
            .advect(timestep, &self._u, &self._v.clone(), &self._bodies);
        self._d.flip();
        self._t.flip();
        self._u.flip();
        self._v.flip();
    }

    fn addInflow(&mut self, x: f64, y: f64, w: f64, h: f64, d: f64, t: f64, u: f64, v: f64) {
        self._d.addInflow(x, y, x + w, y + h, d);
        self._t.addInflow(x, y, x + w, y + h, t);
        self._u.addInflow(x, y, x + w, y + h, u);
        self._v.addInflow(x, y, x + w, y + h, v);
    }

    fn ambientT(&self) -> f64 {
        self._tAmb
    }

    fn toImage(&self, renderHeat: bool) -> Vec<u8> {
        let mut rgba = vec![0u8; (self._h * if renderHeat { 2 } else { 1 } * self._w * 4) as usize];

        for y in 0..self._h {
            for x in 0..self._w {
                let (idxl, idxr) = if (renderHeat) {
                    (
                        4 * (x + y * self._w * 2) as usize,
                        4 * (x + y * self._w * 2 + self._w) as usize,
                    )
                } else {
                    (0, 4 * (x + y * self._w) as usize)
                };

                let volume = self._d.volume(x, y);
                let mut shade = (1.0 - self._d.at(x, y)) * volume;

                shade = f64::min(f64::max(shade, 0.0), 1.0);

                rgba[idxr + 0] = (shade * 255.0) as u8;
                rgba[idxr + 1] = (shade * 255.0) as u8;
                rgba[idxr + 2] = (shade * 255.0) as u8;
                rgba[idxr + 3] = 0xFF;

                if (renderHeat) {
                    let mut t = (self._t.at(x, y) - self._tAmb) / 700.0;
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
    let sizeX = 256;
    let sizeY = 256;
    let densityAir = 0.1;
    let densitySoot = 1.0; /* You can make this smaller to get lighter smoke */
    let diffusion = 0.03;
    let timestep = 0.005;
    let renderHeat = true; /* Set this to true to enable heat rendering */

    let mut bodies: Vec<Box<SolidBody>> = Vec::new();

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

    let mut solver = FluidSolver::new(sizeX, sizeY, densityAir, densitySoot, diffusion, bodies);
    let mut time = 0.0;
    let mut iterations = 0;

    while (time < 8.0) {
        for i in 0..4 {
            solver.addInflow(
                0.45,
                0.2,
                0.1,
                0.05,
                1.0,
                650.0 + solver.ambientT(),
                0.0,
                0.0,
            );
            solver.update(timestep);
            time += timestep;
        }
        let image = solver.toImage(renderHeat);
        let path = format!("frame{:0>5}.png", iterations);

        iterations += 1;

        screenshot(
            &path,
            if renderHeat { sizeX * 2 } else { sizeX } as u32,
            sizeY as u32,
            image,
        );
        // lodepng_encode32_file(path, image, (renderHeat ? sizeX*2 : sizeX), sizeY);

        for i in 0..solver._bodies.len() {
            solver._bodies[i].update(timestep);
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
