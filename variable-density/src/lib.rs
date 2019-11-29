// #define _USE_MATH_DEFINES

// #include <algorithm>
// #include <stdint.h>
// #include <stdio.h>
// #include <vector>
// #include <math.h>
// #include <stack>

// #include "../lodepng/lodepng.h"

// using namespace std;

// template <typename T> int f64::signum(T val) {
//     return (T(0) < val) - (val < T(0));
// }

// template <typename T> int nf64::signum(T val) {
//     return (val < T(0) ? -1 : 1);
// }

use std::mem;

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
    *y = -f64::sin(phi) * tmpX + f64::cos(phi) * tmpY;
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
        let x0 = u32::max(ix - 1, 0);
        let x1 = ix;
        let x2 = ix + 1;
        let x3 = u32::min(ix + 2, self._w - 1);
        let y0 = u32::max(iy - 1, 0);
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

    // fn advect(timestep:f64, u:& FluidQuantity, v:&FluidQuantity,      bodies: &Vec<SolidBody>) {
    //         let idx = 0;

    //         for iy in 0.._h {
    //  for ix in 0.._w {
    //      idx +=1;
    //     for (int iy = 0, idx = 0; iy < _h; iy++) {
    //         for (int ix = 0; ix < _w; ix++, idx++) {
    //             if (_cell[idx] == CELL_FLUID) {
    //                 double x = ix + _ox;
    //                 double y = iy + _oy;
    //                 rungeKutta3(x, y, timestep, u, v);
    //                 backProject(x, y, bodies);
    //                 _dst[idx] = cerp(x, y);
    //             }
    //         }
    //     }
    // }

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

        let mut idx = 0usize;

        for iy in 0..self._h {
            for ix in 0..self._w {
                idx += 1;
                let x = (ix as f64 + self._ox - 0.5) * self._hx;
                let y = (iy as f64 + self._oy - 0.5) * self._hx;
                self._phi[idx] = bodies[0].distance(x, y);

                for i in 1..bodies.len() {
                    self._phi[idx] = f64::min(self._phi[idx], bodies[i].distance(x, y));
                }
            }
        }

        idx = 0;
        for iy in 0..self._h {
            for ix in 0..self._w {
                idx += 1;

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

// class FluidSolver {
//     FluidQuantity *_d;
//     FluidQuantity *_t;
//     FluidQuantity *_u;
//     FluidQuantity *_v;
//     /* Densities at staggered grid locations */
//     double *_uDensity;
//     double *_vDensity;
//     int _w;
//     int _h;
//     double _hx;
//     double _densityAir;
//     double _densitySoot;
//     double _diffusion;
//     double *_r;
//     double *_p;
//     double *_z;
//     double *_s;
//     double *_precon;
//     double *_aDiag;
//     double *_aPlusX;
//     double *_aPlusY;
//     double _tAmb;
//     double _g;
//     const vector<const SolidBody *> &_bodies;
//     void buildRhs() {
//         double scale = 1.0/_hx;
//         const uint8_t *cell = _d->cell();
//         const uint8_t *body = _d->body();
//         for (int y = 0, idx = 0; y < _h; y++) {
//             for (int x = 0; x < _w; x++, idx++) {
//                 if (cell[idx] == CELL_FLUID) {
//                     _r[idx] = -scale*
//                         (_u->volume(x + 1, y)*_u->at(x + 1, y) - _u->volume(x, y)*_u->at(x, y) +
//                          _v->volume(x, y + 1)*_v->at(x, y + 1) - _v->volume(x, y)*_v->at(x, y));
//                     double vol = _d->volume(x, y);
//                     if (_bodies.empty())
//                         continue;
//                     if (x > 0)
//                         _r[idx] -= (_u->volume(x, y) - vol)*_bodies[body[idx -  1]]->velocityX(x*_hx, (y + 0.5)*_hx);
//                     if (y > 0)
//                         _r[idx] -= (_v->volume(x, y) - vol)*_bodies[body[idx - _w]]->velocityY((x + 0.5)*_hx, y*_hx);
//                     if (x < _w - 1)
//                         _r[idx] += (_u->volume(x + 1, y) - vol)*_bodies[body[idx +  1]]->velocityX((x + 1.0)*_hx, (y + 0.5)*_hx);
//                     if (y < _h - 1)
//                         _r[idx] += (_v->volume(x, y + 1) - vol)*_bodies[body[idx + _w]]->velocityY((x + 0.5)*_hx, (y + 1.0)*_hx);
//                 } else
//                     _r[idx] = 0.0;
//             }
//         }
//     }
//     /* Computes densities at the staggered grid locations as a function of
//      * temperature and smoke concentration.
//      */
//     void computeDensities() {
//         double alpha = (_densitySoot - _densityAir)/_densityAir;
//         memset(_uDensity, 0, (_w + 1)*_h*sizeof(double));
//         memset(_vDensity, 0, _w*(_h + 1)*sizeof(double));
//         for (int y = 0; y < _h; y++) {
//             for (int x = 0; x < _w; x++) {
//                 double density = _densityAir*_tAmb/_t->at(x, y)*(1.0 + alpha*_d->at(x, y));
//                 density = max(density, 0.05*_densityAir); /* Clamp dangerously low densities */
//                 _uDensity[_u->idx(x, y)]     += 0.5*density;
//                 _vDensity[_v->idx(x, y)]     += 0.5*density;
//                 _uDensity[_u->idx(x + 1, y)] += 0.5*density;
//                 _vDensity[_v->idx(x, y + 1)] += 0.5*density;
//             }
//         }
//     }
//     /* Instead of constant density per cell, the entries must now be modified
//      * to account for variable density at individual grid cells.
//      */
//     void buildPressureMatrix(double timestep) {
//         double scale = timestep/(_hx*_hx);
//         const uint8_t *cell = _d->cell();
//         memset(_aDiag,  0, _w*_h*sizeof(double));
//         memset(_aPlusX, 0, _w*_h*sizeof(double));
//         memset(_aPlusY, 0, _w*_h*sizeof(double));

//         for (int y = 0, idx = 0; y < _h; y++) {
//             for (int x = 0; x < _w; x++, idx++) {
//                 if (cell[idx] != CELL_FLUID)
//                     continue;
//                 if (x < _w - 1 && cell[idx + 1] == CELL_FLUID) {
//                     double factor = scale*_u->volume(x + 1, y)/_uDensity[_u->idx(x + 1, y)];
//                     _aDiag [idx    ] +=  factor;
//                     _aDiag [idx + 1] +=  factor;
//                     _aPlusX[idx    ]  = -factor;
//                 }
//                 if (y < _h - 1 && cell[idx + _w] == CELL_FLUID) {
//                     double factor = scale*_v->volume(x, y + 1)/_vDensity[_u->idx(x, y + 1)];
//                     _aDiag [idx     ] +=  factor;
//                     _aDiag [idx + _w] +=  factor;
//                     _aPlusY[idx     ]  = -factor;
//                 }
//             }
//         }
//     }
//     void buildHeatDiffusionMatrix(double timestep) {
//         for (int i = 0; i < _w*_h; i++)
//             _aDiag[i] = 1.0;
//         memset(_aPlusX, 0, _w*_h*sizeof(double));
//         memset(_aPlusY, 0, _w*_h*sizeof(double));

//         const uint8_t *cell = _d->cell();
//         double scale = _diffusion*timestep*1.0/(_hx*_hx);

//         for (int y = 0, idx = 0; y < _h; y++) {
//             for (int x = 0; x < _w; x++, idx++) {
//                 if (cell[idx] != CELL_FLUID)
//                     continue;
//                 if (x < _w - 1 && cell[idx + 1] == CELL_FLUID) {
//                     _aDiag [idx    ] +=  scale;
//                     _aDiag [idx + 1] +=  scale;
//                     _aPlusX[idx    ]  = -scale;
//                 }

//                 if (y < _h - 1 && cell[idx + _w] == CELL_FLUID) {
//                     _aDiag [idx     ] +=  scale;
//                     _aDiag [idx + _w] +=  scale;
//                     _aPlusY[idx     ]  = -scale;
//                 }
//             }
//         }
//     }
//     void buildPreconditioner() {
//         const double tau = 0.97;
//         const double sigma = 0.25;
//         const uint8_t *cell = _d->cell();

//         for (int y = 0, idx = 0; y < _h; y++) {
//             for (int x = 0; x < _w; x++, idx++) {
//                 if (cell[idx] != CELL_FLUID)
//                     continue;
//                 double e = _aDiag[idx];

//                 if (x > 0 && cell[idx - 1] == CELL_FLUID) {
//                     double px = _aPlusX[idx - 1]*_precon[idx - 1];
//                     double py = _aPlusY[idx - 1]*_precon[idx - 1];
//                     e = e - (px*px + tau*px*py);
//                 }
//                 if (y > 0 && cell[idx - _w] == CELL_FLUID) {
//                     double px = _aPlusX[idx - _w]*_precon[idx - _w];
//                     double py = _aPlusY[idx - _w]*_precon[idx - _w];
//                     e = e - (py*py + tau*px*py);
//                 }

//                 if (e < sigma*_aDiag[idx])
//                     e = _aDiag[idx];

//                 _precon[idx] = 1.0/sqrt(e);
//             }
//         }
//     }
//     void applyPreconditioner(double *dst, double *a) {
//         const uint8_t *cell = _d->cell();
//         for (int y = 0, idx = 0; y < _h; y++) {
//             for (int x = 0; x < _w; x++, idx++) {
//                 if (cell[idx] != CELL_FLUID)
//                     continue;
//                 double t = a[idx];

//                 if (x > 0 && cell[idx -  1] == CELL_FLUID)
//                     t -= _aPlusX[idx -  1]*_precon[idx -  1]*dst[idx -  1];
//                 if (y > 0 && cell[idx - _w] == CELL_FLUID)
//                     t -= _aPlusY[idx - _w]*_precon[idx - _w]*dst[idx - _w];

//                 dst[idx] = t*_precon[idx];
//             }
//         }

//         for (int y = _h - 1, idx = _w*_h - 1; y >= 0; y--) {
//             for (int x = _w - 1; x >= 0; x--, idx--) {
//                 if (cell[idx] != CELL_FLUID)
//                     continue;
//                 double t = dst[idx];

//                 if (x < _w - 1 && cell[idx +  1] == CELL_FLUID)
//                     t -= _aPlusX[idx]*_precon[idx]*dst[idx +  1];
//                 if (y < _h - 1 && cell[idx + _w] == CELL_FLUID)
//                     t -= _aPlusY[idx]*_precon[idx]*dst[idx + _w];

//                 dst[idx] = t*_precon[idx];
//             }
//         }
//     }
//     double dotProduct(double *a, double *b) {
//         const uint8_t *cell = _d->cell();
//         double result = 0.0;
//         for (int i = 0; i < _w*_h; i++)
//             if (cell[i] == CELL_FLUID)
//                 result += a[i]*b[i];
//         return result;
//     }
//     void matrixVectorProduct(double *dst, double *b) {
//         for (int y = 0, idx = 0; y < _h; y++) {
//             for (int x = 0; x < _w; x++, idx++) {
//                 double t = _aDiag[idx]*b[idx];
//                 if (x > 0)
//                     t += _aPlusX[idx -  1]*b[idx -  1];
//                 if (y > 0)
//                     t += _aPlusY[idx - _w]*b[idx - _w];
//                 if (x < _w - 1)
//                     t += _aPlusX[idx]*b[idx +  1];
//                 if (y < _h - 1)
//                     t += _aPlusY[idx]*b[idx + _w];

//                 dst[idx] = t;
//             }
//         }
//     }
//     void scaledAdd(double *dst, double *a, double *b, double s) {
//         const uint8_t *cell = _d->cell();
//         for (int i = 0; i < _w*_h; i++)
//             if (cell[i] == CELL_FLUID)
//                 dst[i] = a[i] + b[i]*s;
//     }
//     double infinityNorm(double *a) {
//         const uint8_t *cell = _d->cell();
//         double maxA = 0.0;
//         for (int i = 0; i < _w*_h; i++)
//             if (cell[i] == CELL_FLUID)
//                 maxA = max(maxA, fabs(a[i]));
//         return maxA;
//     }
//     void project(int limit) {
//         memset(_p, 0,  _w*_h*sizeof(double));
//         applyPreconditioner(_z, _r);
//         memcpy(_s, _z, _w*_h*sizeof(double));
//         double maxError = infinityNorm(_r);
//         if (maxError < 1e-5) {
//             printf("Initial guess sufficiently small\n");
//             return;
//         }
//         double sigma = dotProduct(_z, _r);
//         for (int iter = 0; iter < limit; iter++) {
//             matrixVectorProduct(_z, _s);
//             double alpha = sigma/dotProduct(_z, _s);
//             scaledAdd(_p, _p, _s, alpha);
//             scaledAdd(_r, _r, _z, -alpha);
//             maxError = infinityNorm(_r);
//             if (maxError < 1e-5) {
//                 printf("Exiting solver after %d iterations, maximum error is %f\n", iter, maxError);
//                 return;
//             }
//             applyPreconditioner(_z, _r);
//             double sigmaNew = dotProduct(_z, _r);
//             scaledAdd(_s, _z, _s, sigmaNew/sigma);
//             sigma = sigmaNew;
//         }
//         printf("Exceeded budget of %d iterations, maximum error was %f\n", limit, maxError);
//     }
//     /* Similar to the pressure matrix, we cannot assume constant density per
//      * cell here either and must modify the equations accordingly.
//      */
//     void applyPressure(double timestep) {
//         double scale = timestep/_hx;
//         const uint8_t *cell = _d->cell();
//         for (int y = 0, idx = 0; y < _h; y++) {
//             for (int x = 0; x < _w; x++, idx++) {
//                 if (cell[idx] != CELL_FLUID)
//                     continue;
//                 _u->at(x, y)     -= scale*_p[idx]/_uDensity[_u->idx(x, y)];
//                 _v->at(x, y)     -= scale*_p[idx]/_vDensity[_v->idx(x, y)];
//                 _u->at(x + 1, y) += scale*_p[idx]/_uDensity[_u->idx(x + 1, y)];
//                 _v->at(x, y + 1) += scale*_p[idx]/_vDensity[_v->idx(x, y + 1)];
//             }
//         }
//     }
//     void addBuoyancy(double timestep) {
//         double alpha = (_densitySoot - _densityAir)/_densityAir;

//         for (int y = 0; y < _h; y++) {
//             for (int x = 0; x < _w; x++) {
//                 double buoyancy = timestep*_g*(alpha*_d->at(x, y) - (_t->at(x, y) - _tAmb)/_tAmb);

//                 _v->at(x, y    ) += buoyancy*0.5;
//                 _v->at(x, y + 1) += buoyancy*0.5;
//             }
//         }
//     }
//     void setBoundaryCondition() {
//         const uint8_t *cell = _d->cell();
//         const uint8_t *body = _d->body();
//         for (int y = 0, idx = 0; y < _h; y++) {
//             for (int x = 0; x < _w; x++, idx++) {
//                 if (cell[idx] == CELL_SOLID) {
//                     const SolidBody &b = *_bodies[body[idx]];
//                     _u->at(x, y) = b.velocityX(x*_hx, (y + 0.5)*_hx);
//                     _v->at(x, y) = b.velocityY((x + 0.5)*_hx, y*_hx);
//                     _u->at(x + 1, y) = b.velocityX((x + 1.0)*_hx, (y + 0.5)*_hx);
//                     _v->at(x, y + 1) = b.velocityY((x + 0.5)*_hx, (y + 1.0)*_hx);
//                 }
//             }
//         }
//         for (int y = 0; y < _h; y++)
//             _u->at(0, y) = _u->at(_w, y) = 0.0;
//         for (int x = 0; x < _w; x++)
//             _v->at(x, 0) = _v->at(x, _h) = 0.0;
//     }
// public:
//     FluidSolver(int w, int h, double rhoAir, double rhoSoot, double diffusion,
//             const vector<const SolidBody *> &bodies) : _w(w), _h(h),
//             _densityAir(rhoAir), _densitySoot(rhoSoot), _diffusion(diffusion),
//             _bodies(bodies) {
//         _tAmb = 294.0;
//         _g    = 9.81;
//         _hx = 1.0/min(w, h);
//         _d = new FluidQuantity(_w,     _h,     0.5, 0.5, _hx);
//         _t = new FluidQuantity(_w,     _h,     0.5, 0.5, _hx);
//         _u = new FluidQuantity(_w + 1, _h,     0.0, 0.5, _hx);
//         _v = new FluidQuantity(_w,     _h + 1, 0.5, 0.0, _hx);
//         for (int i = 0; i < _w*_h; i++)
//             _t->src()[i] = _tAmb;
//         _r = new double[_w*_h];
//         _p = new double[_w*_h];
//         _z = new double[_w*_h];
//         _s = new double[_w*_h];
//         _aDiag  = new double[_w*_h];
//         _aPlusX = new double[_w*_h];
//         _aPlusY = new double[_w*_h];
//         _precon = new double[_w*_h];
//         _uDensity = new double[(_w + 1)*_h];
//         _vDensity = new double[_w*(_h + 1)];
//     }
//     ~FluidSolver() {
//         delete _d;
//         delete _t;
//         delete _u;
//         delete _v;
//         delete[] _r;
//         delete[] _p;
//         delete[] _z;
//         delete[] _s;
//         delete[] _aDiag;
//         delete[] _aPlusX;
//         delete[] _aPlusY;
//         delete[] _precon;
//         delete[] _uDensity;
//         delete[] _vDensity;
//     }
//     void update(double timestep) {
//         _d->fillSolidFields(_bodies);
//         _t->fillSolidFields(_bodies);
//         _u->fillSolidFields(_bodies);
//         _v->fillSolidFields(_bodies);
//         memcpy(_r, _t->src(), _w*_h*sizeof(double));
//         buildHeatDiffusionMatrix(timestep);
//         buildPreconditioner();
//         project(2000);
//         memcpy(_t->src(), _p, _w*_h*sizeof(double));
//         _t->extrapolate();
//         addBuoyancy(timestep);
//         setBoundaryCondition();
//         buildRhs();
//         computeDensities();
//         buildPressureMatrix(timestep);
//         buildPreconditioner();
//         project(2000);
//         applyPressure(timestep);
//         _d->extrapolate();
//         _u->extrapolate();
//         _v->extrapolate();
//         setBoundaryCondition();
//         _d->advect(timestep, *_u, *_v, _bodies);
//         _t->advect(timestep, *_u, *_v, _bodies);
//         _u->advect(timestep, *_u, *_v, _bodies);
//         _v->advect(timestep, *_u, *_v, _bodies);
//         _d->flip();
//         _t->flip();
//         _u->flip();
//         _v->flip();
//     }
//     void addInflow(double x, double y, double w, double h, double d, double t, double u, double v) {
//         _d->addInflow(x, y, x + w, y + h, d);
//         _t->addInflow(x, y, x + w, y + h, t);
//         _u->addInflow(x, y, x + w, y + h, u);
//         _v->addInflow(x, y, x + w, y + h, v);
//     }
//     double ambientT() {
//         return _tAmb;
//     }
//     void toImage(unsigned char *rgba, bool renderHeat) {
//         for (int y = 0; y < _h; y++) {
//             for (int x = 0; x < _w; x++) {
//                 int idxl, idxr;
//                 if (renderHeat) {
//                     idxl = 4*(x + y*_w*2);
//                     idxr = 4*(x + y*_w*2 + _w);
//                 } else
//                     idxr = 4*(x + y*_w);
//                 double volume = _d->volume(x, y);
//                 double shade = (1.0 - _d->at(x, y))*volume;
//                 shade = min(max(shade, 0.0), 1.0);
//                 rgba[idxr + 0] = (int)(shade*255.0);
//                 rgba[idxr + 1] = (int)(shade*255.0);
//                 rgba[idxr + 2] = (int)(shade*255.0);
//                 rgba[idxr + 3] = 0xFF;
//                 if (renderHeat) {
//                     double t = (_t->at(x, y) - _tAmb)/700.0;
//                     t = min(max(t, 0.0), 1.0);
//                     double r = 1.0 + volume*(min(t*4.0, 1.0) - 1.0);
//                     double g = 1.0 + volume*(min(t*2.0, 1.0) - 1.0);
//                     double b = 1.0 + volume*(max(min(t*4.0 - 3.0, 1.0), 0.0) - 1.0);
//                     rgba[idxl + 0] = (int)(r*255.0);
//                     rgba[idxl + 1] = (int)(g*255.0);
//                     rgba[idxl + 2] = (int)(b*255.0);
//                     rgba[idxl + 3] = 0xFF;
//                 }
//             }
//         }
//     }
// };
