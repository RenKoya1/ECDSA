use num_bigint::BigUint;

#[derive(PartialEq, Clone, Debug)]
pub enum Point {
    Coor(BigUint, BigUint),
    Identity,
}
pub struct EllipticCurve {
    // y^1 = x^2 + a*x + b
    pub a: BigUint,
    pub b: BigUint,
    pub p: BigUint,
}

impl EllipticCurve {
    pub fn add(&self, c: &Point, d: &Point) -> Point {
        assert!(self.is_on_curve(c), "p must be on the curve");
        assert!(self.is_on_curve(d), "q must be on the curve");
        assert!(*c != *d, "p and q must be different");
        match (c, d) {
            (Point::Identity, d) => d.clone(),
            (c, Point::Identity) => c.clone(),
            (Point::Coor(x1, y1), Point::Coor(x2, y2)) => {
                let y1plusy2 = FiniteField::add(y1, y2, &self.p);
                if x1 == x2 && y1plusy2 == BigUint::from(0u32) {
                    return Point::Identity;
                }
                // s = (y2 - y1) /( x2 - x1)
                // x3 = s^2 - x1 - x2
                // y3 = s(x1 - x3) - y1
                let y2my1 = FiniteField::substract(y2, y1, &self.p);
                let x2mx1 = FiniteField::substract(x2, x1, &self.p);

                let s = FiniteField::divide(&y2my1, &x2mx1, &self.p);

                let s2 = s.modpow(&BigUint::from(2u32), &self.p);
                let x3 = FiniteField::substract(
                    &FiniteField::substract(&s2, &x1, &self.p),
                    &x2,
                    &self.p,
                );
                let y3 = FiniteField::substract(
                    &FiniteField::mult(&s, &FiniteField::substract(x1, &x3, &self.p), &self.p),
                    y1,
                    &self.p,
                );
                Point::Coor(x3, y3)
            }
        }
    }

    fn double(&self, p: &Point) -> Point {
        assert!(self.is_on_curve(p), "p must be on the curve");
        match p {
            Point::Identity => Point::Identity,
            Point::Coor(x, y) => {
                // s = (3x^2 + a) / (2y)
                // x3 = s^2 - 2x
                // y3 = s(x - x3) - y
                let x2 = x.modpow(&BigUint::from(2u32), &self.p);
                let s_u = FiniteField::mult(&BigUint::from(3u32), &x2, &self.p);
                let s_u = FiniteField::add(&s_u, &self.a, &self.p);

                let s_b = FiniteField::mult(&y, &BigUint::from(2u32), &self.p);

                let s = FiniteField::divide(&s_u, &s_b, &self.p);
                let s2 = s.modpow(&BigUint::from(2u32), &self.p);
                let x3 = FiniteField::substract(
                    &s2,
                    &FiniteField::mult(&x, &BigUint::from(2u32), &self.p),
                    &self.p,
                );
                let y3 = FiniteField::substract(
                    &FiniteField::mult(&s, &FiniteField::substract(&x, &x3, &self.p), &self.p),
                    &y,
                    &self.p,
                );

                Point::Coor(x3, y3)
            }
        }
    }
    pub fn scalar_mul(&self, p: &Point, d: &BigUint) -> Point {
        // double and add algorithm
        // B = d * A
        let mut t = p.clone();
        if *d == BigUint::from(0u32){ 
          return Point::Identity; 
      }
        for i in (0..(d.bits() - 1)).rev() {
            t = self.double(&t);
            if d.bit(i) {
                t = self.add(&t, &p);
            }
        }
        t
    }
    fn is_on_curve(&self, p: &Point) -> bool {
        // y^2 = x^3 + ax + b
        match p {
            Point::Coor(x, y) => {
                let y2 = y.modpow(&BigUint::from(2u32), &self.p);
                let x3 = x.modpow(&BigUint::from(3u32), &self.p);
                let ax = FiniteField::mult(&self.a, x, &self.p);
                y2 == FiniteField::add(&FiniteField::add(&x3, &ax, &self.p), &self.b, &self.p)
            }
            Point::Identity => true,
        }
    }
}

pub struct FiniteField {}

impl FiniteField {
    pub fn add(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint {
        // c + d = r mod p
        let r = c + d;
        r.modpow(&BigUint::from(1u32), p)
    }
    pub fn mult(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint {
        // c * d = r mod p
        let r = c * d;
        r.modpow(&BigUint::from(1u32), p)
    }

    pub fn inv_add(c: &BigUint, p: &BigUint) -> BigUint {
        // -c mod p
        assert!(c < p, "c must be less than p");
        let r = p - c;
        r.modpow(&BigUint::from(1u32), p)
    }

    pub fn substract(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint {
        // -c mod p
        let d_inv = FiniteField::inv_add(d, p);
        FiniteField::add(c, &d_inv, p)
    }

    pub fn inv_mult(c: &BigUint, p: &BigUint) -> BigUint {
        // c^-1 mod p = c^(p-2) mod p
        // Fermat's little theorem
        // c^(p-1) mod p = 1
        c.modpow(&(p - BigUint::from(2u32)), &p)
    }

    pub fn divide(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint {
        let d_inv = FiniteField::inv_mult(d, p);
        FiniteField::mult(c, &d_inv, p)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_add() {
        let c = BigUint::from(4u32);
        let d = BigUint::from(10u32);
        let p = BigUint::from(11u32);
        let r = FiniteField::add(&c, &d, &p);
        assert_eq!(r, BigUint::from(3u32));
    }
    #[test]
    fn test_mult() {
        let c = BigUint::from(4u32);
        let d = BigUint::from(10u32);
        let p = BigUint::from(11u32);
        let r = FiniteField::mult(&c, &d, &p);
        assert_eq!(r, BigUint::from(7u32));
    }
    #[test]
    fn test_inv_add() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(11u32);
        let r = FiniteField::inv_add(&c, &p);
        assert_eq!(
            (c + r).modpow(&BigUint::from(1u32), &p),
            BigUint::from(0u32)
        );
    }

    #[test]
    fn test_inv_mult() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(11u32);
        let r = FiniteField::inv_mult(&c, &p);
        assert_eq!(
            (c * r).modpow(&BigUint::from(1u32), &p),
            BigUint::from(1u32)
        );
    }

    #[test]
    fn test_substract() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(11u32);
        let r = FiniteField::substract(&c, &c, &p);
        assert_eq!(r, BigUint::from(0u32));
    }

    #[test]
    fn test_divide() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(11u32);
        let r = FiniteField::divide(&c, &c, &p);
        assert_eq!(r, BigUint::from(1u32));
    }

    #[test]
    fn test_ec_point_addition() {
        //y^2 = x^3 + 2x + 2 mod 17
        let ec = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            p: BigUint::from(17u32),
        };
        // (6, 3) + (5, 1) = (10, 6)
        let p1 = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));
        let p2 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        let pr = Point::Coor(BigUint::from(10u32), BigUint::from(6u32));
        let res = ec.add(&p1, &p2);
        assert_eq!(res, pr);
    }
    #[test]
    fn test_ec_scalar() {
        //y^2 = x^3 + 2x + 2 mod 17
        let ec = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            p: BigUint::from(17u32),
        };
        // 2(5, 1) = (6, 3)
        let p1 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        let p2 = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));

        let res = ec.scalar_mul(&p1, &BigUint::from(2u32));
        assert_eq!(res, p2);
    }
}
