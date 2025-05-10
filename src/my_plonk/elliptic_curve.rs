use super::field::{FieldElement, FieldElementExt};

#[derive(Clone, Debug, PartialEq)]
pub struct Point {
    pub x: FieldElement,
    pub y: FieldElement,
    pub is_infinity: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PointExt {
    pub x: FieldElementExt,
    pub y: FieldElementExt,
    pub is_infinity: bool,
}

#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct EllipticCurve {
    a: FieldElement,
    b: FieldElement,
    g1: Point,    // Generator for G1
    g2: PointExt, // Generator for G2
    order: u64,   // Subgroup order
}

impl Point {
    #[allow(dead_code)]
    pub fn infinity() -> Self {
        Point {
            x: FieldElement::new(0),
            y: FieldElement::new(0),
            is_infinity: true,
        }
    }

    #[allow(dead_code)]
    pub fn scalar_mul(&self, curve: &EllipticCurve, scalar: FieldElement) -> Point {
        let mut result = Point::infinity();
        let mut temp = self.clone();
        let mut s = scalar.value;
        while s > 0 {
            if s & 1 == 1 {
                result = curve.add(&result, &temp);
            }
            temp = curve.add(&temp, &temp);
            s >>= 1;
        }
        result
    }
}

impl PointExt {
    #[allow(dead_code)]
    pub fn infinity() -> Self {
        PointExt {
            x: FieldElementExt::new(FieldElement::zero(), FieldElement::zero()),
            y: FieldElementExt::new(FieldElement::zero(), FieldElement::zero()),
            is_infinity: true,
        }
    }

    #[allow(dead_code)]
    pub fn scalar_mul(&self, curve: &EllipticCurve, scalar: FieldElement) -> PointExt {
        let mut result = PointExt::infinity();
        let mut temp = self.clone();
        let mut s = scalar.value;
        while s > 0 {
            if s & 1 == 1 {
                result = curve.add_ext(&result, &temp);
            }
            temp = curve.add_ext(&temp, &temp);
            s >>= 1;
        }
        result
    }
}

impl EllipticCurve {
    #[allow(dead_code)]
    pub fn new() -> Self {
        let a = FieldElement::new(0); // y^2 = x^3 + 0x + 3
        let b = FieldElement::new(3);
        let g1 = Point {
            x: FieldElement::new(1),
            y: FieldElement::new(2),
            is_infinity: false,
        };
        let g2 = PointExt {
            x: FieldElementExt::new(FieldElement::new(36), FieldElement::zero()),
            y: FieldElementExt::new(FieldElement::zero(), FieldElement::new(31)), // 31u
            is_infinity: false,
        };
        EllipticCurve {
            a,
            b,
            g1,
            g2,
            order: 17, // there are 17 valid points generated from (1, 2)
        }
    }

    #[allow(dead_code)]
    pub fn generator_g1(&self) -> Point {
        self.g1.clone()
    }

    #[allow(dead_code)]
    pub fn generator_g2(&self) -> PointExt {
        self.g2.clone()
    }

    #[allow(dead_code)]
    pub fn infinity(&self) -> Point {
        Point::infinity()
    }

    #[allow(dead_code)]
    pub fn infinity_ext(&self) -> PointExt {
        PointExt::infinity()
    }

    pub fn add(&self, p1: &Point, p2: &Point) -> Point {
        if p1.is_infinity {
            return p2.clone();
        }
        if p2.is_infinity {
            return p1.clone();
        }
        if p1.x == p2.x && p1.y != p2.y {
            return Point::infinity();
        }

        let m: FieldElement = if p1 == p2 {
            let num: FieldElement =
                p1.x.multiply(&p1.x)
                    .multiply(&FieldElement::new(3))
                    .add(&self.a);
            let den = p1.y.multiply(&FieldElement::new(2));
            num.divide(&den)
        } else {
            let num: FieldElement = p2.y.substract(&p1.y);
            let den: FieldElement = p2.x.substract(&p1.x);
            num.divide(&den)
        };

        let x3: FieldElement = m.multiply(&m).substract(&p1.x).substract(&p2.x);
        let y3: FieldElement = m.multiply(&p1.x.substract(&x3)).substract(&p1.y);
        Point {
            x: x3,
            y: y3,
            is_infinity: false,
        }
    }

    pub fn add_ext(&self, p1: &PointExt, p2: &PointExt) -> PointExt {
        if p1.is_infinity {
            return p2.clone();
        }
        if p2.is_infinity {
            return p1.clone();
        }
        if p1.x == p2.x && p1.y != p2.y {
            return PointExt::infinity();
        }

        let m = if p1 == p2 {
            let num =
                p1.x.multiply(&p1.x)
                    .multiply(&FieldElementExt::new(
                        FieldElement::new(3),
                        FieldElement::zero(),
                    ))
                    .add(&FieldElementExt::new(self.a.clone(), FieldElement::zero()));
            let den = p1.y.multiply(&FieldElementExt::new(
                FieldElement::new(2),
                FieldElement::zero(),
            ));
            // Simplified: This needs proper division in F_101^2
            let m_num = num.a; // Toy approximation
            let m_den = den.a;
            let m_field = m_num.divide(&m_den);
            FieldElementExt::new(m_field, FieldElement::zero())
        } else {
            let num = p1.y.substract(&p2.y);
            let den = p1.x.substract(&p2.x);
            let m_num = num.a;
            let m_den = den.a;
            let m_field = m_num.divide(&m_den);
            FieldElementExt::new(m_field, FieldElement::zero())
        };

        let x3 = m.multiply(&m).substract(&p1.x).substract(&p2.x);
        let y3 = m.multiply(&p1.x.substract(&x3)).substract(&p1.y);
        PointExt {
            x: x3,
            y: y3,
            is_infinity: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    /* These two tests will be later manually calculated as they change with the EC
    #[test]
    fn test_point_addition() {
        let curve = EllipticCurve::new(
            FieldElement::zero(),
            FieldElement::new(7),
            16,
            Point::Affine(FieldElement::zero(), FieldElement::zero()), // doesn't matter the generator for testing
        );
        let p1 = curve.point(FieldElement::new(3), FieldElement::new(0));
        let p2 = curve.point(FieldElement::new(8), FieldElement::new(3));
        let sum = p1 + p2;
        if let Point::Affine(x, y) = sum {
            assert_eq!(x.value, 5);
            assert_eq!(y.value, 9);
        } else {
            panic!("Expected affine point");
        }
    }

    #[test]
    fn test_scalar_mul() {
        let curve = EllipticCurve::new(
            FieldElement::new(0),
            FieldElement::new(7),
            16,
            Point::Affine(FieldElement::zero(), FieldElement::zero()), // doesn't matter the generator for testing
        );
        let p = curve.point(FieldElement::new(8), FieldElement::new(3));
        let result = p.scalar_mul(FieldElement::new(2), &curve);
        if let Point::Affine(x, y) = result {
            assert_eq!(x.value, 5);
            assert_eq!(y.value, 8);
        } else {
            panic!("Expected affine point");
        }
    }
    */

    #[test]
    fn test_curve_order() {
        let curve = EllipticCurve::new();
        let g1 = curve.generator_g1();
        let order = FieldElement::new(curve.order);
        let result = g1.scalar_mul(&curve, order);
        assert_eq!(result, Point::infinity());
    }

    #[test]
    fn test_double_g1() {
        let curve = EllipticCurve::new();
        let g1 = curve.generator_g1(); // (1, 2)
        let double_g1 = curve.add(&g1, &g1);
        assert_eq!(
            double_g1,
            Point {
                x: FieldElement::new(68),
                y: FieldElement::new(74),
                is_infinity: false,
            }
        );
    }
}
