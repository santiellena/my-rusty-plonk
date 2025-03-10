use crate::own_field::FieldElement;

#[derive(Clone, Debug, PartialEq)]
pub struct EllipticCurve {
    // y^2 = x^3 + ax + b
    a: FieldElement,  // Coefficient a
    b: FieldElement,  // Coefficient b
    field_order: u64, // p
    curve_order: u64, // n (number of points)
}

#[derive(Clone, Debug)]
pub enum Point {
    Infinity,                           // Point at infinity (identity)
    Affine(FieldElement, FieldElement), // (x, y)
}

impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Point::Infinity, Point::Infinity) => true,
            (Point::Affine(x1, y1), Point::Affine(x2, y2)) => x1 == x2 && y1 == y2,
            _ => false,
        }
    }
}

impl EllipticCurve {
    #[allow(dead_code)]
    pub fn new(a: FieldElement, b: FieldElement, curve_order: u64) -> Self {
        assert_eq!(a.order, b.order, "Coefficients must be in the same field");
        let field_order: u64 = a.order;
        let sixteen: FieldElement = FieldElement::new(16, field_order);
        let neg_sixteen: FieldElement = sixteen.negate();
        let twenty_seven: FieldElement = FieldElement::new(27, field_order);
        let four: FieldElement = FieldElement::new(4, field_order);

        let discriminant: FieldElement =
            neg_sixteen * ((four * a.pow(3)) + (twenty_seven * b.pow(2)));

        assert!(
            discriminant != FieldElement::zero(field_order),
            "Must be a valid elliptic curve"
        );

        EllipticCurve {
            a,
            b,
            field_order,
            curve_order,
        }
    }

    #[allow(dead_code)]
    pub fn point(&self, x: FieldElement, y: FieldElement) -> Point {
        assert_eq!(
            x.order, self.field_order,
            "Point must be in the curve's field"
        );
        assert_eq!(
            y.order, self.field_order,
            "Point must be in the curve's field"
        );
        // Verify point is on curve: y^2 = x^3 + ax + b
        let lhs = y.pow(2);
        let rhs = (x.pow(3)) + (self.a.multiply(&x)) + self.b.clone();
        assert_eq!(lhs, rhs, "Point ({:?}, {:?}) not on curve", x, y);
        Point::Affine(x, y)
    }

    #[allow(dead_code)]
    pub fn infinity(&self) -> Point {
        Point::Infinity
    }
}

impl std::ops::Add for Point {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        match (self, other) {
            (Point::Infinity, p) => p,
            (p, Point::Infinity) => p,
            (Point::Affine(x1, y1), Point::Affine(x2, y2)) => {
                if x1.order != x2.order {
                    panic!("Points must be in the same field");
                }
                let field_order = x1.order;

                // Check if P = -Q
                if x1 == x2 && y1 == y2.negate() {
                    return Point::Infinity;
                }

                let lambda = if x1 == x2 && y1 == y2 {
                    // Doubling: λ = (3x₁² + a) / (2y₁)
                    let num = (x1.multiply(&x1) * FieldElement::new(3, field_order))
                        + FieldElement::zero(field_order);
                    // Assuming a=0 for simplicity;
                    // Explanation on previous a=0: most used Plonk curves has a=0 and current setup doesn't allow us to access
                    // the Elliptic Curve data of the points. An option here will be creating a custom add function that
                    // receives the elliptic curve in which the points being added belong to. However, I prefer this solution
                    // because I want to use the sintactic sugar this option provides, and for the curves I'll be integrating this
                    // Point structure, it will be more than useful.
                    let denom = y1.clone() * FieldElement::new(2, field_order);
                    num / denom
                } else {
                    // Addition: λ = (y₂ - y₁) / (x₂ - x₁)
                    let num = y2.substract(&y1);
                    let denom = x2.substract(&x1);
                    num / denom
                };

                let x3 = lambda.multiply(&lambda).substract(&x1).substract(&x2);
                let y3 = lambda.multiply(&x1.substract(&x3)).substract(&y1);
                Point::Affine(x3, y3)
            }
        }
    }
}

impl Point {
    #[allow(dead_code)]
    pub fn scalar_mul(&self, scalar: FieldElement, curve: &EllipticCurve) -> Self {
        if let Point::Infinity = self {
            return Point::Infinity;
        }
        let Point::Affine(x, _y) = self else {
            unreachable!()
        };
        assert_eq!(x.order, scalar.order, "Scalar must match point's field");
        let mut result = Point::Infinity;
        let mut base = self.clone();
        let mut exp = scalar.value % curve.curve_order; // Modulo curve order

        while exp > 0 {
            if exp & 1 == 1 {
                result = result + base.clone();
            }
            base = base.clone() + base.clone();
            exp >>= 1;
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_addition() {
        let order = 17;
        let curve = EllipticCurve::new(FieldElement::zero(order), FieldElement::new(7, order), 16);
        let p1 = curve.point(FieldElement::new(3, order), FieldElement::new(0, order));
        let p2 = curve.point(FieldElement::new(8, order), FieldElement::new(3, order));
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
        let order = 17;
        let curve =
            EllipticCurve::new(FieldElement::new(0, order), FieldElement::new(7, order), 16);
        let p = curve.point(FieldElement::new(8, order), FieldElement::new(3, order));
        let result = p.scalar_mul(FieldElement::new(2, order), &curve);
        if let Point::Affine(x, y) = result {
            assert_eq!(x.value, 5);
            assert_eq!(y.value, 8);
        } else {
            panic!("Expected affine point");
        }
    }
}
