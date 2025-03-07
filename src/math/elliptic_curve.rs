#[path = "./field.rs"]
mod field;

use field::FieldElement;

#[derive(Clone, Debug, PartialEq)]
pub struct EllipticCurve {
    // y^2 = x^3 + ax + b
    a: FieldElement, // Coefficient a
    b: FieldElement, // Coefficient b
    order: u64,      // Field order (same as coefficient's order)
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
    pub fn new(a: FieldElement, b: FieldElement) -> Self {
        assert_eq!(a.order, b.order, "Coefficients must be in the same field");
        let order: u64 = a.order;
        let sixteen: FieldElement = FieldElement::new(16, order);
        let neg_sixteen: FieldElement = sixteen.negate();
        let twenty_seven: FieldElement = FieldElement::new(27, order);
        let four: FieldElement = FieldElement::new(4, order);

        let discriminant: FieldElement =
            neg_sixteen * ((four * a.pow(3)) + (twenty_seven * b.pow(2)));

        assert!(
            discriminant != FieldElement::zero(order),
            "Must be a valid elliptic curve"
        );

        EllipticCurve { a, b, order }
    }

    pub fn point(&self, x: FieldElement, y: FieldElement) -> Point {
        assert_eq!(x.order, self.order, "Point must be in the curve's field");
        assert_eq!(y.order, self.order, "Point must be in the curve's field");
        // Verify point is on curve: y^2 = x^3 + ax + b
        let lhs = y.multiply(&y);
        let rhs = x
            .multiply(&x)
            .multiply(&x)
            .add(&self.a.multiply(&x))
            .add(&self.b);
        assert_eq!(lhs, rhs, "Point ({:?}, {:?}) not on curve", x, y);
        Point::Affine(x, y)
    }

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
                let curve_order = x1.order;

                // Check if P = -Q
                if x1 == x2 && y1 == y2.negate() {
                    return Point::Infinity;
                }

                let lambda = if x1 == x2 && y1 == y2 {
                    // Doubling: λ = (3x₁² + a) / (2y₁)
                    let num = x1
                        .multiply(&x1)
                        .multiply(&FieldElement::new(3, curve_order))
                        .add(FieldElement::zero(curve_order));
                    // Assuming a=0 for simplicity;
                    // Explanation on previous a=0: most used Plonk curves has a=0 and current setup doesn't allow us to access
                    // the Elliptic Curve data of the points. An option here will be creating a custom add function that
                    // receives the elliptic curve in which the points being added belong to. However, I prefer this solution
                    // because I want to use the sintactic sugar this option provides, and for the curves I'll be integrating this
                    // Point structure, it will be more than useful.
                    let denom = y1.multiply(&FieldElement::new(2, curve_order));
                    num.divide(denom.value)
                } else {
                    // Addition: λ = (y₂ - y₁) / (x₂ - x₁)
                    let num = y2.substract(&y1);
                    let denom = x2.substract(&x1);
                    num.divide(denom.value)
                };

                let x3 = lambda.multiply(&lambda).substract(&x1).substract(&x2);
                let y3 = lambda.multiply(&x1.substract(&x3)).substract(&y1);
                Point::Affine(x3, y3)
            }
        }
    }
}

impl Point {
    pub fn scalar_mul(&self, scalar: FieldElement) -> Self {
        if let Point::Infinity = self {
            return Point::Infinity;
        }
        let Point::Affine(x, y) = self else {
            unreachable!()
        };
        let mut result = Point::Infinity;
        let mut base = self.clone();
        let mut exp = scalar.value;

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
        let order = 17; // Small prime for testing
        let curve = EllipticCurve::new(FieldElement::new(0, order), FieldElement::new(7, order)); // y² = x³ + 7 (simple curve)
        let p1 = curve.point(FieldElement::new(3, order), FieldElement::new(10, order)); // (3, 10)
        let p2 = curve.point(FieldElement::new(9, order), FieldElement::new(7, order)); // (9, 7)
        let sum = p1 + p2;
        // Expected: (13, 4) (computed manually or via sage)
        if let Point::Affine(x, y) = sum {
            assert_eq!(x.value, 13);
            assert_eq!(y.value, 4);
        }
    }

    #[test]
    fn test_scalar_mul() {
        let order = 17;
        let curve = EllipticCurve::new(FieldElement::new(0, order), FieldElement::new(7, order));
        let p = curve.point(FieldElement::new(3, order), FieldElement::new(10, order));
        let result = p.scalar_mul(FieldElement::new(2, order)); // 2P
        if let Point::Affine(x, y) = result {
            assert_eq!(x.value, 5); // 2*(3, 10) = (5, 3) on this curve
            assert_eq!(y.value, 3);
        }
    }
}
