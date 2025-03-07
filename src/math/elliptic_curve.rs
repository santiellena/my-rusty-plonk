use std::fmt;

#[derive(Clone, Debug, PartialEq)]
pub struct EllipticCurve {
    // y^2 = x^3 + ax + b
    a: FieldElement, // Coefficient a
    b: FieldElement, // Coefficient b
    order: u64,      // Field order (same as FieldElement's order)
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
        // Could check discriminant here, but we'll assume valid a, b for now
        EllipticCurve {
            a,
            b,
            order: a.order,
        }
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
