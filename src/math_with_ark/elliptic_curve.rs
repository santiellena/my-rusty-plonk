use crate::field::FieldElement;
use ark_bn254::g1::G1Affine;
// use ark_bn254::Fr;
use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::Zero;
use std::ops::{Add, Mul};

#[derive(Clone, Debug, PartialEq)]
pub struct EllipticCurve {
    generator: G1Affine,
    // curve_order: FieldElement, // G1 subgroup order won't be keeping this value
    /* Why not keeping curve_order here?

       Well, FieldElement is just a wrapper of ark-fr::Fr type. This type accepts
       values from 0-Fr::MODULUS, non-inclusive. Storing the curve order could be
       helpful for testing properties like the curve order times (scalar multiplication)
       the generator is the point at infinity (identity).
                               [r]G = O
    */
}

#[derive(Clone, Debug, PartialEq)]
pub enum Point {
    Infinity,
    #[allow(dead_code)]
    Affine(G1Affine),
}

impl EllipticCurve {
    #[allow(dead_code)]
    pub fn new() -> Self {
        let generator = G1Affine::generator();
        println!("Generator: {:?}", generator);

        /*let mut modulus = Fr::MODULUS; // Copy BigInteger256
        modulus.0[0] -= 1; // Subtract 1 from the least significant word
        let curve_order = FieldElement::new_bigint(modulus);
        println!("Curve order: {:?}", curve_order);
        */
        EllipticCurve {
            generator,
            //curve_order,
        }
    }

    #[allow(dead_code)]
    pub fn point(&self, scalar: FieldElement) -> Point {
        let result = self.generator.mul(scalar.value());
        if result.is_zero() {
            Point::Infinity
        } else {
            Point::Affine(result.into_affine())
        }
    }

    #[allow(dead_code)]
    pub fn infinity(&self) -> Point {
        Point::Infinity
    }

    #[allow(dead_code)]
    pub fn generator(&self) -> Point {
        Point::Affine(self.generator)
    }
}

impl Add for Point {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        match (self, other) {
            (Point::Infinity, p) => p,
            (p, Point::Infinity) => p,
            (Point::Affine(p1), Point::Affine(p2)) => {
                let sum = p1 + p2; // Returns Projective
                if sum.is_zero() {
                    Point::Infinity
                } else {
                    Point::Affine(sum.into_affine())
                }
            }
        }
    }
}

impl Point {
    #[allow(dead_code)]
    pub fn scalar_mul(&self, scalar: FieldElement, _curve: &EllipticCurve) -> Self {
        match self {
            Point::Infinity => Point::Infinity,
            Point::Affine(p) => {
                let result = p.mul(scalar.value());
                if result.is_zero() {
                    Point::Infinity
                } else {
                    Point::Affine(result.into_affine())
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_addition() {
        let curve = EllipticCurve::new();
        let p1 = curve.generator();
        let p2 = curve.point(FieldElement::new(2));
        let sum = p1 + p2;
        let expected = curve.point(FieldElement::new(3));
        assert_eq!(sum, expected);
    }

    #[test]
    fn test_scalar_mul() {
        let curve = EllipticCurve::new();
        let p = curve.generator();
        let result = p.scalar_mul(FieldElement::new(3), &curve);
        let expected = curve.point(FieldElement::new(3));
        assert_eq!(result, expected);
    }

    #[test]
    fn test_infinity() {
        let curve = EllipticCurve::new();
        let p = curve.generator();
        let inf = curve.infinity();
        let sum = p.clone() + inf;
        assert_eq!(sum, p);
    }

    #[test]
    fn test_basic_mul() {
        let curve = EllipticCurve::new();
        let g = curve.generator();
        let result = g.scalar_mul(FieldElement::new(1), &curve);
        assert_eq!(result, g);
    }

    #[test]
    fn test_new() {
        let _curve = EllipticCurve::new();
        assert!(true);
    }
}
