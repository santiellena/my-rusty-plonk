use ark_bn254::Fr;
use ark_ff::{Field, PrimeField};
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Clone, PartialEq, Copy)]
pub struct FieldElement(Fr);

impl FieldElement {
    #[allow(dead_code)]
    pub fn new(value: u64) -> Self {
        FieldElement(Fr::from(value))
    }

    #[allow(dead_code)]
    pub fn zero() -> Self {
        FieldElement(Fr::from(0u64))
    }

    #[allow(dead_code)]
    pub fn one() -> Self {
        FieldElement(Fr::from(1u64))
    }

    #[allow(dead_code)]
    pub fn multiply(&self, other: &Self) -> Self {
        FieldElement(self.0 * other.0)
    }

    pub fn add(&self, other: &Self) -> Self {
        FieldElement(self.0 + other.0)
    }

    pub fn substract(&self, other: &Self) -> Self {
        FieldElement(self.0 - other.0)
    }

    pub fn divide(&self, other: &Self) -> Self {
        FieldElement(self.0 / other.0)
    }

    #[allow(dead_code)]
    pub fn negate(&self) -> Self {
        FieldElement(-self.0)
    }

    #[allow(dead_code)]
    pub fn inverse(&self) -> Self {
        FieldElement(self.0.inverse().expect("Inverse exists"))
    }

    #[allow(dead_code)]
    pub fn pow(&self, exp: u64) -> Self {
        FieldElement(self.0.pow(&[exp]))
    }
}

impl fmt::Debug for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let bigint = self.0.into_bigint();
        write!(f, "FieldElement({})", bigint)
    }
}

impl Add for FieldElement {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        FieldElement::add(&self, &other)
    }
}

impl Sub for FieldElement {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self.substract(&other)
    }
}

impl Mul for FieldElement {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        self.multiply(&other)
    }
}

impl Div for FieldElement {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        self.divide(&other)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_field_element_add() {
        let a = FieldElement::new(5);
        let b = FieldElement::new(10);
        let c = a + b;
        assert_eq!(c.0, Fr::from(15));
    }

    #[test]
    fn test_field_element_mul() {
        let a = FieldElement::new(3);
        let b = FieldElement::new(4);
        let c = a * b;
        assert_eq!(c.0, Fr::from(12));
    }

    #[test]
    fn test_field_element_pow() {
        let a = FieldElement::new(3);
        let c = a.pow(3);
        assert_eq!(c.0, Fr::from(27));
    }

    #[test]
    fn test_field_element_inverse() {
        let a = FieldElement::new(5);
        let inv = a.inverse();
        assert_eq!(a.multiply(&inv), FieldElement::one());
    }
}
