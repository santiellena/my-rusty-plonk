use crate::field::FieldElement;
use std::ops::{Add, Mul, Sub};

#[derive(Clone, Debug, PartialEq)]
pub struct Polynomial {
    pub coeffs: Vec<FieldElement>,
}

impl Polynomial {
    pub fn new(coeffs: Vec<FieldElement>) -> Self {
        Polynomial { coeffs }
    }

    #[allow(dead_code)]
    pub fn degree(&self) -> Option<usize> {
        self.coeffs.iter().rposition(|c| *c != FieldElement::zero())
    }

    pub fn add(&self, other: &Self) -> Self {
        let max_len: usize = self.coeffs.len().max(other.coeffs.len());
        let mut result: Vec<FieldElement> = Vec::with_capacity(max_len);
        let zero: FieldElement = FieldElement::zero();
        for i in 0..max_len {
            let a: &FieldElement = self.coeffs.get(i).unwrap_or(&zero);
            let b: &FieldElement = other.coeffs.get(i).unwrap_or(&zero);
            result.push(a.add(b));
        }
        Polynomial::new(result)
    }

    #[allow(dead_code)]
    pub fn scalar_mul(&self, scalar: FieldElement) -> Self {
        let coeffs: Vec<FieldElement> = self.coeffs.iter().map(|c| c.multiply(&scalar)).collect();
        Polynomial::new(coeffs)
    }
}

impl Add for Polynomial {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Polynomial::add(&self, &other)
    }
}

impl Sub for Polynomial {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        let max_len: usize = self.coeffs.len().max(other.coeffs.len());
        let mut result: Vec<FieldElement> = Vec::with_capacity(max_len);
        let zero: FieldElement = FieldElement::zero();
        for i in 0..max_len {
            let a: &FieldElement = self.coeffs.get(i).unwrap_or(&zero);
            let b: &FieldElement = other.coeffs.get(i).unwrap_or(&zero);
            result.push(a.substract(b));
        }
        Polynomial::new(result)
    }
}

impl Mul for Polynomial {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let mut result: Vec<FieldElement> =
            vec![FieldElement::zero(); self.coeffs.len() + other.coeffs.len() - 1];
        for (i, a) in self.coeffs.iter().enumerate() {
            for (j, b) in other.coeffs.iter().enumerate() {
                result[i + j] = result[i + j].add(a.multiply(b));
            }
        }
        Polynomial::new(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        let p1 = Polynomial::new(vec![FieldElement::new(1), FieldElement::new(2)]);
        let p2 = Polynomial::new(vec![FieldElement::new(3), FieldElement::new(4)]);
        let sum = p1 + p2;
        let expected = Polynomial::new(vec![FieldElement::new(4), FieldElement::new(6)]);
        assert_eq!(sum, expected);
    }

    #[test]
    fn test_scalar_mul() {
        let p = Polynomial::new(vec![FieldElement::new(1), FieldElement::new(2)]);
        let scalar = FieldElement::new(3);
        let result = p.scalar_mul(scalar);
        let expected = Polynomial::new(vec![FieldElement::new(3), FieldElement::new(6)]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_polynomial_mul() {
        let p1 = Polynomial::new(vec![FieldElement::new(1), FieldElement::new(2)]);
        let p2 = Polynomial::new(vec![FieldElement::new(3), FieldElement::new(4)]);
        let product = p1 * p2;
        let expected = Polynomial::new(vec![
            FieldElement::new(3),  // 1*3
            FieldElement::new(10), // 1*4 + 2*3
            FieldElement::new(8),  // 2*4
        ]);
        assert_eq!(product, expected);
    }

    #[test]
    fn test_degree() {
        let p = Polynomial::new(vec![
            FieldElement::new(0),
            FieldElement::new(2),
            FieldElement::new(3),
        ]);
        assert_eq!(p.degree(), Some(2));
        let p_zero = Polynomial::new(vec![FieldElement::new(0)]);
        assert_eq!(p_zero.degree(), None);
    }
}
