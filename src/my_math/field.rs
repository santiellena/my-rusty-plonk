use crate::gcd;

use std::ops::{Add, AddAssign, Div, Mul};

#[derive(Clone, Debug)]
pub struct FieldElement {
    pub value: u64,
    pub order: u64,
}

impl FieldElement {
    /// Constructor
    pub fn new(value: u64, order: u64) -> Self {
        Self {
            value: value % order,
            order,
        }
    }

    /// Addition in the FieldElement
    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(self.order, other.order, "Fields MUST have the same order");
        Self::new(self.value + other.value, self.order)
    }

    /// Subtraction in the FieldElement
    #[allow(dead_code)]
    pub fn substract(&self, other: &Self) -> Self {
        assert_eq!(self.order, other.order, "Fields MUST have the same order");
        Self::new(self.value + self.order - other.value, self.order)
    }

    /// Multiplication in the FieldElement
    pub fn multiply(&self, other: &Self) -> Self {
        assert_eq!(self.order, other.order, "Fields MUST have the same order");
        Self::new(self.value * other.value, self.order)
    }

    /// Modular inverse using Extended Euclidean Algorithm
    #[allow(dead_code)]
    fn mod_inverse(&self) -> Option<u64> {
        let (s, _, gcd) = gcd::ext_gcd(self.value as i64, self.order as i64);
        if gcd != 1 {
            return None;
        }
        Some(((s % self.order as i64 + self.order as i64) % self.order as i64) as u64)
    }

    /// Modular inverse (for division)
    #[allow(dead_code)]
    pub fn inverse(&self) -> Self {
        if let Some(inv) = self.mod_inverse() {
            Self {
                value: inv,
                order: self.order,
            }
        } else {
            Self::zero(self.order)
        }
    }

    /// Division (a / b = a * b⁻¹ mod p)
    #[allow(dead_code)]
    pub fn divide(&self, other: &Self) -> Self {
        let inverse_other: &FieldElement = &FieldElement {
            value: other.inverse().value,
            order: other.order,
        };

        self.multiply(inverse_other)
    }

    /// Modular exponentiation (a^exp mod p) using fast exponentiation
    #[allow(dead_code)]
    pub fn pow(&self, exp: u64) -> Self {
        if exp == 0 {
            Self::one(self.order)
        } else {
            let mut fsb_found: bool = false; // first significant byte found
            let mut result: Self = self.clone();

            for n in (0..64).rev() {
                let bin: u64 = (exp >> n) & 1;
                if fsb_found {
                    if bin == 1 {
                        result = result.multiply(&result).multiply(&self);
                    } else {
                        result = result.multiply(&result);
                    }
                } else if bin == 1 {
                    fsb_found = true;
                } else {
                    continue;
                }
            }
            result
        }
    }

    /// Negation in the FieldElement
    #[allow(dead_code)]
    pub fn negate(&self) -> Self {
        Self::new(self.order, self.order).substract(&self)
    }

    /// Zero element
    pub fn zero(order: u64) -> Self {
        Self::new(0, order)
    }

    /// One element
    pub fn one(order: u64) -> Self {
        Self::new(1, order)
    }
}

/// Implementations to facilitate writing the code in polynomials
impl Add for FieldElement {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self::add(&self, &other)
    }
}

impl Mul for FieldElement {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        self.multiply(&other)
    }
}

impl Default for FieldElement {
    fn default() -> Self {
        Self::zero(0)
    }
}

impl AddAssign for FieldElement {
    fn add_assign(&mut self, other: Self) {
        *self = Self::add(&self, &other);
    }
}

// use PartialEq trait to check whether two fields are equal or not
impl PartialEq for FieldElement {
    fn eq(&self, other: &FieldElement) -> bool {
        self.value == other.value && self.order == other.order
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
    fn test_negate() {
        let x: FieldElement = FieldElement {
            value: 3,
            order: 14,
        };
        let y: FieldElement = x.negate() + x;
        println!("{:?}", y);
        assert_eq!(y.value, 0);
    }

    #[test]
    fn test_inverse_basic() {
        // Field F_7 (order = 7, a prime)
        let order = 7;

        // Test 1: Inverse of 1 should be 1 (1 * 1 ≡ 1 mod 7)
        let a = FieldElement::new(1, order);
        let inv_a = a.inverse();
        println!("Inverse: {:?}", inv_a);
        assert_eq!(inv_a.value, 1);
        assert_eq!(a.multiply(&inv_a).value, 1);

        // Test 2: Inverse of 2 (2 * 4 ≡ 1 mod 7)
        let b = FieldElement::new(2, order);
        let inv_b = b.inverse();
        assert_eq!(inv_b.value, 4);
        assert_eq!(b.multiply(&inv_b).value, 1);

        // Test 3: Inverse of 3 (3 * 5 ≡ 1 mod 7)
        let c = FieldElement::new(3, order);
        let inv_c = c.inverse();
        assert_eq!(inv_c.value, 5);
        assert_eq!(c.multiply(&inv_c).value, 1);
    }

    #[test]
    fn test_inverse_non_invertible() {
        // Field F_7
        let order = 7;

        // Test: 0 has no inverse, should return 0
        let zero = FieldElement::new(0, order);
        let inv_zero = zero.inverse();
        assert_eq!(inv_zero.value, 0); // By your design, returns zero when no inverse exists

        // Field F_15 (not a prime, so not all elements are invertible)
        let order = 15;
        let d = FieldElement::new(5, order); // gcd(5, 15) = 5 ≠ 1
        let inv_d = d.inverse();
        assert_eq!(inv_d.value, 0); // No inverse exists
    }

    #[test]
    fn test_pow_basic() {
        let order = 7;

        // Test 1: a^0 = 1
        let a = FieldElement::new(3, order);
        let pow_0 = a.pow(0);
        assert_eq!(pow_0.value, 1);

        // Test 2: a^1 = a
        let pow_1 = a.pow(1);
        assert_eq!(pow_1.value, 3);

        // Test 3: a^2 = a * a (3 * 3 = 9 ≡ 2 mod 7)
        let pow_2 = a.pow(2);
        assert_eq!(pow_2.value, 2);
        assert_eq!(pow_2, a.multiply(&a));

        // Test 4: a^3 = a * a * a
        let pow_3 = a.pow(3);
        assert_eq!(pow_3.value, 6);
        assert_eq!(pow_3, a.clone() * a.clone() * a.clone());
    }

    #[test]
    fn test_pow_larger_exponent() {
        let order = 11;

        // Test: 2^5 mod 11
        // 2^1 = 2, 2^2 = 4, 2^3 = 8, 2^4 = 16 ≡ 5, 2^5 = 32 ≡ 10 mod 11
        let b = FieldElement::new(2, order);
        let pow_5 = b.pow(5);
        assert_eq!(pow_5.value, 10);

        // Manual verification
        let b2 = b.multiply(&b); // 4
        let b4 = b2.multiply(&b2); // 16 ≡ 5
        let b5 = b4.multiply(&b); // 5 * 2 = 10
        assert_eq!(b5.value, 10);
    }

    #[test]
    fn test_pow_zero_base() {
        let order = 7;

        // Test: 0^0 = 1 (by convention in your implementation)
        let zero = FieldElement::new(0, order);
        let pow_0 = zero.pow(0);
        assert_eq!(pow_0.value, 1);

        // Test: 0^1 = 0
        let pow_1 = zero.pow(1);
        assert_eq!(pow_1.value, 0);

        // Test: 0^2 = 0
        let pow_2 = zero.pow(2);
        assert_eq!(pow_2.value, 0);
    }
}
