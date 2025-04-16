use crate::gcd;

use std::ops::{AddAssign, Div, Mul};

#[derive(Clone, Debug)]
pub struct FieldElement {
    pub value: u64,
}

impl FieldElement {
    pub const MODULUS: u64 = 101;

    /// Constructor
    pub fn new(value: u64) -> Self {
        Self {
            value: value % Self::MODULUS,
        }
    }

    /// Addition in the FieldElement
    pub fn add(&self, other: &Self) -> Self {
        Self::new(self.value + other.value)
    }

    /// Subtraction in the FieldElement
    #[allow(dead_code)]
    pub fn substract(&self, other: &Self) -> Self {
        Self::new(self.value + Self::MODULUS - other.value)
    }

    /// Multiplication in the FieldElement
    pub fn multiply(&self, other: &Self) -> Self {
        Self::new(self.value * other.value)
    }

    /// Modular inverse using Extended Euclidean Algorithm
    #[allow(dead_code)]
    fn mod_inverse(&self) -> Option<u64> {
        let (s, _, gcd) = gcd::ext_gcd(self.value as i64, Self::MODULUS as i64);
        if gcd != 1 {
            return None;
        }
        Some(((s % Self::MODULUS as i64 + Self::MODULUS as i64) % Self::MODULUS as i64) as u64)
    }

    /// Modular inverse (for division)
    #[allow(dead_code)]
    pub fn inverse(&self) -> Self {
        if let Some(inv) = self.mod_inverse() {
            Self { value: inv }
        } else {
            Self::zero()
        }
    }

    /// Division (a / b = a * b⁻¹ mod p)
    #[allow(dead_code)]
    pub fn divide(&self, other: &Self) -> Self {
        let inverse_other: &FieldElement = &FieldElement {
            value: other.inverse().value,
        };

        self.multiply(inverse_other)
    }

    /// Modular exponentiation (a^exp mod p) using fast exponentiation
    #[allow(dead_code)]
    pub fn pow(&self, exp: u64) -> Self {
        if exp == 0 {
            Self::one()
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
        Self::new(Self::MODULUS).substract(&self)
    }

    /// Zero element
    pub fn zero() -> Self {
        Self::new(0)
    }

    /// One element
    pub fn one() -> Self {
        Self::new(1)
    }
}

/// Implementations to facilitate writing the code in polynomials

impl Mul for FieldElement {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        self.multiply(&other)
    }
}

impl Default for FieldElement {
    fn default() -> Self {
        Self::zero()
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
        self.value == other.value
    }
}

impl Div for FieldElement {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        self.divide(&other)
    }
}

/*
    Field element extension F_(p^k)
*/
#[derive(Clone, Debug, PartialEq)]
pub struct FieldElementExt {
    a: FieldElement, // Real part
    b: FieldElement, // Imaginary part (coefficient of u)
}

impl FieldElementExt {
    #[allow(dead_code)]
    pub fn new(a: FieldElement, b: FieldElement) -> Self {
        FieldElementExt { a, b }
    }

    #[allow(dead_code)]
    pub fn add(&self, other: &Self) -> Self {
        FieldElementExt {
            a: self.a.add(&other.a),
            b: self.b.add(&other.b),
        }
    }

    #[allow(dead_code)]
    pub fn sub(&self, other: &Self) -> Self {
        FieldElementExt {
            a: self.a.substract(&other.a),
            b: self.b.substract(&other.b),
        }
    }

    #[allow(dead_code)]
    pub fn mul(&self, other: &Self) -> Self {
        // (a + bu)(c + du) = (ac + bd u^2) + (ad + bc)u
        // u^2 = -2, so bd u^2 = bd(-2)
        let ac = self.a.multiply(&other.a);
        let ad_u = self.a.multiply(&other.b);
        let bc_u = self.b.multiply(&other.a);
        let bd_u2 = self
            .b
            .multiply(&other.b)
            .multiply(&FieldElement::new(101 - 2)); // u^2 = -2

        let real = ac.add(&bd_u2);
        let imag = ad_u.add(&bc_u);

        FieldElementExt { a: real, b: imag }
    }

    #[allow(dead_code)]
    pub fn scalar_mul(&self, scalar: FieldElement) -> Self {
        FieldElementExt {
            a: self.a.multiply(&scalar),
            b: self.b.multiply(&scalar),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_negate() {
        let x: FieldElement = FieldElement { value: 3 };
        let y: FieldElement = x.negate().add(&x);
        println!("{:?}", y);
        assert_eq!(y.value, 0);
    }

    #[test]
    fn test_inverse_basic() {
        // Test 1: Inverse of 1 should be 1 (1 * 1 ≡ 1 mod 101)
        let a = FieldElement::new(1);
        let inv_a = a.inverse();
        println!("Inverse: {:?}", inv_a);
        assert_eq!(inv_a.value, 1);
        assert_eq!(a.multiply(&inv_a).value, 1);

        // Test 2: Inverse of 2 (2 * 51 ≡ 1 mod 101)
        let b = FieldElement::new(2);
        let inv_b = b.inverse();
        assert_eq!(inv_b.value, 51);
        assert_eq!(b.multiply(&inv_b).value, 1);

        // Test 3: Inverse of 3 (3 * 34 ≡ 1 mod 101)
        let c = FieldElement::new(3);
        let inv_c = c.inverse();
        assert_eq!(inv_c.value, 34);
        assert_eq!(c.multiply(&inv_c).value, 1);
    }

    // #[test]
    // Irreproducible test with fixed modulus changes but leaving it here because it is helpful
    /*fn test_inverse_non_invertible() {
        // Field F_7
        let order = 7;

        // Test: 0 has no inverse, should return 0
        let zero = FieldElement::new(0);
        let inv_zero = zero.inverse();
        assert_eq!(inv_zero.value, 0); // By your design, returns zero when no inverse exists

        // Field F_15 (not a prime, so not all elements are invertible)
        let order = 15;
        let d = FieldElement::new(5, order); // gcd(5, 15) = 5 ≠ 1
        let inv_d = d.inverse();
        assert_eq!(inv_d.value, 0); // No inverse exists
    }*/

    #[test]
    fn test_pow_basic() {
        // Test 1: a^0 = 1
        let a = FieldElement::new(3);
        let pow_0 = a.pow(0);
        assert_eq!(pow_0.value, 1);

        // Test 2: a^1 = a
        let pow_1 = a.pow(1);
        assert_eq!(pow_1.value, 3);

        // Test 3: a^2 = a * a (3 * 3 = 9 ≡ 9 mod 101)
        let pow_2 = a.pow(2);
        assert_eq!(pow_2.value, 9);
        assert_eq!(pow_2, a.multiply(&a));

        // Test 4: a^3 = a * a * a
        let pow_3 = a.pow(3);
        assert_eq!(pow_3.value, 27);
        assert_eq!(pow_3, a.clone() * a.clone() * a.clone());
    }

    #[test]
    fn test_pow_larger_exponent() {
        // Test: 2^5 mod 101
        // 2^1 = 2, 2^2 = 4, 2^3 = 8, 2^4 = 16 ≡ 5, 2^8 = 256 ≡ 54 mod 101
        let b = FieldElement::new(2);
        let pow_8 = b.pow(8);
        assert_eq!(pow_8.value, 54);

        // Manual verification
        let b2 = b.multiply(&b); // 4
        let b4 = b2.multiply(&b2); // 16
        let b8 = b4.multiply(&b4); // 54 == 256
        assert_eq!(b8.value, 54);
    }

    #[test]
    fn test_pow_zero_base() {
        // Test: 0^0 = 1 (by convention in your implementation)
        let zero = FieldElement::new(0);
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

#[cfg(test)]
mod tests_ext {
    use super::*;

    #[test]
    fn test_field_ext_operations() {
        let a = FieldElementExt::new(FieldElement::new(50), FieldElement::new(60));
        let b = FieldElementExt::new(FieldElement::new(30), FieldElement::new(40));
        let sum = a.add(&b);
        assert_eq!(sum.a.value, 80); // 50 + 30
        assert_eq!(sum.b.value, 100); // 60 + 40
        let prod = a.mul(&b);
        // (50 + 60u)(30 + 40u) = 1500 + 2000u + 1800u + 2400u^2
        // = (86 + 77(99) + (81 + 83)u (mod 101)
        // = (86 + 48) + 63u = 33 + 63u (mod 101)
        assert_eq!(prod.a.value, 33);
        assert_eq!(prod.b.value, 63);
    }
}
