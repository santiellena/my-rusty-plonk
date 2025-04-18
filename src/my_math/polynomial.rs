use crate::my_field::FieldElement;

use std::ops::{Add, Mul, Sub};

/// Note: A polynomial built with the following coefficients:
///     [1, 2, 3, 4] will represent y = 1 + 2x + 3x^2 + 4x^3
///     I want to make this clear because we tend to read polynomials
///     from the highest degree to the lowest degree, but this implementation
///     does as shown avobe. Thanks for reading!

#[derive(Debug, Clone)]
pub struct Polynomial {
    pub coeffs: Vec<FieldElement>,
}

impl Polynomial {
    pub fn new(coeffs: Vec<FieldElement>) -> Self {
        let mut p: Polynomial = Polynomial {
            coeffs: coeffs
                .into_iter()
                .map(|c| FieldElement::new(c.value))
                .collect(),
        };
        p.trim();
        p
    }

    fn trim(&mut self) {
        while self.coeffs.len() > 1 && self.coeffs.last().unwrap() == &FieldElement::zero() {
            self.coeffs.pop();
        }
    }

    #[allow(dead_code)]
    pub fn evaluate(&self, x: FieldElement) -> FieldElement {
        self.coeffs
            .iter()
            .rev()
            .fold(FieldElement::zero(), |acc, coeff| {
                acc.multiply(&x).add(&coeff)
            })
    }

    #[allow(dead_code)]
    fn scalar_mul(&self, scalar: FieldElement) -> Self {
        Polynomial {
            coeffs: self.coeffs.iter().map(|c| c.multiply(&scalar)).collect(),
        }
    }

    #[allow(dead_code)]
    pub fn divide(&self, divisor: &Self) -> (Self, Self) {
        assert!(
            !divisor.coeffs.is_empty(),
            "Cannot divide by zero polynomial"
        );

        // If divisor is a constant (degree 0), handle scalar division
        if divisor.coeffs.len() == 1 {
            let scalar: &FieldElement = &divisor.coeffs[0];
            let inv_scalar: FieldElement = scalar.inverse();
            let quotient_coeffs: Vec<FieldElement> = self
                .coeffs
                .iter()
                .map(|c| c.multiply(&inv_scalar))
                .collect();
            return (
                Polynomial::new(quotient_coeffs),
                Polynomial::new(vec![]), // Remainder is 0
            );
        }

        let mut dividend: Vec<FieldElement> = self.coeffs.clone();
        let divisor_deg: usize = divisor.coeffs.len() - 1;
        let divisor_lead: &FieldElement = divisor.coeffs.last().unwrap(); // Leading coefficient
        let inv_divisor_lead: FieldElement = divisor_lead.inverse();

        // If dividend degree < divisor degree, quotient is 0, remainder is dividend
        if dividend.len() <= divisor_deg {
            return (Polynomial::new(vec![]), Polynomial::new(dividend));
        }

        let mut quotient: Vec<FieldElement> =
            vec![FieldElement::zero(); dividend.len() - divisor_deg];
        for i in (0..quotient.len()).rev() {
            let dividend_deg = i + divisor_deg;
            let term = dividend[dividend_deg].clone().multiply(&inv_divisor_lead);
            quotient[i] = term.clone();

            // Subtract term * divisor from dividend
            for j in 0..divisor.coeffs.len() {
                let prod: FieldElement = divisor.coeffs[j].clone().multiply(&term);
                dividend[dividend_deg - (divisor_deg - j)] = dividend
                    [dividend_deg - (divisor_deg - j)]
                    .clone()
                    .substract(&prod);
            }
        }

        // Trim leading zeros from remainder
        while dividend.len() > divisor_deg && dividend.last().unwrap() == &FieldElement::zero() {
            dividend.pop();
        }

        (Polynomial::new(quotient), Polynomial::new(dividend))
    }

    #[allow(dead_code)]
    pub fn degree(&self) -> usize {
        if self.coeffs.is_empty() || self.coeffs.iter().all(|c| c == &FieldElement::zero()) {
            return 0; // Convention: zero polynomial has degree 0 (or -∞ in some contexts)
        }
        self.coeffs.len() - 1
    }

    /// This vanishing polynomial calculation instead of the classic: ∏(x -y_i)
    /// is possible because we are working in a cyclic subgroup.
    /// This way is more efficient.
    #[allow(dead_code)]
    pub fn vanishing_polynomial(n: usize) -> Self {
        let mut coeffs = vec![FieldElement::zero(); n + 1];
        coeffs[0] = FieldElement::new(1).negate(); // -1
        coeffs[n] = FieldElement::one(); // 1
        Polynomial::new(coeffs)
    }

    /// Resource I recommend to understand Lagrange Interpolation:
    /// LambdaClass YT video: https://www.youtube.com/watch?v=REnFOKo9gXs
    ///
    /// P(x) = ∑_i y_i ⋅ l_i(x)
    /// where l_i(x)= ∏_j≠i (x−xj)/(xi−xj)
    #[allow(dead_code)]
    pub fn lagrange_interpolate(points: &[(FieldElement, FieldElement)]) -> Self {
        assert!(!points.is_empty(), "Need at least one point");
        let mut result = Polynomial::new(vec![]);

        for (i, &(ref xi, ref yi)) in points.iter().enumerate() {
            let mut term = Polynomial::new(vec![FieldElement::one()]); // Start with 1
            let mut denominator = FieldElement::one();

            for (j, &(ref xj, _)) in points.iter().enumerate() {
                if i != j {
                    // Numerator: (x - xj)
                    let num = Polynomial::new(vec![xj.clone().negate(), FieldElement::one()]);
                    term = term * num;
                    // Denominator: (xi - xj)
                    denominator = denominator.multiply(&xi.substract(xj));
                }
            }

            let inv_denominator = denominator.inverse();
            term = term.scalar_mul(inv_denominator);
            term = term.scalar_mul(yi.clone());
            result = result + term;
        }

        result
    }
}

impl Add for Polynomial {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let max_len: usize = self.coeffs.len().max(other.coeffs.len());
        let mut result: Vec<FieldElement> = Vec::with_capacity(max_len);

        for i in 0..max_len {
            let a: &FieldElement = if i < self.coeffs.len() {
                &self.coeffs[i]
            } else {
                &FieldElement::zero()
            };
            let b: &FieldElement = if i < other.coeffs.len() {
                &other.coeffs[i]
            } else {
                &FieldElement::zero()
            };
            result.push(a.add(b));
        }

        Polynomial::new(result)
    }
}

impl Mul for Polynomial {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let n = self.coeffs.len();
        let m = other.coeffs.len();
        let mut result: Vec<FieldElement> = vec![FieldElement::zero(); n + m - 1];

        for i in 0..n {
            for j in 0..m {
                let prod: FieldElement = self.coeffs[i].clone().multiply(&other.coeffs[j]);
                result[i + j] = result[i + j].clone().add(&prod);
            }
        }

        Polynomial::new(result)
    }
}

impl Sub for Polynomial {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        let negated: Polynomial = Polynomial {
            coeffs: other.coeffs.into_iter().map(|c| c.negate()).collect(),
        };
        self + negated
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polys() {
        // P(x) = 1 + 2x + 3x^2 (mod 101)
        let p1 = Polynomial::new(vec![
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
        ]);
        // Q(x) = 5 + 4x (mod 101)
        let p2 = Polynomial::new(vec![FieldElement::new(5), FieldElement::new(4)]);

        // Addition: (1 + 2x + 3x^2) + (5 + 4x) = 6 + 6x + 3x^2 (mod 101)
        let sum = p1.clone() + p2.clone();
        println!("Sum: {:?}", sum.coeffs); // Should be [6, 6, 3]
        assert_eq!(
            sum.coeffs,
            vec![
                FieldElement::new(6),
                FieldElement::new(6),
                FieldElement::new(3)
            ]
        );

        // Multiplication: (1 + 2x + 3x^2)(5 + 4x) = 5 + 14x + 23x^2 + 12x^3 (mod 101)
        let product = p1.clone() * p2.clone();
        println!("Product: {:?}", product.coeffs); // Should be [5, 14, 23, 12]
        assert_eq!(
            product.coeffs,
            vec![
                FieldElement::new(5),
                FieldElement::new(14),
                FieldElement::new(23),
                FieldElement::new(12)
            ]
        );

        // Evaluate P1 at x = 2: 1 + 2*2 + 3*4 = 1 + 4 + 12 = 17 ≡ 17 (mod 101)
        let value = p1.evaluate(FieldElement::new(2));
        println!("P1(2) = {:?}", value); // Should be FieldElement { value: 17 }
        assert_eq!(value, FieldElement::new(17));
    }

    #[test]
    fn test_polynomial_add() {
        let p1 = Polynomial::new(vec![FieldElement::new(1), FieldElement::new(2)]);
        let p2 = Polynomial::new(vec![FieldElement::new(3), FieldElement::new(4)]);
        let sum = p1 + p2;
        assert_eq!(sum.coeffs, vec![FieldElement::new(4), FieldElement::new(6)]);
    }

    #[test]
    fn test_polynomial_divide() {
        // P(x) = 1 + 2x + 3x^2
        let p = Polynomial::new(vec![
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
        ]);
        // D(x) = 2 + x
        let d = Polynomial::new(vec![FieldElement::new(2), FieldElement::new(1)]);
        let (q, r) = p.divide(&d);
        // Expected: Q(x) = 3x - 4 (3 mod 101), R(x) = 2 (since 3x^2 + 2x + 1 = (x + 2)(3x - 4) + 2)

        // WILL PROBABLY FAIL
        assert_eq!(q.coeffs, vec![FieldElement::new(97), FieldElement::new(3)]); // -4 ≡ 97 mod 101
        assert_eq!(r.coeffs, vec![FieldElement::new(9)]);
    }

    #[test]
    fn test_polynomial_degree() {
        let p = Polynomial::new(vec![
            FieldElement::new(1),
            FieldElement::new(2),
            FieldElement::new(3),
        ]);
        assert_eq!(p.degree(), 2);
        let zero = Polynomial::new(vec![]);
        assert_eq!(zero.degree(), 0);
    }

    #[test]
    fn test_vanishing_polynomial() {
        let n = 2; // x^2 - 1
        let z = Polynomial::vanishing_polynomial(n);
        assert_eq!(
            z.coeffs,
            vec![
                FieldElement::new(100), // -1 ≡ 100 mod 101
                FieldElement::new(0),
                FieldElement::new(1),
            ]
        );
        // Roots at x = 1 and x = -1 (100 mod 101)
        assert_eq!(z.evaluate(FieldElement::new(1)).value, 0);
        assert_eq!(z.evaluate(FieldElement::new(100)).value, 0);
    }

    #[test]
    fn test_lagrange_interpolate() {
        let points = vec![
            (FieldElement::new(0), FieldElement::new(1)),
            (FieldElement::new(1), FieldElement::new(2)),
            (FieldElement::new(2), FieldElement::new(4)),
        ];
        let p = Polynomial::lagrange_interpolate(&points);
        for (x, y) in points {
            assert_eq!(p.evaluate(x.clone()), y);
        }
    }
}
