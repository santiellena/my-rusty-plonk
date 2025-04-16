use crate::my_field::FieldElement;

use std::ops::{Add, Mul, Sub};

/// Note: A polynomial built with the following coefficients:
///     [1, 2, 3, 4] will represent y = 1 + 2x + 3x^2 + 4x^3
///     I want to make this clear because we tend to read polynomials
///     from the highest degree to the lowest degree, but this implementation
///     does as shown avobe. Thanks for reading!

#[derive(Debug, Clone)]
pub struct Polynomial {
    coeffs: Vec<FieldElement>,
    order: u64, // so I make sure all field elements have the same order than the expected in the poly
}

impl Polynomial {
    fn new(coeffs: Vec<FieldElement>, order: u64) -> Self {
        let mut p: Polynomial = Polynomial {
            coeffs: coeffs
                .into_iter()
                .map(|c| FieldElement::new(c.value, order))
                .collect(),
            order,
        };
        p.trim();
        p
    }

    fn trim(&mut self) {
        while self.coeffs.len() > 1
            && self.coeffs.last().unwrap() == &FieldElement::zero(self.order)
        {
            self.coeffs.pop();
        }
    }

    #[allow(dead_code)]
    fn evaluate(&self, x: FieldElement) -> FieldElement {
        assert_eq!(
            x.order, self.order,
            "Evaluation point must be in the same field"
        );
        self.coeffs
            .iter()
            .rev()
            .fold(FieldElement::zero(self.order), |acc, coeff| {
                acc.multiply(&x).add(coeff.clone())
            })
    }

    #[allow(dead_code)]
    fn scalar_mul(&self, scalar: FieldElement) -> Self {
        assert_eq!(scalar.order, self.order, "Scalar must be in the same field");
        Polynomial {
            coeffs: self.coeffs.iter().map(|c| c.multiply(&scalar)).collect(),
            order: self.order,
        }
    }

    #[allow(dead_code)]
    pub fn divide(&self, divisor: &Self) -> (Self, Self) {
        assert_eq!(
            self.order, divisor.order,
            "Polynomials must be over the same field"
        );
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
                Polynomial::new(quotient_coeffs, self.order),
                Polynomial::new(vec![], self.order), // Remainder is 0
            );
        }

        let mut dividend: Vec<FieldElement> = self.coeffs.clone();
        let divisor_deg: usize = divisor.coeffs.len() - 1;
        let divisor_lead: &FieldElement = divisor.coeffs.last().unwrap(); // Leading coefficient
        let inv_divisor_lead: FieldElement = divisor_lead.inverse();

        // If dividend degree < divisor degree, quotient is 0, remainder is dividend
        if dividend.len() <= divisor_deg {
            return (
                Polynomial::new(vec![], self.order),
                Polynomial::new(dividend, self.order),
            );
        }

        let mut quotient: Vec<FieldElement> =
            vec![FieldElement::zero(self.order); dividend.len() - divisor_deg];
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
        while dividend.len() > divisor_deg
            && dividend.last().unwrap() == &FieldElement::zero(self.order)
        {
            dividend.pop();
        }

        (
            Polynomial::new(quotient, self.order),
            Polynomial::new(dividend, self.order),
        )
    }

    #[allow(dead_code)]
    pub fn degree(&self) -> usize {
        if self.coeffs.is_empty()
            || self
                .coeffs
                .iter()
                .all(|c| c == &FieldElement::zero(self.order))
        {
            return 0; // Convention: zero polynomial has degree 0 (or -∞ in some contexts)
        }
        self.coeffs.len() - 1
    }

    /// This vanishing polynomial calculation instead of the classic: ∏(x -y_i)
    /// is possible because we are working in a cyclic subgroup.
    /// This way is more efficient.
    #[allow(dead_code)]
    pub fn vanishing_polynomial(n: usize, order: u64) -> Self {
        let mut coeffs = vec![FieldElement::zero(order); n + 1];
        coeffs[0] = FieldElement::new(1, order).negate(); // -1
        coeffs[n] = FieldElement::one(order); // 1
        Polynomial::new(coeffs, order)
    }

    /// Resource I recommend to understand Lagrange Interpolation:
    /// LambdaClass YT video: https://www.youtube.com/watch?v=REnFOKo9gXs
    ///
    /// P(x) = ∑_i y_i ⋅ l_i(x)
    /// where l_i(x)= ∏_j≠i (x−xj)/(xi−xj)
    #[allow(dead_code)]
    pub fn lagrange_interpolate(points: &[(FieldElement, FieldElement)], order: u64) -> Self {
        assert!(!points.is_empty(), "Need at least one point");
        let mut result = Polynomial::new(vec![], order);

        for (i, &(ref xi, ref yi)) in points.iter().enumerate() {
            let mut term = Polynomial::new(vec![FieldElement::one(order)], order); // Start with 1
            let mut denominator = FieldElement::one(order);

            for (j, &(ref xj, _)) in points.iter().enumerate() {
                if i != j {
                    // Numerator: (x - xj)
                    let num =
                        Polynomial::new(vec![xj.clone().negate(), FieldElement::one(order)], order);
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
        assert_eq!(
            self.order, other.order,
            "Polynomials must be over the same field"
        );
        let max_len: usize = self.coeffs.len().max(other.coeffs.len());
        let mut result: Vec<FieldElement> = Vec::with_capacity(max_len);

        for i in 0..max_len {
            let a: &FieldElement = if i < self.coeffs.len() {
                &self.coeffs[i]
            } else {
                &FieldElement::zero(self.order)
            };
            let b: &FieldElement = if i < other.coeffs.len() {
                &other.coeffs[i]
            } else {
                &FieldElement::zero(self.order)
            };
            result.push(a.add(b));
        }

        Polynomial::new(result, self.order)
    }
}

impl Mul for Polynomial {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        assert_eq!(
            self.order, other.order,
            "Polynomials must be over the same field"
        );
        let n = self.coeffs.len();
        let m = other.coeffs.len();
        let mut result: Vec<FieldElement> = vec![FieldElement::zero(self.order); n + m - 1];

        for i in 0..n {
            for j in 0..m {
                let prod: FieldElement = self.coeffs[i].clone().multiply(&other.coeffs[j]);
                result[i + j] = result[i + j].clone().add(prod);
            }
        }

        Polynomial::new(result, self.order)
    }
}

impl Sub for Polynomial {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        let negated: Polynomial = Polynomial {
            coeffs: other.coeffs.into_iter().map(|c| c.negate()).collect(),
            order: other.order,
        };
        self + negated
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polys() {
        let order = 7; // Example: Field F_7 (prime modulus)

        // P(x) = 1 + 2x + 3x^2 (mod 7)
        let p1 = Polynomial::new(
            vec![
                FieldElement::new(1, order),
                FieldElement::new(2, order),
                FieldElement::new(3, order),
            ],
            order,
        );
        // Q(x) = 5 + 4x (mod 7)
        let p2 = Polynomial::new(
            vec![FieldElement::new(5, order), FieldElement::new(4, order)],
            order,
        );

        // Addition: (1 + 2x + 3x^2) + (5 + 4x) = 6 + 6x + 3x^2 (mod 7)
        let sum = p1.clone() + p2.clone();
        println!("Sum: {:?}", sum.coeffs); // Should be [6, 6, 3]
        assert_eq!(
            sum.coeffs,
            vec![
                FieldElement::new(6, order),
                FieldElement::new(6, order),
                FieldElement::new(3, order)
            ]
        );

        // Multiplication: (1 + 2x + 3x^2)(5 + 4x) = 5 + 14x + 23x^2 + 12x^3 ≡ 5 + 0x + 2x^2 + 5x^3 (mod 7)
        let product = p1.clone() * p2.clone();
        println!("Product: {:?}", product.coeffs); // Should be [5, 0, 2, 5]
        assert_eq!(
            product.coeffs,
            vec![
                FieldElement::new(5, order),
                FieldElement::new(0, order),
                FieldElement::new(2, order),
                FieldElement::new(5, order)
            ]
        );

        // Evaluate P1 at x = 2: 1 + 2*2 + 3*4 = 1 + 4 + 12 = 17 ≡ 3 (mod 7)
        let value = p1.evaluate(FieldElement::new(2, order));
        println!("P1(2) = {:?}", value); // Should be FieldElement { value: 3, order: 7 }
        assert_eq!(value, FieldElement::new(3, order));
    }

    #[test]
    fn test_polynomial_add() {
        let order = 7;
        let p1 = Polynomial::new(
            vec![FieldElement::new(1, order), FieldElement::new(2, order)],
            order,
        );
        let p2 = Polynomial::new(
            vec![FieldElement::new(3, order), FieldElement::new(4, order)],
            order,
        );
        let sum = p1 + p2;
        assert_eq!(
            sum.coeffs,
            vec![FieldElement::new(4, order), FieldElement::new(6, order)]
        );
    }

    #[test]
    fn test_polynomial_divide() {
        let order = 7;
        // P(x) = 1 + 2x + 3x^2
        let p = Polynomial::new(
            vec![
                FieldElement::new(1, order),
                FieldElement::new(2, order),
                FieldElement::new(3, order),
            ],
            order,
        );
        // D(x) = 2 + x
        let d = Polynomial::new(
            vec![FieldElement::new(2, order), FieldElement::new(1, order)],
            order,
        );
        let (q, r) = p.divide(&d);
        // Expected: Q(x) = 3x - 4 (3 mod 7), R(x) = 2 (since 3x^2 + 2x + 1 = (x + 2)(3x - 4) + 2)
        assert_eq!(
            q.coeffs,
            vec![FieldElement::new(3, order), FieldElement::new(3, order)]
        ); // -4 ≡ 3 mod 7
        assert_eq!(r.coeffs, vec![FieldElement::new(2, order)]);
    }

    #[test]
    fn test_polynomial_degree() {
        let order = 7;
        let p = Polynomial::new(
            vec![
                FieldElement::new(1, order),
                FieldElement::new(2, order),
                FieldElement::new(3, order),
            ],
            order,
        );
        assert_eq!(p.degree(), 2);
        let zero = Polynomial::new(vec![], order);
        assert_eq!(zero.degree(), 0);
    }

    #[test]
    fn test_vanishing_polynomial() {
        let order = 7;
        let n = 2; // x^2 - 1
        let z = Polynomial::vanishing_polynomial(n, order);
        assert_eq!(
            z.coeffs,
            vec![
                FieldElement::new(6, order), // -1 ≡ 6 mod 7
                FieldElement::new(0, order),
                FieldElement::new(1, order),
            ]
        );
        // Roots at x = 1 and x = -1 (6 mod 7)
        assert_eq!(z.evaluate(FieldElement::new(1, order)).value, 0);
        assert_eq!(z.evaluate(FieldElement::new(6, order)).value, 0);
    }

    #[test]
    fn test_lagrange_interpolate() {
        let order = 7;
        let points = vec![
            (FieldElement::new(0, order), FieldElement::new(1, order)),
            (FieldElement::new(1, order), FieldElement::new(2, order)),
            (FieldElement::new(2, order), FieldElement::new(4, order)),
        ];
        let p = Polynomial::lagrange_interpolate(&points, order);
        for (x, y) in points {
            assert_eq!(p.evaluate(x.clone()), y);
        }
    }
}
