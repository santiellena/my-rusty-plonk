use super::field::FieldElement;
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

    /// This vanishing polynomial calculation instead of the classic: ∏(x -y_i)
    /// is possible because we are working in a cyclic subgroup.
    /// This way is more efficient.
    #[allow(dead_code)]
    pub fn vanishing_polynomial(n: usize) -> Self {
        let mut coeffs = vec![FieldElement::zero(); n + 1];
        coeffs[0] = FieldElement::one().negate(); // -1
        coeffs[n] = FieldElement::one(); // 1
        Polynomial::new(coeffs)
    }

    #[allow(dead_code)]
    pub fn scalar_mul(&self, scalar: FieldElement) -> Self {
        let coeffs: Vec<FieldElement> = self.coeffs.iter().map(|c| c.multiply(&scalar)).collect();
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

    #[allow(dead_code)]
    pub fn evaluate(&self, x: FieldElement) -> FieldElement {
        self.coeffs
            .iter()
            .rev()
            .fold(FieldElement::zero(), |acc, coeff| {
                acc.multiply(&x).add(coeff.clone())
            })
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

    #[test]
    fn test_vanishing_polynomial() {
        let n = 2; // x^2 - 1
        let z = Polynomial::vanishing_polynomial(n);
        assert_eq!(
            z.coeffs,
            vec![
                FieldElement::one().negate(),
                FieldElement::zero(),
                FieldElement::new(1),
            ]
        );
        // Roots at x = 1 and x = -1
        assert_eq!(z.evaluate(FieldElement::one()), FieldElement::zero());
        assert_eq!(
            z.evaluate(FieldElement::one().negate()),
            FieldElement::zero()
        );
    }

    #[test]
    fn test_lagrange_interpolate() {
        let points = vec![
            (FieldElement::zero(), FieldElement::one()),
            (FieldElement::one(), FieldElement::new(2)),
            (FieldElement::new(2), FieldElement::new(4)),
        ];
        let p = Polynomial::lagrange_interpolate(&points);
        for (x, y) in points {
            assert_eq!(p.evaluate(x.clone()), y);
        }
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
        let d = Polynomial::new(vec![FieldElement::new(2), FieldElement::one()]);
        let (q, r) = p.divide(&d);
        // Expected: Q(x) = 3x - 4, R(x) = 9
        // Since 3x^2 + 2x + 1 = (x + 2)(3x - 4) + 9
        assert_eq!(
            q.coeffs,
            vec![FieldElement::new(4).negate(), FieldElement::new(3)]
        ); // -4 = 13 mod 17 if p=17, else adjust
        assert_eq!(r.coeffs, vec![FieldElement::new(9)]);
    }
}
