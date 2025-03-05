#[path = "./field.rs"]
mod field_element;
use field_element::FieldElement;

use std::ops::{Add, Mul};

#[derive(Debug, Clone)]
pub struct Polynomial {
    coeffs: Vec<FieldElement>,
    order: u32, // so I make sure all field elements have the same order than the expected in the poly
}

impl Polynomial {
    fn new(coeffs: Vec<FieldElement>, order: u32) -> Self {
        let mut p = Polynomial {
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
}

impl Add for Polynomial {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        assert_eq!(
            self.order, other.order,
            "Polynomials must be over the same field"
        );
        let max_len = self.coeffs.len().max(other.coeffs.len());
        let mut result = Vec::with_capacity(max_len);

        for i in 0..max_len {
            let a = if i < self.coeffs.len() {
                self.coeffs[i].clone()
            } else {
                FieldElement::zero(self.order)
            };
            let b = if i < other.coeffs.len() {
                other.coeffs[i].clone()
            } else {
                FieldElement::zero(self.order)
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
        let mut result = vec![FieldElement::zero(self.order); n + m - 1];

        for i in 0..n {
            for j in 0..m {
                let prod = self.coeffs[i].clone().multiply(&other.coeffs[j]);
                result[i + j] = result[i + j].clone().add(prod);
            }
        }

        Polynomial::new(result, self.order)
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
}
