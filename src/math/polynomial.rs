#[derive(Debug, Clone)]
struct Polynomial {
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
                acc.multiply(&x).add(coeff)
            })
    }
}
