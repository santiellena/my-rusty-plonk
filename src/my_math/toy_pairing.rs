use crate::my_field::FieldElement;

/*
    Disclaimer: This is obviously a toy pairing implementation.
    Given my lack of knowledge on actual pairings and the time it
    will take me to understand them to then implement some simple
    version, I decided to keep this as simple as possible.

    The check still works, but it lacks of key security features of
    real pairings (affecting soundness):
    - Outputs a F_101 element ignoring the embedding degree (k = 2)
    - Lacks bilinearity ( e([a]P, b[Q]) != e(P, Q)^ab )
*/

#[derive(Clone, Debug, PartialEq)]
pub struct Pairing {
    pub e: FieldElement,
}

impl Pairing {
    pub fn new(a: FieldElement, b: FieldElement) -> Self {
        Pairing { e: a.multiply(&b) }
    }

    #[allow(dead_code)]
    pub fn mul(&self, other: &Self) -> Self {
        Pairing {
            e: self.e.multiply(&other.e),
        }
    }

    #[allow(dead_code)]
    pub fn pow(&self, n: u64) -> Self {
        let mut result = Pairing {
            e: FieldElement::new(1),
        };
        let mut base = self.clone();
        let mut exp = n;
        while exp > 0 {
            if exp & 1 == 1 {
                result = result.mul(&base);
            }
            base = base.mul(&base);
            exp >>= 1;
        }
        result
    }
}
