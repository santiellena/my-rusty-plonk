#[path = "ext_euclidean_algo.rs"]
mod gcd;

#[derive(Clone, Debug)]
pub struct FieldElement {
    pub value: u32,
    pub order: u32,
}

// use PartialEq trait to check whether two fields are equal or not
impl PartialEq for FieldElement {
    fn eq(&self, other: &FieldElement) -> bool {
        self.value == other.value && self.order == other.order
    }
}

impl FieldElement {
    /// Constructor
    pub fn new(value: u32, order: u32) -> Self {
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
    fn mod_inverse(&self) -> Option<u32> {
        let (s, gcd) = gcd::ext_gcd(self.value as i32, self.order as i32); // Compute s such that s * a + t * m = gcd(a, m)

        if gcd != 1 {
            return None; // No modular inverse exists if value and order are not coprime
        }

        // Ensure the result is in the range [0, order)
        Some(((s % self.order as i32 + self.order as i32) % self.order as i32) as u32)
    }

    /// Modular inverse (for division)
    pub fn inverse(&self) -> Self {
        if let Some(inv) = Self::mod_inverse(self) {
            Self {
                value: inv,
                order: self.order,
            }
        } else {
            Self::zero(self.order)
        }
    }

    /// Division (a / b = a * b⁻¹ mod p)
    pub fn divide(&self, b: u32) -> Self {
        let b_field: FieldElement = FieldElement {
            value: b,
            order: self.order,
        };

        let inverse_b: &FieldElement = &FieldElement {
            value: b_field.inverse().value,
            order: b_field.order,
        };

        self.multiply(inverse_b)
    }

    /// Modular exponentiation (a^exp mod p) using fast exponentiation
    pub fn pow(&self, exp: u32) -> Self {
        let mut fsb_found: bool = false; // first significant byte found
        let mut result: Self = self.clone();

        for n in (0..32).rev() {
            let bin: u32 = (exp >> n) & 1;
            print!("{}", bin);
            if (fsb_found) {
                if bin == 1 {
                    result = result.multiply(&result).multiply(&Self {
                        value: 2,
                        order: self.order,
                    });
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

    /// Negation in the FieldElement
    pub fn negate(&self) -> Self {
        Self::new(self.order, self.order).substract(&self)
    }

    /// Zero element
    pub fn zero(order: u32) -> Self {
        Self::new(0, order)
    }

    /// One element
    pub fn one(order: u32) -> Self {
        Self::new(1, order)
    }
}
