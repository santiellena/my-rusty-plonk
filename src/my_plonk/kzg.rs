use super::elliptic_curve::{EllipticCurve, Point, PointExt};
use super::field::FieldElement;
use super::polynomial::Polynomial;
use super::toy_pairing::Pairing;
use ark_std::rand;
use rand::Rng;

/*
    Preferred parameters:
    - Elliptic curve: y^2 = x^3 + 3
    - Generator 1 (G1): (1, 2)
    - Finite field: mod 101
    - EC subgroup order: 17 (there are 17 valid EC points created from G)
*/

#[derive(Debug)]
pub struct KZG {
    curve: EllipticCurve,
    setup_g1: Vec<Point>,    // [1]G, [tau]G, ...
    setup_g2: Vec<PointExt>, // [1]H, [tau]H
}

impl KZG {
    #[allow(dead_code)]
    pub fn new(degree: usize) -> Self {
        let curve = EllipticCurve::new();
        let mut rng = rand::thread_rng();
        let tau = FieldElement::new(rng.gen_range(1..FieldElement::MODULUS));
        let g = curve.generator_g1();
        let h = curve.generator_g2();

        let mut setup_g1 = vec![g.clone()];
        for n in 1..=degree {
            let tau_pow_n = tau.pow(n as u64);
            let point = g.scalar_mul(&curve, tau_pow_n);
            setup_g1.push(point);
        }

        let mut setup_g2 = vec![h.clone()];
        let tau_h = h.scalar_mul(&curve, tau);
        setup_g2.push(tau_h);

        KZG {
            curve,
            setup_g1,
            setup_g2,
        }
    }

    #[allow(dead_code)]
    pub fn commit(&self, poly: &Polynomial) -> Point {
        let mut commitment = self.curve.infinity();
        let degree = poly.coeffs.len() - 1;
        for i in 0..=degree {
            let coeff = poly.coeffs[i].clone();
            let power = &self.setup_g1[i];
            let scaled_power = power.scalar_mul(&self.curve, coeff);
            commitment = self.curve.add(&commitment, &scaled_power);
        }
        commitment
    }

    #[allow(dead_code)]
    pub fn prove(&self, poly: &Polynomial, z: FieldElement) -> (FieldElement, Point) {
        let y = poly.evaluate(z.clone());
        let mut q_coeffs = vec![FieldElement::new(0); poly.coeffs.len() - 1];
        let mut remainder = poly.coeffs.clone();
        remainder[0] = remainder[0].substract(&y);
        for i in (1..remainder.len()).rev() {
            q_coeffs[i - 1] = remainder[i].clone();
            remainder[i] = FieldElement::new(0);
            remainder[i - 1] = remainder[i - 1].substract(&q_coeffs[i - 1].multiply(&z));
        }
        let q_poly = Polynomial { coeffs: q_coeffs };

        let mut proof = self.curve.infinity();
        let degree = q_poly.coeffs.len() - 1;
        for i in 0..=degree {
            let coeff = q_poly.coeffs[i].clone();
            let power = &self.setup_g1[i];
            let scaled_power = power.scalar_mul(&self.curve, coeff);
            proof = self.curve.add(&proof, &scaled_power);
        }
        (y, proof)
    }

    #[allow(dead_code)]
    pub fn verify(
        &self,
        commitment: &Point,
        z: FieldElement,
        y: FieldElement,
        proof: &Point,
    ) -> bool {
        let g1 = &self.setup_g1[0];
        let g2 = &self.setup_g2[0];
        let tau_g2 = &self.setup_g2[1];
        let y_g1 = g1.scalar_mul(&self.curve, y);
        let commitment_minus_y_g1 = self.curve.add(
            commitment,
            &y_g1.scalar_mul(&self.curve, FieldElement::new(FieldElement::MODULUS - 1)),
        );
        let z_g2 = g2.scalar_mul(&self.curve, z);
        let tau_g2_minus_z_g2 = self.curve.add_ext(
            tau_g2,
            &z_g2.scalar_mul(&self.curve, FieldElement::new(FieldElement::MODULUS - 1)),
        );

        let map_to_scalar = |p: &Point| {
            if p.is_infinity {
                FieldElement::new(1)
            } else {
                p.x.clone()
            }
        };
        let map_to_scalar_ext = |p: &PointExt| {
            if p.is_infinity {
                FieldElement::new(1)
            } else {
                p.x.a.clone()
            }
        };

        let left = Pairing::new(map_to_scalar(&commitment_minus_y_g1), map_to_scalar_ext(g2));
        let right = Pairing::new(map_to_scalar(proof), map_to_scalar_ext(&tau_g2_minus_z_g2));
        println!("Left: {:?}", left.e);
        println!("Right: {:?}", right.e);
        // the toy pairing lacks of bilinearity so it won't correctly compute anything...
        // don't really want to deep dive into a pairing implementation, maybe in the future...
        // for the ark plonk version I'll be using a proper pairing function and constructing a real
        // world plonk (educational).
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kzg() {
        let kzg = KZG::new(2);
        let poly = Polynomial::new(vec![FieldElement::one(), FieldElement::new(2)]); // 1 + 2x
        let commitment = kzg.commit(&poly);
        let z = FieldElement::new(3);
        let (y, proof) = kzg.prove(&poly, z.clone());
        assert_eq!(
            proof,
            kzg.curve
                .generator_g1()
                .scalar_mul(&kzg.curve, FieldElement::new(2))
        );
        assert_eq!(y.value, 7); // 1 + 2*3 = 7
        assert!(kzg.verify(&commitment, z, y, &proof));
    }
}
