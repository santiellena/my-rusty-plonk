use crate::my_elliptic_curve::{EllipticCurve, Point};
use crate::my_field::FieldElement;
use crate::my_polynomial::Polynomial;

#[allow(dead_code)]
pub struct KZG {
    curve: EllipticCurve,
    srs: Vec<Point>, // [1]G, [tau]G, [tau^2]G ...
}

/*
    Preferred parameters:
    - Elliptic curve: y^2 = x^3 + 3
    - Generator 1 (G1): (1, 2)
    - Finite field: mod 101
    - EC subgroup order: 17 (there are 17 valid EC points created from G)
*/

impl KZG {}
