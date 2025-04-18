/// MY OWN MATH MODULES

#[path = "my_math/polynomial.rs"]
mod my_polynomial;

#[path = "my_math/elliptic_curve.rs"]
mod my_elliptic_curve;

#[path = "my_math/field.rs"]
mod my_field;

#[path = "my_math/ext_euclidean_algo.rs"]
mod gcd;

#[path = "my_math/toy_pairing.rs"]
mod my_toy_pairing;

#[path = "my_math/kzg.rs"]
mod my_kzg;

/// MATH WITH ARKWORKS MODULES

#[path = "math_with_ark/field.rs"]
mod field;

#[path = "math_with_ark/polynomial.rs"]
mod polynomial;

#[path = "math_with_ark/elliptic_curve.rs"]
mod elliptic_curve;

fn main() {}
