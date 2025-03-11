#[path = "my_math/polynomial.rs"]
mod own_polynomial;

#[path = "my_math/elliptic_curve.rs"]
mod own_elliptic_curve;

#[path = "my_math/field.rs"]
mod own_field;

#[path = "my_math/ext_euclidean_algo.rs"]
mod gcd;

#[path = "math_with_ark/field.rs"]
mod field;

#[path = "math_with_ark/polynomial.rs"]
mod polynomial;

#[path = "math_with_ark/elliptic_curve.rs"]
mod elliptic_curve;

fn main() {}
