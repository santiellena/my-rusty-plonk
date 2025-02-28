#[path = "math/field.rs"]
mod field;

use field::Field;

fn main() {
    let x: Field = Field {
        value: 3,
        order: 14,
    };
    let y: Field = x.negate().add(&x);
    println!("{:?}", y);
}
