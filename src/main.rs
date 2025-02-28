#[path = "math/field.rs"]
mod field;

use field::Field;

fn main() {
    let x: Field = Field {
        value: 240,
        order: 14,
    };

    println!("{:?}", x.pow(262));
}
