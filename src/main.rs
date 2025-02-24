#[path = "math/field.rs"]
mod field;

use field::Field;

fn main() {
    let x: Field = Field { value: 3, order: 7 };

    println!("{:?}", x.inverse());
}
