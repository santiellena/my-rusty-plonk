#[path = "math/field.rs"]
mod field_element;

use field_element::FieldElement;

fn main() {
    let x: FieldElement = FieldElement {
        value: 3,
        order: 14,
    };
    let y: FieldElement = x.negate().add(&x);
    println!("{:?}", y);
}
