/*pub fn gcd(mut x: u32, mut y: u32) -> u32 {
    while y != 0 {
        let temp = y;
        y = x % y;
        x = temp;
    }
    x
}*/

pub fn ext_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    let (mut old_r, mut r) = (a, b);
    let (mut old_s, mut s) = (1, 0);
    let (mut old_t, mut t) = (0, 1);

    while r != 0 {
        let quotient = old_r / r;
        (old_r, r) = (r, old_r - quotient * r);
        (old_s, s) = (s, old_s - quotient * s);
        (old_t, t) = (t, old_t - quotient * t);
    }

    (old_s, old_t, old_r) // Return (s, t, gcd)
}
