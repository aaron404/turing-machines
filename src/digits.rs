use image::{self, ImageBuffer, Luma};
pub const SCALE: u32 = 3;

pub fn load_digits() -> Vec<ImageBuffer<Luma<u8>, Vec<u8>>> {
    let img =
        image::open(format!("assets/digits{SCALE}x.png")).expect("Failed to open digits image");
    let img = img.to_luma8();
    let mut digits = Vec::new();
    for i in 0..10 {
        let x = i;
        let digit =
            image::imageops::crop_imm(&img, x * 4 * SCALE, 0, 3 * SCALE, 5 * SCALE).to_image();
        digits.push(digit);
    }
    digits
}
