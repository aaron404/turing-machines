#![feature(cfg_boolean_literals)]
#![feature(random)]
#![feature(slice_as_chunks)]

use core::f32;
use image::{ImageBuffer, Luma};
use rle::print_rle;
use std::collections::HashMap;
use std::env;
use std::fmt::Display;
use std::fs::read_to_string;
use std::ops::Shl;
use std::sync::atomic::AtomicU32;
use std::sync::atomic::Ordering::Relaxed;
use std::time::Instant;

mod beaver;
mod rle;
mod types;

use beaver::Beaver;
use beaver::Progress;
use types::*;

struct Tape {
    left: Vec<Bin>,
    right: Vec<Bin>,
}

impl Tape {
    fn new() -> Self {
        Self {
            left: vec![0],
            right: vec![0],
        }
    }

    fn pos_to_bin(&self, i: i32) -> &Bin {
        if i < 0 {
            &self.left[(i.abs() - 1) as usize / SIZE]
        } else {
            &self.right[i as usize / SIZE]
        }
    }

    fn get(&self, i: i32) -> u8 {
        let s = SIZE as i32;
        let bin = if i < 0 { (i + 1) / s - 1 } else { i / s };
        let bit = i.rem_euclid(s);
        let side = if bin < 0 { &self.left } else { &self.right };

        let bin = if bin < 0 { -bin - 1 } else { bin } as usize;
        if bin >= side.len() {
            return 0;
        }
        (side[bin] >> bit) as u8 & 1
    }

    fn set(&mut self, i: i32, v: bool) {
        let s = SIZE as i32;
        let bin = if i < 0 { (i + 1) / s - 1 } else { i / s };
        let bit = i.rem_euclid(s);
        let side = if bin < 0 {
            &mut self.left
        } else {
            &mut self.right
        };

        let bin = if bin < 0 { -bin - 1 } else { bin } as usize;
        while bin >= side.len() {
            side.push(0);
        }

        if v {
            side[bin] |= 1 << bit;
        } else {
            side[bin] &= !(1 << bit);
        }
    }

    fn get_neighborhood(&self, i: i32) -> CacheBucket {
        let mut n = 0;
        let l_pos = i - CACHE_BUCKET_HALF_SIZE as i32;
        let left_span = 0;

        for j in 0..CACHE_BUCKET_SIZE {
            let index = i + j as i32 - CACHE_BUCKET_HALF_SIZE as i32;
            n |= (self.get(index) as CacheBucket) << j;
        }
        n
    }

    fn get_neighborhood2(&self, i: i32) -> CacheBucket {
        let l_pos = i - CACHE_BUCKET_HALF_SIZE as i32;
        let r_pos = l_pos + CACHE_BUCKET_SIZE as i32;
        println!("l_pos: {l_pos}, r_pos: {r_pos}");

        let l_bin = self.pos_to_bin(l_pos);
        let r_bin = self.pos_to_bin(r_pos);
        println!("l_bin: {:016b}, r_bin: {:016b}", l_bin, r_bin);

        if l_bin == r_bin {
            let start_bit = l_pos.rem_euclid(SIZE as i32) as i8;
            let shift = SIZE - start_bit as usize - CACHE_BUCKET_SIZE;
            return (l_bin >> shift) as CacheBucket;
        } else {
            let start_bit = l_pos.rem_euclid(SIZE as i32) as i8;
            let end_bit = r_pos.rem_euclid(SIZE as i32) as i8;
            println!("start_bit: {start_bit}, end_bit: {end_bit}");
            let low_bits = r_bin.reverse_bits() >> (SIZE - end_bit as usize);
            // let low_bits = (low_bits & (Bin::MAX >> (SIZE - end_bit as usize))) as CacheBucket;
            println!("low_bits: {:016b}", low_bits);
            let high_bits = l_bin;
            // let low_bits = signed_shift(r_bin, CACHE_BUCKET_SIZE as i8 - start_bit); // off by one maybe
            // let low_bits = (low_bits & (Bin::MAX >> end_bit as usize)) as CacheBucket;
            // let high_bits = signed_shift(l_bin, CACHE_BUCKET_SIZE as i8 - start_bit);
            // let high_bits = (high_bits & (Bin::MAX << start_bit as usize)) as CacheBucket;
            // return low_bits | high_bits;
            return 0;
        }
    }

    fn get_strings(&self, pos: i32) -> (String, String) {
        let mut left_str = String::new();
        let mut right_str = String::new();

        let left_extent = self.left.len() * SIZE;
        let right_extent = self.right.len() * SIZE;

        for i in (-(left_extent as isize)..pos as isize).rev() {
            let bit = self.get(i as i32);
            left_str.push(if bit == 0 { '0' } else { '1' });
        }
        for i in pos as isize..right_extent as isize {
            let bit = self.get(i as i32);
            right_str.push(if bit == 0 { '0' } else { '1' });
        }

        (
            left_str.trim_end_matches('0').to_string(),
            right_str.trim_end_matches('0').to_string(),
        )
    }
}

impl Display for Tape {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in (0..self.left.len()).rev() {
            for j in 0..SIZE {
                write!(f, "{}", (self.left[i] >> j) & 1)?;
            }
        }
        write!(f, " ")?;
        for i in 0..self.right.len() {
            for j in 0..SIZE {
                write!(f, "{}", (self.right[i] >> j) & 1)?;
            }
        }
        Ok(())
    }
}

use rayon::prelude::*;

mod digits;
use digits::*;

fn main() {
    // for i in -32..=32i16 {
    //     // println!("{}: {}", i, i.rem_euclid(16));
    //     println!("  {}: {:016b}", i, 0b0000000100000000u16 << i)
    // }

    let args = env::args().collect::<Vec<_>>();
    let holdouts = read_to_string(args.get(1).expect("missing holdouts file"))
        .expect("failed to read holdouts file");
    let holdouts = holdouts.lines().collect::<Vec<&str>>();

    let max_steps = args
        .get(2)
        .expect("missing max steps")
        .parse::<u64>()
        .expect("max_steps must be an integer");

    // const w: u32 = 128;
    // const h: u32 = 512;
    // for (i, holdout) in holdouts.iter().take(25).enumerate() {
    //     draw(holdout, w, h, 0, format!("out_{:03}.png", i));
    // }

    // enumerate_beavers();
    holdout_thumbnails(holdouts);
    // run_parallel(holdouts, max_steps);
    // for holdout in holdouts.iter().take(100) {
    //     run_cached(holdout, max_steps);
    // }
    // rle_test2(holdouts);
}

#[allow(unused)]
fn holdout_thumbnails(holdouts: Vec<&str>) {
    const OUT_DIR: &str = "out/thumbs";

    const BORDER: u32 = 1;
    const NUM_PAGES: u32 = 50;
    const TW: u32 = 240;
    const TH: u32 = 128;
    const SKIP: u32 = 100;
    const W: u32 = 2460;
    const H: u32 = 1340;
    const NX: u32 = W / (TW + BORDER);
    const NY: u32 = H / (TH + BORDER);
    const THUMBS_PER_PAGE: u32 = NX * NY;
    println!("{}x{} grid, {}x{} thumbnails", NX, NY, TW, TH);
    // clear output directory
    std::fs::remove_dir_all(OUT_DIR).ok();
    std::fs::create_dir_all(OUT_DIR).expect("Failed to create output directory");

    let digits = digits::load_digits();

    let t0 = Instant::now();
    for page in 0..NUM_PAGES {
        let mut grid = image::ImageBuffer::from_pixel(NX * TW - 1, NY * TH - 1, Luma([128u8]));
        for y in 0..NY {
            for x in 0..NX {
                let i = page * THUMBS_PER_PAGE + y * NX + x;
                if i >= holdouts.len() as u32 {
                    break;
                }
                let holdout = holdouts[i as usize];
                let thumb = Beaver::from_string(&holdout.to_string()).thumbnail(TW, TH, SKIP);

                let offset_x = (x * (TW + BORDER)) as i64;
                let offset_y = (y * (TH + BORDER)) as i64;
                image::imageops::overlay(&mut grid, &thumb, offset_x, offset_y);

                for (k, digit) in i.to_string().chars().enumerate() {
                    let digit = digits
                        .get((digit as u8 - b'0') as usize)
                        .expect("Invalid digit");
                    image::imageops::overlay(
                        &mut grid,
                        digit,
                        1 + offset_x + k as i64 * 4 * digits::SCALE as i64,
                        1 + offset_y,
                    );
                }
            }
        }
        let index = page * THUMBS_PER_PAGE;
        grid.save(format!("{OUT_DIR}/beavers_{:03}.png", index))
            .expect("Failed to save grid image");
    }

    let elapsed = t0.elapsed().as_secs_f32();
    println!(
        "Generated {} pages of thumbnails in {:.2}s",
        NUM_PAGES, elapsed
    );
    println!("  {} pages per sec", NUM_PAGES as f32 / elapsed);
    println!(
        "  {} thumbs per sec",
        (THUMBS_PER_PAGE * NUM_PAGES) as f32 / elapsed
    );
}

#[allow(unused)]
fn enumerate_beavers() {
    let mut beavers = Beaver::enumerate(3);
    println!("enumerated {} beavers", beavers.len());
    let max_steps = 10;
    for (i, b) in beavers.iter_mut().enumerate() {
        // println!("{}: {}", i, b.get_rule());
        // b.thumbnail(32, 128, format!("beaver_{:03}.png", i));
        let mut count = 0;
        let mut halted = false;
        loop {
            if !b.step() {
                // println!("  halted at step {}", count);
                halted = true;
                break;
            }
            count += 1;
            if count >= max_steps {
                // println!("  did not halt after {} steps", count);
                break;
            }
        }
    }
}

#[allow(unused)]
fn draw(holdout: &str, w: u32, h: u32, skip: u64, fname: String) {
    println!("holdout: {:?}", holdout);
    let mut b = Beaver::from_string(&holdout.to_string());
    let mut img = image::GrayImage::new(w, h);
    let max_steps = 10000000u64;
    let mut count = 0u64;
    loop {
        count += 1;

        if !b.draw(&mut img) {
            break;
        }

        if !b.step() || count >= h as u64 {
            // println!("halted at step {}", count);
            break;
        }
    }
    img.save(fname).unwrap();
}

fn run_cached(holdout: &str, max_steps: u64) {
    let mut b = Beaver::from_string(&holdout.to_string());
    let mut steps = 0u64;
    let mut halted = true;

    let t0 = Instant::now();
    loop {
        match b.step_with_cache() {
            Progress::Halted(s) => {
                steps += s as u64;
                // println!("tape@{}: {}", steps, b.tape);
                halted = true;
                break;
            }
            Progress::DidNotHalt(s) => {
                steps += s as u64;
                // println!("tape@{}: {}", steps, b.tape);
                if steps > max_steps {
                    break;
                }
            }
        }
    }

    let elapsed = t0.elapsed().as_secs_f32();
    print!("holdout: {:?} ", holdout);
    print!(
        "{: >7.2} M steps/sec ",
        steps as f32 / 1_000_000.0 / elapsed
    );
    if halted {
        println!("(halted at step {})", steps);
    } else {
        println!("(did not halt after {} steps) ", steps);
    }
}

#[allow(unused)]
fn run_parallel(holdouts: Vec<&str>, max_steps: u64) {
    let halt_count = AtomicU32::new(0);
    let processed_count = AtomicU32::new(0);

    holdouts.par_iter().for_each(|holdout| {
        // println!("holdout: {:?}", holdout);

        let mut b = Beaver::from_string(&holdout.to_string());
        let mut step = 0;
        while b.step() {
            step += 1;
            if step >= max_steps {
                break;
            }
        }
        processed_count.fetch_add(1, Relaxed);
        if step == max_steps {
            // println!("  did not halt");
        } else {
            let h = halt_count.fetch_add(1, Relaxed) + 1;
            let p = processed_count.fetch_add(0, Relaxed);
            let percent = h as f32 / p as f32;
            println!("  [{}/{}] ({}) halted at step {}", h, p, percent, step);
        }
    });

    println!("halt count: {}", halt_count.fetch_add(0, Relaxed));
}

use std::random::random;
#[allow(unused)]
fn experiment() {
    // let mut b = Beaver::new("1RB1LC_1RC1RB_1RD0LE_1LA1LD_---0LA".to_string());
    // let mut b = Beaver::new("1RB1RD_1RC1RE_0LA0RC_1LB---_0LD1LE".to_string());
    // let mut b = Beaver::new("1RB0LD_1RC0RA_1RD0RF_1LE1LA_1RA1LE_0RB---".to_string()); // bb6 antihydra
    // let mut b = Beaver::new("1RB---_1LC1RA_1RD1LC_1LF1RE_1RC0RD_0LB0LF".to_string());
    // let mut b = Beaver::new("1RB0LD_0RC0RD_1RD0LA_1RE1RE_1LF1RB_1RZ1LG_1LC1LE".to_string());
    let mut b =
        Beaver::from_string(&"1RB1RG_1RC0RB_1RD0RE_1RE0LD_1LF1RA_---1LD_1LC1RG".to_string());

    // let mut img = image::GrayImage::new(256, 512);
    let mut img = image::GrayImage::new(1024, 1024);

    let max_steps = 10000000u64;
    // let max_steps = 10000000u64;
    // let mut snapshot_count = 500000;
    let mut snapshot_count = 100u64;
    let t0 = Instant::now();
    let mut count = 0u64;
    loop {
        count += 1;

        // if count == snapshot_count {
        //     if !b.draw(&mut img) {
        //         // if !b.draw_centered(&mut img) {
        //         // snapshot_count = u64::MAX;
        //         break;
        //     }
        //     // snapshot_count = ((snapshot_count * 1000) / 990).max(snapshot_count + 1);
        //     snapshot_count += 1;
        //     // snapshot_count *= 2;
        // }

        if !b.step() || count >= max_steps {
            println!("halted at step {}", count);
            break;
        }
    }
    let elapsed = t0.elapsed().as_secs_f32();

    img.save("out.png").unwrap();

    // println!("position: {}", b.position);
    println!(
        "{} steps in {:.2}s. {:.2} steps/s",
        count,
        elapsed,
        count as f32 / elapsed
    );

    // rle test
    // rle_test();

    // rle_bits_test();
}

#[allow(unused)]
fn rle_test(holdouts: Vec<&str>) {
    const N: usize = 100000;
    let mut data = vec![0u8; N];
    for i in 0..N {
        data[i] = b'0' + (random::<u8>() % 2) as u8;
    }

    let data = "B 111111110011111100111111001100110011111001111100111111001111111110011111111100111110011111100111110011111001100111110011001111100110011111001111110011111111100111111111001111111110011111111100111111111001111111110011111111100111111111001111111110011111001111110011111001111100110011111001100111110011001111100110011111001100111110011001111100110011111001100111110011001111100110011111001100111110011001111100110011111001100111110011001111100110011111001100111110011001111100110011111001111110011111111100111111111001111111110011111111100111111111001111111110011111111100111111111001111111110011111111100111111111001111111110011111111100111111111001111111110011111111100111111111001111100 >1 010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010101001010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101010010101010100101010101001010101001".as_bytes().to_vec();
    // data = "00000001000".as_bytes().to_vec();

    println!("data: {}", String::from_utf8_lossy(&data));

    let t0 = Instant::now();
    let res = rle::rle(&data);

    let elapsed = t0.elapsed().as_secs_f32();
    println!(
        "RLE: {} bits in {:.2}s. {:.2} steps/s",
        data.len(),
        elapsed,
        res.len() as f32 / elapsed
    );

    rle::print_rle(res);
}

#[allow(unused)]
fn rle_test2(holdouts: Vec<&str>) {
    let holdout = holdouts[0];
    println!("holdout: {:?}", holdout);
    let mut b = Beaver::from_string(&holdout.to_string());
    for i in 0..100000 {
        if !b.step() {
            println!("halted at step {}", i);
            break;
        }
    }

    let (l, r) = b.get_tape_strings();
    println!("beaver: {}", b.tape);
    println!("left: {}", l);
    println!("right: {}", r);

    print!("left RLE: ");
    print_rle(rle::rle(l.as_bytes()));
    print!("right RLE: ");
    print_rle(rle::rle(r.as_bytes()));
}

#[allow(unused)]
fn rle_bits_test() {
    const N: usize = 6;
    const NUM_BINS: usize = (N - 1) / SIZE + 1;

    let mut data: Vec<Bin> = vec![0; NUM_BINS];
    for i in 0..N {
        data[i / SIZE] |= (random::<Bin>() % 2) << (i % SIZE);
    }

    for i in 0..NUM_BINS {
        println!("data[{}]: {:064b}", i, data[i]);
    }

    let t0 = Instant::now();
    // let res = rle::rle_bits(&data.to_vec(), N);
    let elapsed = t0.elapsed().as_secs_f32();
}
