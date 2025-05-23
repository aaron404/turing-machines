#![feature(cfg_boolean_literals)]
#![feature(random)]

use core::f32;
use image::{ImageBuffer, Luma};
use std::collections::HashMap;
use std::env;
use std::fmt::Display;
use std::fs::read_to_string;
use std::ops::Shl;
use std::sync::atomic::AtomicU32;
use std::sync::atomic::Ordering::Relaxed;
use std::time::Instant;

mod rle;
mod types;

use types::*;

struct Tape {
    left: Vec<Bin>,
    right: Vec<Bin>,
}

fn signed_shift(x: CacheBucket, shift: i8) -> CacheBucket {
    if shift < 0 { x << shift } else { x >> shift }
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
        let side = if i < 0 { &self.left } else { &self.right };
        let index = if i < 0 {
            (i.abs() - 1) as usize
        } else {
            i as usize
        };
        let bin = index / SIZE;
        let bit = index % SIZE;

        if bin >= side.len() {
            return 0;
        }

        (side[bin] >> bit) as u8 & 1
    }

    fn set(&mut self, i: i32, v: bool) {
        let side = if i < 0 {
            &mut self.left
        } else {
            &mut self.right
        };
        let index = if i < 0 {
            (i.abs() - 1) as usize
        } else {
            i as usize
        };
        let bin = index / SIZE;
        let bit = index % SIZE;

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
}

impl Display for Tape {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in (0..self.left.len()).rev() {
            for j in 0..SIZE {
                write!(f, "{}", (self.left[i] >> (SIZE - j - 1)) & 1)?;
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

#[derive(Clone, Copy, Debug)]
struct Transition(u8);
impl Transition {
    fn from_str(s: &str) -> Self {
        let s = s.as_bytes();
        if s[0] == b'-' {
            return Transition(0xff);
        }
        let write = (s[0] == b'1') as u8;
        let step = (s[1] == b'R') as u8;
        let next_state = s[2] - b'A';
        Transition(write | (step << 1) | (next_state << 2))
    }

    fn write(self) -> bool {
        self.0 & 1 > 0
    }

    fn step(self) -> i8 {
        (self.0 & 2) as i8 - 1
    }

    fn next_state(self) -> u8 {
        (self.0 & 0b01111100) >> 2
    }
}

struct CacheResult {
    steps: u16,
    new_tape: CacheBucket,
    offset: i8,
    new_state: u8,
}

struct Beaver {
    state: u8,
    cache: Vec<HashMap<CacheBucket, CacheResult>>,
    position: i32,
    transitions: Vec<[Transition; 2]>,
    tape: Tape,
    row: u32,
}

enum Progress {
    Halted(u16),
    DidNotHalt(u16),
}

const MAX_CACHED_STEPS: u16 = u16::MAX;

impl Beaver {
    fn new(rule: &String) -> Self {
        let mut transitions = Vec::new();
        let mut cache = Vec::new();
        for t in rule.split("_") {
            let (left, right) = t.split_at(3);
            transitions.push([Transition::from_str(left), Transition::from_str(right)]);
            cache.push(HashMap::new()); // one cache per state
        }
        Self {
            cache: cache,
            state: 0,
            position: 0,
            transitions: transitions,
            tape: Tape::new(),
            row: 0,
        }
    }

    fn step(&mut self) -> bool {
        let side = if self.position < 0 {
            &mut self.tape.left
        } else {
            &mut self.tape.right
        };

        let pos = if self.position < 0 {
            (self.position.abs() - 1) as usize
        } else {
            self.position as usize
        };

        let bin = pos / SIZE;
        let bit = pos % SIZE;

        if bin >= side.len() {
            side.push(0);
        }

        let tape_state = (side[bin] >> bit) & 1;
        let transition = self.transitions[self.state as usize][tape_state as usize];

        if transition.0 == 0xff {
            return false;
        } else {
            self.state = transition.next_state();
        }

        if transition.write() {
            side[bin] |= 1 << bit;
        } else {
            side[bin] &= !(1 << bit);
        }

        self.position += transition.step() as i32;

        true
    }

    fn step_with_cache(&mut self) -> Progress {
        // get neighborhood
        let n = self.tape.get_neighborhood(self.position);
        // check cache
        if let Some(cached) = self.cache[self.state as usize].get(&n) {
            // println!("cache hit");
            // println!("  advancing {} steps", cached.steps);
            // println!(
            //     "  moving from {} to {}",
            //     self.position,
            //     self.position + cached.offset as i64
            // );
            // println!("  old_tape: {:016b}", n);
            // println!("  new_tape: {:08b}", cached.new_tape);
            // println!("  offset: {}", cached.offset);
            // println!("  new_state: {}", cached.new_state);

            let l_pos = self.position - CACHE_BUCKET_HALF_SIZE as i32;

            for i in 0..CACHE_BUCKET_SIZE {
                // println!(
                //     "    setting {} to {}",
                //     l_pos + i,
                //     (cached.new_tape >> i) & 1 > 0
                // );
                self.tape
                    .set(l_pos + i as i32, (cached.new_tape >> i) & 1 > 0);
            }
            self.position += cached.offset as i32;
            self.state = cached.new_state;
            return Progress::DidNotHalt(cached.steps);
        } else {
            // run normally until we leave the neighborhood
            let mut steps = 0;
            let start_pos = self.position;
            let start_state = self.state;
            while self.position.abs_diff(start_pos) < CACHE_BUCKET_SIZE as u32 / 2
                && steps < MAX_CACHED_STEPS
            {
                // println!("  [{}] self.position: {}", steps, self.position);
                if self.step() {
                    steps += 1;
                } else {
                    // halted
                    return Progress::Halted(steps);
                }
            }
            let mut new_tape = 0;
            let l_pos = start_pos - CACHE_BUCKET_HALF_SIZE as i32;

            for i in 0..CACHE_BUCKET_SIZE {
                if self.tape.get(l_pos + i as i32) > 0 {
                    new_tape |= 1 << i;
                }
            }
            // cache the result
            // println!("cache miss, inserting");
            // println!("  steps: {}", steps);
            // println!("  old_tape: {:016b}", n);
            // println!("  new_tape: {:08b}", new_tape);
            // println!("  offset: {}", (self.position - start_pos) as i8);
            // println!("  new_state: {}", self.state);
            self.cache[start_state as usize].insert(
                n,
                CacheResult {
                    steps: steps,
                    new_tape: new_tape,
                    offset: (self.position - start_pos) as i8,
                    new_state: self.state,
                },
            );
            return Progress::DidNotHalt(steps);
        }
    }

    #[allow(unused)]
    fn draw(&mut self, img: &mut ImageBuffer<Luma<u8>, Vec<u8>>) -> bool {
        if self.row >= img.height() {
            return false;
        }

        let halfw = img.width() / 2;
        for x in 0..halfw - 1 {
            let bin = x / SIZE as u32;
            let bit = x % SIZE as u32;
            if bin >= self.tape.right.len() as u32 {
                break;
            }
            img.put_pixel(
                x + halfw,
                self.row,
                Luma([(((self.tape.right[bin as usize] >> bit) as u8) & 1) * 255]),
            );
        }

        for x in 0..halfw - 2 {
            let bin = x / SIZE as u32;
            let bit = x % SIZE as u32;
            if bin >= self.tape.left.len() as u32 {
                break;
            }
            img.put_pixel(
                halfw - x - 1,
                self.row,
                Luma([(((self.tape.left[bin as usize] >> bit) as u8) & 1) * 255]),
            )
        }

        self.row += 1;

        true
    }

    #[allow(unused)]
    fn draw_centered(&mut self, img: &mut ImageBuffer<Luma<u8>, Vec<u8>>) -> bool {
        if self.row >= img.height() {
            return false;
        };

        let halfw = img.width() as i32 / 2;
        let pos = self.position;

        for x in 0..img.width() {
            // let col = self.tape.get(x as i64 - halfw) * 255;
            let col = self.tape.get(pos - halfw + x as i32) * 255;
            img.put_pixel(x, self.row, Luma([col]));
        }

        self.row += 1;

        true
    }
}

use rayon::prelude::*;

fn main() {
    for i in -32..=32i16 {
        // println!("{}: {}", i, i.rem_euclid(16));
        println!("  {}: {:016b}", i, 0b0000000100000000u16 << i)
    }

    let args = env::args().collect::<Vec<_>>();
    let holdouts = read_to_string(args.get(1).expect("missing holdouts file"))
        .expect("failed to read holdouts file");
    let holdouts = holdouts.lines().collect::<Vec<&str>>();

    let max_steps = args
        .get(2)
        .expect("missing max steps")
        .parse::<u64>()
        .expect("max_steps must be an integer");

    // run_parallel(holdouts, max_steps);
    for holdout in holdouts.iter().take(10) {
        run_cached(holdout, max_steps);
    }
    // draw(holdouts[0], max_steps);
}

#[allow(unused)]
fn draw(holdout: &str, max_steps: u64) {
    println!("holdout: {:?}", holdout);
    let mut b = Beaver::new(&holdout.to_string());
    let mut img = image::GrayImage::new(1024, 1024);
    let max_steps = 10000000u64;
    let mut count = 0u64;
    loop {
        count += 1;

        if !b.draw(&mut img) {
            break;
        }

        if !b.step() || count >= max_steps {
            println!("halted at step {}", count);
            break;
        }
    }
    img.save("out.png").unwrap();
}

fn run_cached(holdout: &str, max_steps: u64) {
    let mut b = Beaver::new(&holdout.to_string());
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

    if halted {
        println!("halted at step {}", steps);
    } else {
        println!("did not halt after {} steps", steps);
    }
    println!(
        "simulated {} steps in {:.2}s ({} steps/sec)",
        steps,
        t0.elapsed().as_secs_f32(),
        steps as f32 / t0.elapsed().as_secs_f32()
    );
}

#[allow(unused)]
fn run_parallel(holdouts: Vec<&str>, max_steps: u64) {
    let halt_count = AtomicU32::new(0);
    let processed_count = AtomicU32::new(0);

    holdouts.par_iter().for_each(|holdout| {
        // println!("holdout: {:?}", holdout);

        let mut b = Beaver::new(&holdout.to_string());
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
    let mut b = Beaver::new(&"1RB1RG_1RC0RB_1RD0RE_1RE0LD_1LF1RA_---1LD_1LC1RG".to_string());

    println!("{:?}", b.transitions);

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

    println!("position: {}", b.position);
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
fn rle_test() {
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_get_neighborhood2() {
        let mut b = Beaver::new(&"1RB1LE_0LC1LG_0LD1RD_1RC1LA_1LF0RE_0LD1LB_---0LA".to_string());
        for _ in 0..1000 {
            b.step();
        }
        println!("{}", b.tape);
        println!("position: {}", b.position);
        let n = b.tape.get_neighborhood2(b.position);
        println!("neighborhood1: {:016b}", n);
        b.position = 7;
        let n = b.tape.get_neighborhood2(b.position);
        println!("neighborhood2: {:016b}", n);
        assert!(false);
    }
}
