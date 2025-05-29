use std::collections::HashMap;

use image::{ImageBuffer, Luma};
use itertools::{self, Itertools};

use crate::{Tape, types::*};

use strum::IntoEnumIterator;
use strum_macros::EnumIter;

#[derive(EnumIter, Clone, Copy)]
enum Write {
    Zero,
    One,
}

#[derive(EnumIter, Clone, Copy)]
enum Action {
    Read,
    Write,
}

#[derive(EnumIter, Clone, Copy)]
enum Dir {
    Left,
    Right,
}

/// Transitions are stored as a u8, where
/// - bit 0 is the write action (0 for zero, 1 for one)
/// - bit 1 is the direction (0 for left, 1 for right)
/// - bits 2-7 are the next state (0-63)
#[derive(Clone, Copy, Debug)]
struct Transition(u8);
impl Transition {
    fn new(transition: (Write, Dir, u8)) -> Self {
        let (write, dir, next_state) = transition;
        if next_state >= 0b111111 {
            return Transition(0xff);
        }
        let write = match write {
            Write::Zero => 0,
            Write::One => 1,
        };
        let step = match dir {
            Dir::Left => 0,
            Dir::Right => 1,
        };
        Transition(write | (step << 1) | (next_state << 2))
    }

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

    fn to_string(&self) -> String {
        if self.0 == 0xff {
            return "---".to_string();
        }
        let write = if self.write() { '1' } else { '0' };
        let dir = if self.dir() == -1 { 'L' } else { 'R' };
        let next_state = (self.next_state() + b'A') as char;
        format!("{}{}{}", write, dir, next_state)
    }

    fn write(self) -> bool {
        self.0 & 1 > 0
    }

    fn dir(self) -> i8 {
        (self.0 & 2) as i8 - 1
    }

    fn next_state(self) -> u8 {
        (self.0 & 0b01111100) >> 2
    }

    fn valid(self) -> bool {
        self.0 != 0xff
    }
}

struct CacheResult {
    steps: u16,
    new_tape: CacheBucket,
    offset: i8,
    new_state: u8,
}

pub struct Beaver {
    state: u8,
    cache: Vec<HashMap<CacheBucket, CacheResult>>,
    position: i32,
    transitions: Vec<[Transition; 2]>,
    pub tape: Tape,
    row: u32,
    #[cfg(feature = "track_states")]
    states_visitied: u64,
}

pub enum Progress {
    Halted(u16),
    DidNotHalt(u16),
}

const MAX_CACHED_STEPS: u16 = u16::MAX;

impl Beaver {
    pub fn new(transition_list: Vec<&Transition>) -> Self {
        let mut transitions = Vec::new();
        for chunk in transition_list.chunks_exact(2) {
            transitions.push([chunk[0].clone(), chunk[1].clone()]);
        }

        Self {
            cache: Vec::new(),
            state: 0,
            position: 0,
            transitions,
            tape: Tape::new(),
            row: 0,
            #[cfg(feature = "track_states")]
            states_visitied: 1, // initial state is visited
        }
    }

    pub fn from_string(rule: &String) -> Self {
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

    pub fn get_tape_strings(&self) -> (String, String) {
        self.tape.get_strings(self.position)
    }

    pub fn get_rule(&self) -> String {
        self.transitions
            .iter()
            .map(|[l, r]| {
                let mut s = r.to_string();
                s.push_str(l.to_string().as_str());
                s
            })
            .join("_")
    }

    pub fn step(&mut self) -> bool {
        let tape_state = self.tape.get(self.position);
        let transition = self.transitions[self.state as usize][tape_state as usize];

        if transition.0 == 0xff {
            return false;
        } else {
            self.state = transition.next_state();
            #[cfg(feature = "track_states")]
            {
                self.states_visitied |= 1 << self.state;
            }
        }

        self.tape.set(self.position, transition.write());
        self.position += transition.dir() as i32;

        true
    }

    pub fn step_with_cache(&mut self) -> Progress {
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
    pub fn draw(&mut self, img: &mut ImageBuffer<Luma<u8>, Vec<u8>>) -> bool {
        if self.row >= img.height() {
            return false;
        }

        let halfw = img.width() as i32 / 2;
        for x in 0..img.width() {
            let bit = self.tape.get(x as i32 - halfw);
            img.put_pixel(x, self.row, Luma([bit * 255]));
        }

        self.row += 1;

        true
    }

    pub fn thumbnail(&mut self, w: u32, h: u32, skip: u32) -> ImageBuffer<Luma<u8>, Vec<u8>> {
        let mut img = ImageBuffer::new(w, h);
        self.row = 0;

        let mut jump = skip;

        while self.draw(&mut img) {
            for _ in 0..jump {
                if !self.step() {
                    break;
                }
            }
            // jump += skip;
        }

        img
    }

    pub fn save_thumbnail(&mut self, w: u32, h: u32, skip: u32, fname: String) {
        self.thumbnail(w, h, skip)
            .save(fname)
            .expect("Failed to save image");
    }

    #[allow(unused)]
    pub fn draw_centered(&mut self, img: &mut ImageBuffer<Luma<u8>, Vec<u8>>) -> bool {
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

    pub fn enumerate(num_states: u8) -> Vec<Beaver> {
        assert!(num_states < 63, "Number of states must be less than 63");
        let mut transitions = itertools::iproduct!(Write::iter(), Dir::iter(), 0..num_states,)
            .map(Transition::new)
            .collect::<Vec<_>>();
        transitions.push(Transition(0xff)); // add invalid transition

        let permutations = transitions.iter().permutations(num_states as usize * 2);
        // .combinations_with_replacement(num_states as usize * 2);

        let permutations = permutations.map(Beaver::new).collect::<Vec<_>>();

        permutations
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_get_neighborhood2() {
        let mut b =
            Beaver::from_string(&"1RB1LE_0LC1LG_0LD1RD_1RC1LA_1LF0RE_0LD1LB_---0LA".to_string());
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
