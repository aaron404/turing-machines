// use crate::types::*;

pub fn rle(s: &[u8]) -> Vec<(usize, &[u8])> {
    let length = s.len();
    let mut rle_result = vec![];

    let mut best_best = 0;
    let mut i = 0;
    loop {
        // score is calculate based on the total number of bytes that would be
        // consumed by a chunk of this size
        let mut best_score = 0;
        let mut best_chunk_size = 1;
        let mut best_chunk_count = 1;
        let mut best_chunk = &s[i..i + 1];
        for chunk_size in 1..(length - i) {
            let mut iter = s[i..].chunks_exact(chunk_size);
            let chunk = iter.next().unwrap();
            let copies = iter.take_while(|&c| c == chunk).count();
            let score = chunk_size * (copies);
            if score > best_score {
                best_score = score;
                best_chunk_size = chunk_size;
                best_chunk_count = copies + 1;
                best_chunk = chunk;
                if chunk_size > best_best {
                    best_best = chunk_size;
                    // println!(
                    //     "best chunk: {} ({})",
                    //     String::from_utf8_lossy(chunk),
                    //     chunk_size
                    // );
                }
            }
        }

        rle_result.push((best_chunk_count, best_chunk));
        let advance = if best_score == 0 {
            1
        } else {
            best_chunk_size * (1 + best_score / best_chunk_size)
        };
        i += advance;

        if i >= length {
            break;
        }
    }

    rle_result
}

// start, len
// pub struct Span(usize, usize);

// pub fn rle_bits(s: &[Bin], len: usize) -> Vec<(usize, Span)> {
//     let mut rle_result = vec![];

//     let mut tmp = Vec::<Bin>::with_capacity(len);

//     let mut best_best = 0;
//     let mut start = 0;
//     loop {
//         // score is calculate based on the total number of bytes that would be
//         // consumed by a chunk of this size
//         let mut best_score = 0;
//         let mut best_chunk_size = 1;
//         let mut best_chunk_count = 1;
//         let mut best_chunk = (0, 1);

//         println!("start: {start}");
//         for chunk_size in 1..(len - start) {
//             // clear the tmp buffer
//             let num_bins = chunk_size / SIZE;
//             tmp.iter_mut().take(num_bins).for_each(|x| *x = 0);

//             let mut iter = (start..len).step_by(chunk_size);
//             for i in iter {
//                 println!("i: {i}");
//                 // get first chunk to make comparisons against, store in tmp
//             }

//             // let mut iter = s[start..].chunks_exact(chunk_size);
//             // let chunk = iter.next().unwrap();
//             // let copies = iter.take_while(|&c| c == chunk).count();
//             // let score = chunk_size * (copies);
//             // if score > best_score {
//             //     best_score = score;
//             //     best_chunk_size = chunk_size;
//             //     best_chunk_count = copies + 1;
//             //     best_chunk = chunk;
//             //     if chunk_size > best_best {
//             //         best_best = chunk_size;
//             //         println!(
//             //             "best chunk: {} ({})",
//             //             String::from_utf8_lossy(chunk),
//             //             chunk_size
//             //         );
//             //     }
//             // }
//         }

//         rle_result.push((best_chunk_count, best_chunk));
//         let advance = if best_score == 0 {
//             1
//         } else {
//             best_chunk_size * (1 + best_score / best_chunk_size)
//         };
//         start += advance;

//         if start >= len {
//             break;
//         }
//     }

//     todo!()
// }

pub fn print_rle(rle: Vec<(usize, &[u8])>) {
    print!("0∞ ");

    for (count, chunk) in rle {
        let chunk_str = String::from_utf8_lossy(chunk);
        print!("{chunk_str}{} ", superscript(count));
    }

    println!("0∞");
}

fn superscript(n: usize) -> String {
    let mut result = String::new();
    for c in n.to_string().chars() {
        let c = match c {
            '0' => '⁰',
            '1' => '¹',
            '2' => '²',
            '3' => '³',
            '4' => '⁴',
            '5' => '⁵',
            '6' => '⁶',
            '7' => '⁷',
            '8' => '⁸',
            '9' => '⁹',
            _ => c,
        };
        result.push(c);
    }
    result
}
