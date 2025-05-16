use std::{
    collections::BTreeMap, 
    fs::create_dir,
    fs::File, 
    io::{BufRead, BufReader}, 
    path::Path, 
    string::String, 
    sync::{Arc, Mutex}, thread
};
use clap::Parser;
use std::sync::mpsc::Sender;
use memmap2::MmapOptions;
use flate2::read::GzDecoder;
use std::collections::HashMap;
// use needletail::parse_fastx_file;
use log::{debug,warn,info,error};
use env_logger;
use std::time::Instant;
use tokio::io::AsyncWriteExt;

mod utils;
use utils::*;

#[derive(Parser, Debug)]
#[command(
    author = "BENM",
    version = "version 0.1.15",
    about = "SegPE: a simple program for classifing and separating PE FASTQ after removing adapters via adapter-index seq",
    long_about = "Algorithmic ideas of SegPE:
    - Exact match search: This is a more direct method for finding exact matches of artificial sequences.
    Regular expression matching and Hamming distance:
    This is suitable for detecting and locating index sequences, especially when mismatches of a certain length are taken into account.

    - Process and classify PE sequences: After removing the adapter and index sequences, classify the PE sequences and create new PE FASTQ files.

    - Use Multi-threads, SIMD and AsyncIO to handle large amounts of data.
    "
)]

struct Args {
    /// Path of 5' artificial fasta file
    #[arg(long)]
    five_art_fa: String,

    /// Path of 3' artificial fasta file
    #[arg(long)]
    three_art_fa: String,

    /// Path of 5‘ index fasta file
    #[arg(long)]
    five_idx_fa: String,

    /// Path of 3’ index fasta file
    #[arg(long, default_value_t = String::from(""))]
    three_idx_fa: String,

    /// Location of index,
    /// 1: PE1
    /// 2: PE2
    /// 3: both
    #[arg(long, default_value_t = 1)]
    idx_loc: u8,

    /// Path of PE1 fastq file
    #[arg(long)]
    pe1_fastq: String,

    /// Path of PE1 fastq file
    #[arg(long, default_value_t = String::from(""))]
    pe2_fastq: String,

    /// Number of seed length, not allow to longer than index length
    #[arg(short, long, default_value_t = 6)]
    seed_len: usize,

    /// Merged overlapped PE
    #[arg(long, default_value_t = false)]
    merge_pe: bool,

    /// Number of error tolerance
    #[arg(long, default_value_t = 1)]
    error_tolerance: u8,

    /// Number of alignment macth score
    #[arg(short, long, default_value_t = 1)]
    match_score: i8,

    /// Number of alignment mismacth score
    #[arg(long, default_value_t = -1)]
    error_score: i8,

    /// Number of alignment gap open score
    #[arg(long, default_value_t = -5)]
    gap_open_score: i8,

    /// Number of alignment gap extend score
    #[arg(long, default_value_t = -1)]
    gap_extend_score: i8,

     /// Number of cucurrency threads
     #[arg(short, long, default_value_t = 8)]
    num_threads: usize,

    /// batch size of reads, which every thread need to handle
    #[arg(short, long, default_value_t = 10_000)]
    batch: u32,

    /// pretrain size of reads
    #[arg(long, default_value_t = 0)]
    train: usize,

    /// add trim info in reads_name
    #[arg(long, default_value_t = false)]
    trim_name: bool,

    /// Path of output directory
    #[arg(short, long, default_value_t = String::from("./"))]
    outdir: String,
}

#[derive(Clone)]
struct PEreads {
    read1_name: String,
    read1_vec: Vec<String>,
    read2_name: String,
    read2_vec: Vec<String>
}

// mrege PE
fn merge_pe(pair: (Vec<String>, Vec<String>), seed_length: usize,
    error_tolerance: u8,match_score: i8, mismatch_score: i8,
    gap_open_score: i8, gap_extend_score: i8) -> (Vec<String>, Vec<String>) {
    let reads1 = pair.0.clone();
    let reads2 = pair.1.clone();
    let mut new_read1: Vec<String> = Vec::new();
    let mut new_read2: Vec<String> = Vec::new();
    for i in 0..reads2.len()/4 {
        new_read1.push(reads1[i*4+0].clone());
        new_read1.push(reads1[i*4+1].clone());
        new_read1.push(reads1[i*4+2].clone());
        new_read1.push(reads1[i*4+3].clone());
        new_read2.push(reads2[i*4+0].clone());
        new_read2.push(reads2[i*4+1].clone());
        new_read2.push(reads2[i*4+2].clone());
        new_read2.push(reads2[i*4+3].clone());

        let read1 = &reads1[i*4+1];
        // 使用优化的反向互补
        let read2_rc = fast_reverse_complement(&reads2[i*4+1]);
        let qual2_r: String = reverse(&reads2[i*4+3]);
        // 从read2_rc 5'端第一个碱基开始到第二个碱基（只允许一个错配）取seed_length长度的序列在read1中查找是否有匹配的位置
        let mut score = 0;
        for j in 0..2 {
            let seed = &read2_rc[j..j+seed_length];
            // 使用优化的序列搜索
            let match_pos = fast_sequence_search(&read1, seed);
            // 如果找到了匹配的位置，则将read1中的从第pos-j碱基到3‘端的序列，标记长度为overlapping_len，与read2_rc从开始到overlapping_len序列进行Smith-Waterman比对,
            // 如果比对的错误率小于阈值，则将read1中的第pos-j碱基与read2_rc中的第0碱基对齐合并
            if match_pos.len() > 0 {
                for pos in match_pos {
                    if pos != -1i16 && pos >= j as i16 {
                        let pos_usize = pos as usize;
                        let overlapping_len = read1.len() - pos_usize + j;
                        // score = sw_align(&read1[pos_usize-j..], &read2_rc[..overlapping_len], match_score, mismatch_score, gap_open_score, gap_extend_score);
                        if overlapping_len<read2_rc.len(){
                            let res = block_align(&read1[pos_usize-j..], &read2_rc[..overlapping_len], match_score, mismatch_score, gap_open_score, gap_extend_score);
                            if res.score > overlapping_len as i32 - error_tolerance as i32 {
                                score = res.score;
                                let (seq1_new, qual1_new) = merge_pe_by_qual(&read1, &read2_rc, &reads1[i*4+3], &qual2_r, pos_usize-j);
                                new_read1[i*4+0]=new_read1[i*4+0].trim_end().to_owned() + " merged";
                                new_read1[i*4+1] = seq1_new;
                                new_read1[i*4+3] = qual1_new;
                                new_read2[i*4+0] = "".to_string();
                                new_read2[i*4+1] = "".to_string();
                                new_read2[i*4+2] = "".to_string();
                                new_read2[i*4+3] = "".to_string();
                                break;
                            }
                        }

                    }
                }
                if score>0{
                    break;
                }
            }
        }
    }

    (new_read1, new_read2)
}

fn fastmap_five_art_idx(five_artificial_map: &BTreeMap<Vec<u8>, String>, five_index_map: &BTreeMap<Vec<u8>, String>) -> (usize, BTreeMap<Vec<u8>, String>) {
    let mut five_artificial_idx_map: BTreeMap<Vec<u8>, String> = BTreeMap::new();
    let mut stat_artfical_length: HashMap<usize, i32> = HashMap::new();
    for fa_map in five_artificial_map.iter() {
        let art_seq = String::from_utf8(fa_map.0.clone()).unwrap();
        let count = stat_artfical_length.entry(art_seq.len()).or_insert(1);
        *count += 1;
        for idx_map in five_index_map.iter() {
            let idx_seq = String::from_utf8(idx_map.0.clone()).unwrap();
            let match_pos = regrex_sarch(&art_seq, &idx_seq);
            if match_pos.len() > 0 {
                five_artificial_idx_map.insert(fa_map.0.clone(), idx_map.1.clone());
                break;
            }
            // let score = sparse_align(&art_seq, &idx_seq, idx_seq.len());
            // if score>=idx_seq.len() as i32 - 1 {
            //     five_artificial_idx_map.insert(fa_map.0.clone(), idx_map.1.clone());
            //     break;
            // }
        }
        if !five_artificial_idx_map.contains_key(fa_map.0){
            five_artificial_idx_map.insert(fa_map.0.clone(), String::from("unknown"));
        }
    }
    let mut max_count = 0;
    let mut max_length = 0;
    for (length, count) in stat_artfical_length.iter() {
        if *count > max_count {
            max_count = *count;
            max_length = *length;
        }
    }
    (max_length, five_artificial_idx_map)

}

fn fastmap_art_seed(artificial_map: &BTreeMap<Vec<u8>, String>, seed_len: usize) -> BTreeMap<String, Vec<(String,u8)>> {
    let mut art_seed_map: BTreeMap<String, Vec<(String,u8)>> = BTreeMap::new();
    for fa_map in artificial_map.iter() {
        let art_seq = String::from_utf8(fa_map.0.clone()).unwrap();
        let id = fa_map.1.clone();
        for i in 0..art_seq.len()-seed_len {
            let seed = &art_seq[i..i+seed_len];
            if art_seed_map.contains_key(seed) {
                art_seed_map.get_mut(seed).unwrap().push((id.clone(),i as u8));
            } else {
                art_seed_map.insert(seed.to_string(), vec![(id.clone(),i as u8)]);
            }
        }
    }
    art_seed_map

}

fn classify_reads(five_artificial_idx_map: &BTreeMap<Vec<u8>, String>, three_artificial_map: &BTreeMap<Vec<u8>, String>,
    five_index_map: &BTreeMap<Vec<u8>, String>, three_index_map: &BTreeMap<Vec<u8>, String>, 
    five_artificial_id_map: &BTreeMap<String, Vec<u8>>, three_artificial_id_map: &BTreeMap<String, Vec<u8>>,
    five_art_seed_map: &BTreeMap<String, Vec<(String,u8)>>, three_art_seed_map: &BTreeMap<String, Vec<(String,u8)>>,
    five_art_length: usize, index_map_bool: bool,
    reads: &mut Vec<String>, seed_length: usize,
    error_tolerance: u8, match_score: i8, mismatch_score: i8,
    gap_open_score: i8, gap_extend_score: i8, trim_name: bool) -> (Vec<(String, String)>,BTreeMap<Vec<u8>, String>) {
    // mutable borrow from five_artificial_idx_map
    let mut five_artificial_idx_map_ext: BTreeMap<Vec<u8>, String> = five_artificial_idx_map.clone();

    let mut cls = Vec::new();
    // Implement the exact matching, regular expression matching, and Hamming distance calculation here
    for i in 0..reads.len()/4 {  // 遍历每一个reads
        let mut cls5 = String::from("");
        let mut cls3 = String::from("");
        let read_len = reads[i*4+1].len();
        let mut left_trim_pos = 0;
        let mut right_trim_pos = read_len;
        let mut cls5_bool = false;

        // Check if read is long enough for five_art_length processing
        if read_len > five_art_length {
            let trimed_seq_bytes = (&reads[i*4+1][0..five_art_length]).as_bytes();
            if five_artificial_idx_map_ext.contains_key(trimed_seq_bytes) {
                left_trim_pos = five_art_length;
                cls5 = five_artificial_idx_map_ext[trimed_seq_bytes].to_owned();
                cls5_bool = true;
            }
        }

        let mut break_bool = false;
        if !cls5_bool {
            if read_len > seed_length {
                //filter artificial seq from 5'-read

                // let mut map_id5 = String::new();
                let mut max_score5 = 0;
                
                if five_art_seed_map.len() > 0 {
                    for j in 0..1+error_tolerance as usize { // 从5'端开始，每次向右移动一个碱基至error_tolerance，如果error_tolerance=1, 即0,1
                        if j+seed_length < read_len{
                            let seed = &reads[i*4+1][j..j+seed_length];
                            if five_art_seed_map.contains_key(seed) {
                                let vec_map = five_art_seed_map.get(seed); 
                                for map in vec_map.unwrap() {
                                    let map_id5 = map.0.clone(); 
                                    if five_artificial_id_map.contains_key(&map_id5) {
                                        let art_seq = String::from_utf8(five_artificial_id_map.get(&map_id5).unwrap().to_owned()).unwrap();
                                        let pos = map.1 as usize;
                                        let overlapping_len = art_seq.len() - pos + j;                                                                                    
                                        if overlapping_len < read_len && overlapping_len >= seed_length {
                                            let res: block_aligner::scan_block::AlignResult = block_align(&art_seq, &reads[i * 4 + 1][..overlapping_len], match_score, mismatch_score, gap_open_score, gap_extend_score); //TODO
                                            if res.score > max_score5 && res.score >= overlapping_len as i32 - error_tolerance as i32 {
                                                max_score5 = res.score;
                                                cls5 = map_id5.clone();
                                                cls5_bool = true;
                                                left_trim_pos = overlapping_len as usize;
                                                if max_score5 == overlapping_len as i32 {
                                                    break_bool = true;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                                if break_bool {
                                    break;
                                }
                            }
                        }
                        if break_bool {
                            break;
                        }
                    }
                }
                if !cls5_bool && index_map_bool && max_score5 > 0  && five_index_map.len()>0{
                    let map_artseq = &reads[i*4+1][..left_trim_pos]; //String::from_utf8(artificial_map[&map_id5].clone()).unwrap();
                    if map_artseq.len() >= seed_length {
                        for idx_map in five_index_map {
                            // let idx_id = idx_map.0;
                            let idx_id = idx_map.1;
                            // let idx_seq = String::from_utf8(idx_map.1.clone()).unwrap();
                            let idx_seq = String::from_utf8(idx_map.0.clone()).unwrap();
                            let match_pos = regrex_sarch(&map_artseq, &idx_seq);
                            if match_pos.len() > 0 {
                                cls5 = idx_id.to_owned();
                                break;
                            }
                            let score = sparse_align(&map_artseq, &idx_seq, idx_seq.len());
                            if score>=idx_seq.len() as i32 - error_tolerance as i32 {
                                // cls.push(idx_id.to_string());
                                cls5 = idx_id.to_owned();
                                break;
                            }
                        }

                    }
                }
            }
        }
        //filter artificial seq from 3'-read
        let mut max_score3 = 0;
        break_bool = false;
        if three_artificial_map.len() > 0 && read_len > seed_length {
            for j in 0..1+error_tolerance as usize {
                // Check if we have enough sequence for the seed
                if read_len > seed_length + j {
                    let seed = &reads[i*4+1][read_len-seed_length-j..read_len-j];
                    if three_art_seed_map.contains_key(seed) {
                        let vec_map = three_art_seed_map.get(seed); 
                        for map in vec_map.unwrap() {
                            let map_id3 = map.0.clone(); 
                            // Fix this check - we should look in three_artificial_id_map not three_art_seed_map
                            if three_artificial_id_map.contains_key(&map_id3) {
                                let art_seq = String::from_utf8(three_artificial_id_map.get(&map_id3).unwrap().to_owned()).unwrap();
                                let pos = map.1 as usize;
                                let overlapping_len = pos+seed_length+j;                                                              
                                
                                // Additional safety checks for indexing
                                if overlapping_len >= seed_length && 
                                   pos + overlapping_len < art_seq.len() && 
                                   read_len > seed_length + j + pos && 
                                   read_len - seed_length - j - pos < reads[i*4+1].len() {
                                    
                                    let res = block_align(&art_seq[0..pos+overlapping_len], 
                                                         &reads[i*4+1][read_len-seed_length-j-pos..], 
                                                         match_score, mismatch_score, gap_open_score, gap_extend_score);
                                    
                                    if res.score > max_score3 && res.score >= overlapping_len as i32 - error_tolerance as i32 {
                                        max_score3 = res.score;
                                        right_trim_pos = read_len-seed_length-j-pos;
                                        if max_score3 == overlapping_len as i32 {
                                            break_bool = true;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        if break_bool {
                            break;
                        }
                    }
                }
                if break_bool {
                    break;
                }
            }
        }
        if index_map_bool && max_score3 > 0 && three_index_map.len()>0{
            let map_artseq = &reads[i*4+1][right_trim_pos..]; //String::from_utf8(artificial_map[&map_id3].clone()).unwrap();
            if map_artseq.len() >= seed_length {
                for idx_map in three_index_map {
                    let idx_id = idx_map.1;
                    let idx_seq = String::from_utf8(idx_map.0.clone()).unwrap();
                    let match_pos = regrex_sarch(&map_artseq, &idx_seq);
                    if match_pos.len() > 0 {
                        cls3 = idx_id.to_owned();
                        break;
                    }
                    let score = sparse_align(&map_artseq, &idx_seq, idx_seq.len());
                    if score>=idx_seq.len() as i32 - error_tolerance as i32 {
                        cls3 = idx_id.to_owned();
                        break;
                    }
                }
            }
        }
            // if cls5 != "" || cls3 != "" {
        if trim_name {
            if left_trim_pos > 0 && right_trim_pos < read_len {
                reads[i*4+0] = format!("{} |trimed:0..{},{}..end", &reads[i*4+0].trim_end(), left_trim_pos, right_trim_pos);
            } else if left_trim_pos > 0 {
                reads[i*4+0] = format!("{} |trimed:0..{}", &reads[i*4+0].trim_end(), left_trim_pos);
            } else if right_trim_pos < read_len {
                reads[i*4+0] = format!("{} |trimed:{}..end", &reads[i*4+0].trim_end(), right_trim_pos);
            }
        }
                
        if cls5 != "" {
            let left_trim_seq = (&reads[i * 4 + 1][0..five_art_length]).as_bytes();
            if !five_artificial_idx_map.contains_key(left_trim_seq) {
                five_artificial_idx_map_ext.insert(left_trim_seq.to_owned(), cls5.to_string());
            }
        }
        cls.push((cls5.to_string(),cls3.to_string()));

        //trim artificial seq
        reads[i*4+1] = reads[i*4+1][left_trim_pos..right_trim_pos].to_string();
        reads[i*4+3] = reads[i*4+3][left_trim_pos..right_trim_pos].to_string();
    }
    (cls, five_artificial_idx_map_ext)
}

// 定义处理函数
#[allow(dead_code)]
fn process_reads(five_artificial_idx_map: &BTreeMap<Vec<u8>, String>, three_artificial_map: &BTreeMap<Vec<u8>, String>,
        five_index_map: &BTreeMap<Vec<u8>, String>, three_idx_map: &BTreeMap<Vec<u8>, String>, 
        five_artificial_id_map: &BTreeMap<String, Vec<u8>>, three_artificial_id_map: &BTreeMap<String, Vec<u8>>,
        five_art_seed_map: &BTreeMap<String, Vec<(String,u8)>>, three_art_seed_map: &BTreeMap<String, Vec<(String,u8)>>,
        five_art_length: usize, idx_loc: u8,
        prefix1: &str, prefix2: &str, pair: (Vec<String>, Vec<String>), seed_length: usize, merge_pe_bool: bool,
        error_tolerance: u8, match_score: i8, mismatch_score: i8,
        gap_open_score: i8, gap_extend_score: i8, io_tx: Sender<Vec<PEreads>>, trim_name: bool) {

    let start = Instant::now();

    let mut reads1: Vec<String> = pair.0.clone();
    let mut reads2: Vec<String> = pair.1.clone();
    
    // Check if we have empty reads2 (single-end mode)
    let single_end_mode = reads2.is_empty();
    
    if merge_pe_bool && !single_end_mode {
        let new_pair = merge_pe(pair, seed_length, error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score);
        reads1 = new_pair.0.clone();
        reads2 = new_pair.1.clone();
    } 
    let (cls1, _) = classify_reads(five_artificial_idx_map, three_artificial_map,
        five_index_map, three_idx_map, five_artificial_id_map, three_artificial_id_map, five_art_seed_map, three_art_seed_map,
        five_art_length, idx_loc & 0b1 == 0b1, &mut reads1, seed_length,
        error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score, trim_name);
    
    // Only process reads2 if it's not empty
    let cls2 = if !single_end_mode {
        let (cls, _) = classify_reads(five_artificial_idx_map, three_artificial_map,
            five_index_map, three_idx_map, five_artificial_id_map, three_artificial_id_map, five_art_seed_map, three_art_seed_map, 
            five_art_length, idx_loc & 0b10 == 0b10, &mut reads2, seed_length,
            error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score, trim_name);
        cls
    } else {
        // Create an empty classification vector of the same length as cls1
        vec![(String::new(), String::new()); cls1.len()]
    };

    for i in 0..cls1.len() {
        let mut cls = "";

        // Check if reads2 has content at this position (could be empty for merged reads or single-end)
        let has_read2 = !single_end_mode && i*4+1 < reads2.len() && reads2[i*4+1] != "";
        
        if !has_read2 {
            if idx_loc & 0b1 == 0b1 && cls1[i].0 != "" {
                cls = &cls1[i].0;
            } else if idx_loc & 0b1 == 0b1 && cls1[i].1 != "" {
                cls = &cls1[i].1;
            }
        } else {
            if cls1[i].0 !="" && cls2[i].0 != ""  && (idx_loc & 0b01 == 0b01 || idx_loc & 0b10 == 0b10){
                if idx_loc & 0b01 == 0b01 {
                    cls = &cls1[i].0;
                } else if idx_loc & 0b10 == 0b10 {
                    cls = &cls2[i].0;
                }
            } else if cls1[i].0 !="" && cls2[i].0=="" && idx_loc & 0b01 == 0b01 {
                cls = &cls1[i].0;
            } else if cls1[i].0 =="" && cls2[i].0 !="" && idx_loc & 0b10 == 0b10 {
                cls = &cls2[i].0;
            }
        }
        
        if cls !="" {
            if has_read2 {
                let fq1_name = format!("{}_{}.fastq", prefix1.to_string(), cls);
                let fq2_name  = format!("{}_{}.fastq", prefix2.to_string(), cls);
                let fq1_data = reads1[i*4+0..i*4+4].to_vec();
                let fq2_data = reads2[i*4+0..i*4+4].to_vec();
                io_tx.send(vec![PEreads{read1_name: fq1_name, read1_vec: fq1_data, read2_name: fq2_name, read2_vec: fq2_data}])
                     .expect("Failed to send data through channel");
            } else {
                let suffix = cls.to_owned() + "_single";
                let fq1_name = format!("{}_{}.fastq", prefix1.to_string(), suffix);
                let fq1_data = reads1[i*4+0..i*4+4].to_vec();
                io_tx.send(vec![PEreads{read1_name: fq1_name, read1_vec: fq1_data, read2_name: "".to_string(), read2_vec: Vec::new()}])
                     .expect("Failed to send data through channel");
            }
        } else {
            if has_read2 {
                let suffix =  "unclassified";
                let fq1_name = format!("{}_{}.fastq", prefix1.to_string(), suffix);
                let fq2_name  = format!("{}_{}.fastq", prefix2.to_string(), suffix);
                let fq1_data = reads1[i*4+0..i*4+4].to_vec();
                let fq2_data = reads2[i*4+0..i*4+4].to_vec();
                io_tx.send(vec![PEreads{read1_name: fq1_name, read1_vec: fq1_data, read2_name: fq2_name, read2_vec: fq2_data}])
                     .expect("Failed to send data through channel");
            } else {
                let suffix =  "unclassified";
                let fq1_name = format!("{}_{}.fastq", prefix1.to_string(), suffix);
                let fq1_data = reads1[i*4+0..i*4+4].to_vec();
                io_tx.send(vec![PEreads{read1_name: fq1_name, read1_vec: fq1_data, read2_name: "".to_string(), read2_vec: Vec::new()}])
                     .expect("Failed to send data through channel");
            }
        }
    }

    let duration = start.elapsed();
    info!("Time elapsed in process_reads() is: {:?}", duration);
}

// 单端处理函数 - 专门处理单端FASTQ数据
fn process_reads_se(five_artificial_idx_map: &BTreeMap<Vec<u8>, String>, three_artificial_map: &BTreeMap<Vec<u8>, String>,
        five_index_map: &BTreeMap<Vec<u8>, String>, three_idx_map: &BTreeMap<Vec<u8>, String>, 
        five_artificial_id_map: &BTreeMap<String, Vec<u8>>, three_artificial_id_map: &BTreeMap<String, Vec<u8>>,
        five_art_seed_map: &BTreeMap<String, Vec<(String,u8)>>, three_art_seed_map: &BTreeMap<String, Vec<(String,u8)>>,
        five_art_length: usize, idx_loc: u8,
        prefix1: &str, reads1: Vec<String>, seed_length: usize,
        error_tolerance: u8, match_score: i8, mismatch_score: i8,
        gap_open_score: i8, gap_extend_score: i8, io_tx: Sender<Vec<PEreads>>, trim_name: bool) {
    
    let start = Instant::now();

    let mut reads1 = reads1;  // 允许修改reads1
    
    // 直接进行分类，无需检查merge_pe标志
    let (cls1, _) = classify_reads(five_artificial_idx_map, three_artificial_map,
        five_index_map, three_idx_map, five_artificial_id_map, three_artificial_id_map, five_art_seed_map, three_art_seed_map,
        five_art_length, idx_loc & 0b1 == 0b1, &mut reads1, seed_length,
        error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score, trim_name);
    
    // 预分配结果向量以减少重新分配
    let mut results = Vec::with_capacity(cls1.len());
    
    // 单端处理逻辑 - 简化版，无需检查read2是否存在
    for i in 0..cls1.len() {
        let mut cls = "";

        // 只需检查read1的分类结果
        if idx_loc & 0b1 == 0b1 && !cls1[i].0.is_empty() {
            cls = &cls1[i].0;
        } else if idx_loc & 0b1 == 0b1 && !cls1[i].1.is_empty() {
            cls = &cls1[i].1;
        }
        
        // 根据分类结果创建输出文件名
        let suffix = if cls.is_empty() { "unclassified".to_string() } else { format!("{}", cls) };
        let fq1_name = format!("{}_{}.fastq", prefix1, suffix);
        
        // 仅当索引i*4+3在reads1范围内时才继续
        if i*4+3 < reads1.len() {
            let fq1_data = reads1[i*4..i*4+4].to_vec();
            results.push(PEreads{
                read1_name: fq1_name,
                read1_vec: fq1_data,
                read2_name: String::new(),
                read2_vec: Vec::new()
            });
        }
    }
    
    // 批量发送结果以减少channel开销
    if !results.is_empty() {
        io_tx.send(results).expect("Failed to send data through channel");
    }

    let duration = start.elapsed();
    debug!("Time elapsed in process_reads_se() is: {:?}", duration);
}

// 双端处理函数 - 专门处理双端FASTQ数据
fn process_reads_pe(five_artificial_idx_map: &BTreeMap<Vec<u8>, String>, three_artificial_map: &BTreeMap<Vec<u8>, String>,
        five_index_map: &BTreeMap<Vec<u8>, String>, three_idx_map: &BTreeMap<Vec<u8>, String>, 
        five_artificial_id_map: &BTreeMap<String, Vec<u8>>, three_artificial_id_map: &BTreeMap<String, Vec<u8>>,
        five_art_seed_map: &BTreeMap<String, Vec<(String,u8)>>, three_art_seed_map: &BTreeMap<String, Vec<(String,u8)>>,
        five_art_length: usize, idx_loc: u8,
        prefix1: &str, prefix2: &str, pair: (Vec<String>, Vec<String>), seed_length: usize, merge_pe_bool: bool,
        error_tolerance: u8, match_score: i8, mismatch_score: i8,
        gap_open_score: i8, gap_extend_score: i8, io_tx: Sender<Vec<PEreads>>, trim_name: bool) {

    let start = Instant::now();

    // 避免不必要的clone，直接接受所有权
    let (mut reads1, mut reads2) = pair;
    
    // 当启用merge_pe时，执行PE合并
    if merge_pe_bool {
        // 只在需要合并时克隆数据
        let merged_result = merge_pe((reads1.clone(), reads2.clone()), seed_length, 
                              error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score);
        reads1 = merged_result.0;
        reads2 = merged_result.1;
    }
    
    // 预估结果容量，避免动态扩容
    let reads_count = reads1.len() / 4;
    let mut results = Vec::with_capacity(reads_count);

    // 处理reads1和reads2
    let (cls1, _) = classify_reads(five_artificial_idx_map, three_artificial_map,
        five_index_map, three_idx_map, five_artificial_id_map, three_artificial_id_map, five_art_seed_map, three_art_seed_map,
        five_art_length, idx_loc & 0b1 == 0b1, &mut reads1, seed_length,
        error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score, trim_name);
    
    let (cls2, _) = classify_reads(five_artificial_idx_map, three_artificial_map,
        five_index_map, three_idx_map, five_artificial_id_map, three_artificial_id_map, five_art_seed_map, three_art_seed_map, 
        five_art_length, idx_loc & 0b10 == 0b10, &mut reads2, seed_length,
        error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score, trim_name);

    // 处理分类结果并收集输出
    for i in 0..cls1.len() {
        let cls = get_classification(&cls1[i], &cls2[i], idx_loc, reads2[i*4+1].is_empty());
        
        if reads2[i*4+1].is_empty() {
            // 处理合并后的reads
            let suffix = if cls.is_empty() { "unclassified".to_string() } else { format!("{}_merged", cls) };
            let fq1_name = format!("{}_{}.fastq", prefix1, suffix);
            
            // 避免不必要的克隆，使用切片引用
            let fq1_data = reads1[i*4..i*4+4].to_vec();
            
            results.push(PEreads{
                read1_name: fq1_name,
                read1_vec: fq1_data,
                read2_name: String::new(),
                read2_vec: Vec::new()
            });
        } else {
            // 处理标准双端reads
            let suffix = if cls.is_empty() { "unclassified" } else { cls };
            let fq1_name = format!("{}_{}.fastq", prefix1, suffix);
            let fq2_name = format!("{}_{}.fastq", prefix2, suffix);
            
            // 直接使用切片
            let fq1_data = reads1[i*4..i*4+4].to_vec();
            let fq2_data = reads2[i*4..i*4+4].to_vec();
            
            results.push(PEreads{
                read1_name: fq1_name,
                read1_vec: fq1_data,
                read2_name: fq2_name,
                read2_vec: fq2_data
            });
        }
    }
    
    // 批量发送结果
    if !results.is_empty() {
        io_tx.send(results).expect("Failed to send data through channel");
    }

    let duration = start.elapsed();
    debug!("Time elapsed in process_reads_pe() is: {:?}", duration);
}

// 辅助函数，获取分类结果
#[inline]
fn get_classification<'a>(cls1: &'a (String, String), cls2: &'a (String, String), idx_loc: u8, is_merged: bool) -> &'a str {
    if is_merged {
        if idx_loc & 0b1 == 0b1 && !cls1.0.is_empty() {
            &cls1.0
        } else if idx_loc & 0b1 == 0b1 && !cls1.1.is_empty() {
            &cls1.1
        } else {
            ""
        }
    } else {
        if !cls1.0.is_empty() && !cls2.0.is_empty() && (idx_loc & 0b01 == 0b01 || idx_loc & 0b10 == 0b10) {
            if idx_loc & 0b01 == 0b01 {
                &cls1.0
            } else if idx_loc & 0b10 == 0b10 {
                &cls2.0
            } else {
                ""
            }
        } else if !cls1.0.is_empty() && cls2.0.is_empty() && idx_loc & 0b01 == 0b01 {
            &cls1.0
        } else if cls1.0.is_empty() && !cls2.0.is_empty() && idx_loc & 0b10 == 0b10 {
            &cls2.0
        } else {
            ""
        }
    }
}

fn cache_channel(io_rx: std::sync::mpsc::Receiver<Vec<PEreads>>, tx: std::sync::mpsc::Sender<BTreeMap<String, Vec<String>>>, batch_size: usize) {
    // 预估常见文件名，为它们预分配空间
    let mut file_buffers: HashMap<String, Vec<String>> = HashMap::with_capacity(64);
    let mut buffer_sizes: HashMap<String, usize> = HashMap::with_capacity(64);
    
    // 跟踪最常用的文件以优化缓存
    let mut file_usage_counter: HashMap<String, usize> = HashMap::with_capacity(64);
    
    // 使用预分配的缓冲区
    for pereads in io_rx.iter() {
        // If we receive an empty Vec, treat it as a termination signal
        if pereads.is_empty() {
            break;
        }
        for peread in pereads {
            // 处理PE1数据
            if !peread.read1_vec.is_empty() {
                let file_name = &peread.read1_name;
                
                // 更新使用计数
                *file_usage_counter.entry(file_name.clone()).or_insert(0) += 1;
                
                // 获取或创建缓冲区
                if !file_buffers.contains_key(file_name) {
                    // 为常用文件预分配更大的缓冲区
                    let initial_capacity = if file_usage_counter.get(file_name).unwrap_or(&0) > &5 {
                        batch_size * 8
                    } else {
                        batch_size * 4
                    };
                    file_buffers.insert(file_name.clone(), Vec::with_capacity(initial_capacity));
                    buffer_sizes.insert(file_name.clone(), 0);
                }
                
                // 追加数据到缓冲区
                let buffer = file_buffers.get_mut(file_name).unwrap();
                buffer.extend(peread.read1_vec.iter().cloned());
                *buffer_sizes.get_mut(file_name).unwrap() += peread.read1_vec.len() / 4;
                
                // 当缓冲区足够大时发送数据
                if buffer_sizes.get(file_name).unwrap_or(&0) >= &batch_size {
                    let mut temp_map = BTreeMap::new();
                    temp_map.insert(file_name.clone(), std::mem::replace(buffer, Vec::with_capacity(batch_size * 4)));
                    tx.send(temp_map).unwrap();
                    buffer_sizes.insert(file_name.clone(), 0);
                }
            }
            
            // 处理PE2数据 (类似逻辑)
            if !peread.read2_name.is_empty() && !peread.read2_vec.is_empty() {
                let file_name = &peread.read2_name;
                
                // 更新使用计数
                *file_usage_counter.entry(file_name.clone()).or_insert(0) += 1;
                
                // 获取或创建缓冲区
                if !file_buffers.contains_key(file_name) {
                    // 为常用文件预分配更大的缓冲区
                    let initial_capacity = if file_usage_counter.get(file_name).unwrap_or(&0) > &5 {
                        batch_size * 8
                    } else {
                        batch_size * 4
                    };
                    file_buffers.insert(file_name.clone(), Vec::with_capacity(initial_capacity));
                    buffer_sizes.insert(file_name.clone(), 0);
                }
                
                // 追加数据到缓冲区
                let buffer = file_buffers.get_mut(file_name).unwrap();
                buffer.extend(peread.read2_vec.iter().cloned());
                *buffer_sizes.get_mut(file_name).unwrap() += peread.read2_vec.len() / 4;
                
                // 当缓冲区足够大时发送数据
                if buffer_sizes.get(file_name).unwrap_or(&0) >= &batch_size {
                    let mut temp_map = BTreeMap::new();
                    temp_map.insert(file_name.clone(), std::mem::replace(buffer, Vec::with_capacity(batch_size * 4)));
                    tx.send(temp_map).unwrap();
                    buffer_sizes.insert(file_name.clone(), 0);
                }
            }
        }
    }
    
    // 发送剩余的批处理数据
    for (file_name, buffer) in file_buffers.drain() {
        if !buffer.is_empty() {
            let mut temp_map = BTreeMap::new();
            temp_map.insert(file_name, buffer);
            let _ = tx.send(temp_map);
        }
    }

    // Send termination signal to async IO
    let _ = tx.send(BTreeMap::new());
}

// Aync IO线程
async fn process_io(rx: std::sync::mpsc::Receiver<BTreeMap<String, Vec<String>>>) {
    // 文件缓存，避免重复打开/关闭文件
    let mut file_handles: HashMap<String, tokio::fs::File> = HashMap::with_capacity(64);
    
    // 使用缓冲池减少内存分配
    let mut buffer_pool: Vec<Vec<u8>> = Vec::with_capacity(8);
    for _ in 0..8 {
        buffer_pool.push(Vec::with_capacity(131_072)); // 预分配128KB的缓冲区
    }
    
    loop {
        match rx.recv() {
            Ok(reads_batch) => {
                // If we receive an empty BTreeMap, treat it as a termination signal
                if reads_batch.is_empty() {
                    break;
                }
                for (fq_name, fq_data) in reads_batch.iter() {
                    // 获取或创建文件句柄
                    let file = if let Some(f) = file_handles.get_mut(fq_name) {
                        f
                    } else {
                        let new_file = match tokio::fs::OpenOptions::new()
                            .create(true)
                            .append(true)
                            .open(fq_name)
                            .await {
                                Ok(f) => f,
                                Err(e) => {
                                    error!("Failed to open file {}: {}", fq_name, e);
                                    continue;
                                }
                            };
                        file_handles.insert(fq_name.clone(), new_file);
                        file_handles.get_mut(fq_name).unwrap()
                    };
                    
                    // 从池中获取缓冲区
                    let mut buffer = if let Some(buf) = buffer_pool.pop() {
                        buf
                    } else {
                        Vec::with_capacity(131_072)
                    };
                    buffer.clear();
                    
                    // 将数据写入缓冲区
                    for line in fq_data {
                        buffer.extend_from_slice(line.as_bytes());
                        buffer.push(b'\n');
                    }
                    
                    // 写入文件并处理错误
                    if let Err(e) = file.write_all(&buffer).await {
                        error!("Failed to write to file {}: {}", fq_name, e);
                    }
                    
                    // 归还缓冲区到池中
                    buffer_pool.push(buffer);
                }
            }
            Err(_) => break, // 接收通道已关闭
        }
    }
    
    // 关闭所有文件句柄
    for (name, mut file) in file_handles.drain() {
        if let Err(e) = file.flush().await {
            error!("Error flushing file {}: {}", name, e);
        }
    }
}

fn proceed_artidx_map(five_art_fa: &str, three_art_fa: &str,five_idx_fa: &str, three_idx_fa: &str, seed_len: usize) -> 
                                        Result<(usize, BTreeMap<Vec<u8>, String>, BTreeMap<Vec<u8>, String>, 
                                        BTreeMap<Vec<u8>, String>, BTreeMap<Vec<u8>, String>, 
                                        BTreeMap<String, Vec<u8>>, BTreeMap<String, Vec<u8>>,
                                        BTreeMap<String, Vec<(String,u8)>>,  BTreeMap<String, Vec<(String,u8)>>), std::io::Error>{
    debug!("run proceed_artidx_map");
    let start = Instant::now();
    let five_artificial_map = readfa_to_btreemap(&five_art_fa);
    let three_artificial_map = readfa_to_btreemap(&three_art_fa);
    let five_artificial_id_map = readfa_to_id_btreemap(&five_art_fa);
    let three_artificial_id_map = readfa_to_id_btreemap(&three_art_fa);
    let five_index_map = readfa_to_btreemap(&five_idx_fa);
    let three_index_map = readfa_to_btreemap(&three_idx_fa);



    let (five_art_length, five_artificial_idx_map) = fastmap_five_art_idx(&five_artificial_map, &five_index_map);

    let five_art_seed_map = fastmap_art_seed(&five_artificial_idx_map,seed_len);
    let three_art_seed_map = fastmap_art_seed(&three_artificial_map,seed_len);

    info!("valid five_artificial_idx_map number: {}",five_artificial_idx_map.len());
    info!("valid five_index_map number: {}",five_artificial_idx_map.len());
    info!("valid three_artificial_map number: {}",three_artificial_map.len());
    info!("valid three_index_map number: {}",five_artificial_idx_map.len());
    let duration = start.elapsed();
    info!("Time elapsed in proceed_artidx_map() is: {:?}", duration);
    Ok((five_art_length, five_artificial_idx_map, five_index_map, three_artificial_map, three_index_map, 
        five_artificial_id_map, three_artificial_id_map, five_art_seed_map, three_art_seed_map))
}

fn pretrain(pretrain_reads_number: usize, five_art_length: usize, five_artificial_idx_map: BTreeMap<Vec<u8>, String>,
    five_index_map: BTreeMap<Vec<u8>, String>, 
    three_artificial_map: BTreeMap<Vec<u8>, String>, three_index_map: BTreeMap<Vec<u8>, String>,
    five_artificial_id_map: BTreeMap<String, Vec<u8>>, three_artificial_id_map: BTreeMap<String, Vec<u8>>,
    five_art_seed_map: BTreeMap<String, Vec<(String,u8)>>, three_art_seed_map: BTreeMap<String, Vec<(String,u8)>>,
    idx_loc: u8, pe1_file: &str, pe2_file: &str,
    seed_length: usize, error_tolerance: u8,
    match_score: i8, mismatch_score: i8,
    gap_open_score: i8, gap_extend_score: i8) -> Result<BTreeMap<Vec<u8>, String>, std::io::Error> {
    
    if pretrain_reads_number==0 {
        return Ok(five_artificial_idx_map)
    }
    debug!("run pretrain");
    let start = Instant::now();

    let mut pe1_reads: Vec<String> = read_fastq(&pe1_file, pretrain_reads_number);
    let (_, five_artificial_idx_map_ext1) = classify_reads(&five_artificial_idx_map, &three_artificial_map,
        &five_index_map, &three_index_map, &five_artificial_id_map, &three_artificial_id_map, 
        &five_art_seed_map, &three_art_seed_map, five_art_length, idx_loc & 0b1 == 0b1, 
        &mut pe1_reads, seed_length, error_tolerance, match_score, mismatch_score, 
        gap_open_score, gap_extend_score, false);
        
    let mut five_artificial_idx_map_ext2 = BTreeMap::new();
    
    // Only process PE2 if it exists
    if !pe2_file.is_empty() && Path::new(pe2_file).exists() {
        let mut pe2_reads: Vec<String> = read_fastq(&pe2_file, pretrain_reads_number);
        let (_, map_ext2) = classify_reads(&five_artificial_idx_map, &three_artificial_map,
            &five_index_map, &three_index_map, &five_artificial_id_map, &three_artificial_id_map, 
            &five_art_seed_map, &three_art_seed_map, five_art_length, idx_loc & 0b10 == 0b10, 
            &mut pe2_reads, seed_length, error_tolerance, match_score, mismatch_score, 
            gap_open_score, gap_extend_score, false);
        five_artificial_idx_map_ext2 = map_ext2;
    }

    // Combine maps
    let mut five_artificial_idx_map_ext: BTreeMap<Vec<u8>, String> = BTreeMap::new();
    for (k, v) in five_artificial_idx_map.iter() {
        five_artificial_idx_map_ext.insert(k.clone(), v.clone());
    }
    for (k, v) in five_artificial_idx_map_ext1.iter() {
        five_artificial_idx_map_ext.insert(k.clone(), v.clone());
    }
    for (k, v) in five_artificial_idx_map_ext2.iter() {
        five_artificial_idx_map_ext.insert(k.clone(), v.clone());
    }

    info!("pretrain art-idx-map number: {}", five_artificial_idx_map_ext.len());

    let duration = start.elapsed();
    info!("Time elapsed in pretrain() is: {:?}", duration);
    Ok(five_artificial_idx_map_ext)
} 
        

fn seg_pe(five_art_length: usize, five_artificial_idx_map: BTreeMap<Vec<u8>, String>,five_index_map: BTreeMap<Vec<u8>, String>, 
    three_artificial_map: BTreeMap<Vec<u8>, String>, three_index_map: BTreeMap<Vec<u8>, String>,
    five_artificial_id_map: BTreeMap<String, Vec<u8>>, three_artificial_id_map: BTreeMap<String, Vec<u8>>,
    five_art_seed_map: BTreeMap<String, Vec<(String,u8)>>, three_art_seed_map: BTreeMap<String, Vec<(String,u8)>>,
    idx_loc: u8, pe1_file: &str, pe2_file: &str,
    seed_length: usize, merge_pe_bool: bool, error_tolerance: u8,
    match_score: i8, mismatch_score: i8,
    gap_open_score: i8, gap_extend_score: i8, num_threads: usize, batch: u32,
    trim_name: bool, outdir: &str) -> std::io::Result<()> {
    debug!("run seg_pe");
    let start = Instant::now();
    let mut thread_error: Option<std::io::Error> = None;

    let out_path = Path::new(outdir);
    if !out_path.exists() {
        match create_dir(&out_path) {
            Ok(_) => info!("Directory created {}", out_path.display()),
            Err(e) => error!("Error creating directory: {}", e),
        }
    }

    // 设定每个单元包含的reads数量
    let reads_per_unit = batch as usize;

    // 检查PE2文件是否存在
    let pe2_exists = !pe2_file.is_empty() && Path::new(pe2_file).exists();
    
    // 打开FASTQ文件
    let file1 = File::open(pe1_file)?;
    
    // 获取file1 prefix_name
    let prefix1 = outdir.to_string() + "/" + &get_fastq_prefix(pe1_file);
    
    // 仅在PE2存在时获取file2和prefix2
    let file2_option = if pe2_exists {
        Some(File::open(pe2_file)?)
    } else {
        None
    };
    
    let prefix2 = if pe2_exists {
        outdir.to_string() + "/" + &get_fastq_prefix(pe2_file)
    } else {
        String::new()
    };

    // Merge PE标记 - 仅当PE2存在时才开启
    let actual_merge_pe = merge_pe_bool && pe2_exists;

    // 设置worker线程的channels
    let mut channels = Vec::new();
    let (io_tx, io_rx) = std::sync::mpsc::channel();

    // 共享数据结构
    let five_artificial_idx_map = Arc::new(five_artificial_idx_map);
    let three_artificial_map = Arc::new(three_artificial_map);
    let five_index_map = Arc::new(five_index_map);
    let three_index_map = Arc::new(three_index_map);
    let five_artificial_id_map = Arc::new(five_artificial_id_map);
    let three_artificial_id_map = Arc::new(three_artificial_id_map);
    let five_art_seed_map = Arc::new(five_art_seed_map);
    let three_art_seed_map = Arc::new(three_art_seed_map);

    // 创建工作线程 - 根据是否为PE模式创建不同类型的线程
    let mut handles = Vec::new();
    for _ in 0..num_threads {
        let five_artificial_idx_map = Arc::clone(&five_artificial_idx_map);
        let three_artificial_map = Arc::clone(&three_artificial_map);
        let five_index_map = Arc::clone(&five_index_map);
        let three_index_map = Arc::clone(&three_index_map);
        let five_artificial_id_map = Arc::clone(&five_artificial_id_map);
        let three_artificial_id_map = Arc::clone(&three_artificial_id_map);
        let five_art_seed_map = Arc::clone(&five_art_seed_map);
        let three_art_seed_map = Arc::clone(&three_art_seed_map);
        let (tx, rx) = std::sync::mpsc::channel::<(Vec<String>, Vec<String>)>();
        let rx = Arc::new(Mutex::new(rx));
        channels.push((tx, Arc::clone(&rx)));

        let prefix1_clone = prefix1.clone();
        let prefix2_clone = prefix2.clone();
        let io_tx_clone = Arc::new(Mutex::new(io_tx.clone()));
        
        // 根据模式创建不同类型的处理线程
        if pe2_exists {
            // 双端模式处理线程
            let handle = thread::spawn(move || {
                for pair in rx.lock().unwrap().iter() {
                    if pair.0.is_empty() && pair.1.is_empty() {
                        break; // 收到终止信号，退出线程
                    }
                    process_reads_pe(&five_artificial_idx_map, &three_artificial_map, &five_index_map, &three_index_map,  
                        &five_artificial_id_map, &three_artificial_id_map,
                        &five_art_seed_map, &three_art_seed_map,
                        five_art_length, idx_loc,
                        &prefix1_clone, &prefix2_clone,
                        pair, seed_length, actual_merge_pe, 
                        error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score,
                        io_tx_clone.lock().unwrap().clone(), trim_name);
                }
            });
            handles.push(handle);
        } else {
            // 单端模式处理线程
            let handle = thread::spawn(move || {
                for reads1 in rx.lock().unwrap().iter() {
                    // 在单端模式下，我们接收的是(reads1, empty_vec)对，只取第一个元素
                    process_reads_se(&five_artificial_idx_map, &three_artificial_map, &five_index_map, &three_index_map,  
                        &five_artificial_id_map, &three_artificial_id_map,
                        &five_art_seed_map, &three_art_seed_map,
                        five_art_length, idx_loc,
                        &prefix1_clone, reads1.0, seed_length,
                        error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score,
                        io_tx_clone.lock().unwrap().clone(), trim_name);
                }
            });
            handles.push(handle);
        }
    }

    // Create cache channel thread
    let (tx, rx) = std::sync::mpsc::channel();
    let cache_handle = thread::spawn(move || {
        cache_channel(io_rx, tx, (batch/10) as usize);
    });

    // Create async IO thread
    let rt = tokio::runtime::Runtime::new().unwrap();
    let io_handle = rt.spawn(async move {
        process_io(rx).await;
    });

    // Memory-map the PE1 file
    let mmap1 = unsafe { MmapOptions::new().map(&file1)? };
    
    // 为不同的输入格式优化读取策略
    let reader1: Box<dyn std::io::BufRead> = if pe1_file.ends_with(".gz") || pe1_file.ends_with(".gzip") {
        // 为gzip文件使用大缓冲区以提高解压效率
        Box::new(BufReader::with_capacity(131_072, GzDecoder::new(&mmap1[..])))
    } else {
        // 对于未压缩文件，尝试直接处理内存映射
        Box::new(BufReader::with_capacity(262_144, &mmap1[..]))
    };

    let mut t: usize = 0;
    let mut m: usize = 0;

    if pe2_exists {
        // 双端模式优化
        let file2 = file2_option.unwrap();
        let mmap2 = unsafe { MmapOptions::new().map(&file2)? };
        let reader2: Box<dyn std::io::BufRead> = if pe2_file.ends_with(".gz") || pe2_file.ends_with(".gzip") {
            Box::new(BufReader::with_capacity(131_072, GzDecoder::new(&mmap2[..])))
        } else {
            Box::new(BufReader::with_capacity(262_144, &mmap2[..]))
        };

        // 为PE读取预分配缓冲区
        let mut pe1_buffer: Vec<String> = Vec::with_capacity(reads_per_unit * 4);
        let mut pe2_buffer: Vec<String> = Vec::with_capacity(reads_per_unit * 4);
        
        // 追踪已分配的线程工作
        let mut thread_load = vec![0; num_threads];

        for (line1, line2) in reader1.lines().zip(reader2.lines()) {
            let line1 = line1.unwrap();
            let line2 = line2.unwrap();
            pe1_buffer.push(line1);
            pe2_buffer.push(line2);
            m += 1;
            
            if m % (reads_per_unit * 4) == 0 {
                // 选择负载最低的线程
                t = thread_load.iter().enumerate()
                     .min_by_key(|&(_, load)| load)
                     .map(|(idx, _)| idx)
                     .unwrap_or(0);
                
                channels[t].0.send((pe1_buffer, pe2_buffer)).unwrap();
                thread_load[t] += 1;
                
                // 重新分配新的缓冲区
                pe1_buffer = Vec::with_capacity(reads_per_unit * 4);
                pe2_buffer = Vec::with_capacity(reads_per_unit * 4);
            }
        }
        
        // 处理剩余数据
        if !pe1_buffer.is_empty() {
            channels[t].0.send((pe1_buffer, pe2_buffer)).unwrap();
        }
        
        info!("Total handle paired reads: {} x 2", m/4);
    } else {
        // 单端模式优化 - 类似双端模式的优化
        // 预分配单端模式的缓冲区
        let mut pe1_buffer: Vec<String> = Vec::with_capacity(reads_per_unit * 4);
        let empty_buffer: Vec<String> = Vec::new();
        
        // 追踪已分配的线程工作以实现负载均衡
        let mut thread_load = vec![0; num_threads];

        // 使用优化的缓冲区读取方式
        let mut lines = reader1.lines();
        let mut batch_size = 0;
        
        while let Some(Ok(line)) = lines.next() {
            pe1_buffer.push(line);
            m += 1;
            batch_size += 1;
            
            // 每读取完一个完整的reads批次，就分配给负载最低的线程
            if batch_size == reads_per_unit * 4 {
                // 选择负载最低的线程处理当前批次
                t = thread_load.iter()
                    .enumerate()
                    .min_by_key(|&(_, load)| load)
                    .map(|(idx, _)| idx)
                    .unwrap_or(0);
                
                // 直接移动所有权，不再复制数据
                let buffer_to_send = std::mem::replace(&mut pe1_buffer, Vec::with_capacity(reads_per_unit * 4));
                
                // 发送数据到处理线程
                channels[t].0.send((buffer_to_send, empty_buffer.clone())).unwrap();
                thread_load[t] += 1;
                
                // 重置批次计数
                batch_size = 0;
            }
        }
        
        // 处理剩余的不完整批次
        if !pe1_buffer.is_empty() {
            // 优先发送给负载最低的线程
            t = thread_load.iter()
                .enumerate()
                .min_by_key(|&(_, load)| load)
                .map(|(idx, _)| idx)
                .unwrap_or(0);
                
            channels[t].0.send((pe1_buffer, empty_buffer)).unwrap();
        }
        
        info!("Total handled single-end reads: {}", m/4);
    }

    // Close channels
    for channel in &channels {
        let _ = channel.0.send((Vec::new(), Vec::new())); 
    }

    // 关闭所有发送端
    for channel in channels {
        drop(channel); // 关闭发送端
    }

    // 等待所有线程完成
    for handle in handles {
        let _ = handle.join();
    }

    drop(io_tx);

    cache_handle.join().unwrap();
    match rt.block_on(io_handle) {
        Ok(_) => {},
        Err(e) => {
            error!("Error in async IO: {:?}", e);
            thread_error = Some(std::io::Error::new(std::io::ErrorKind::Other, "Error in async IO"));
        }
    }

    // 检查线程是否出现错误
    if let Some(error) = thread_error {
        return Err(error);
    }

    let duration = start.elapsed();
    info!("Time elapsed in seg_pe() is: {:?}", duration);
    Ok(())
}

// 添加SIMD优化的序列处理函数

#[allow(unused_imports)]
#[cfg(target_arch = "x86_64")]
mod simd_ops {
    // No need to import all, just what we need
    use std::arch::x86_64::{__m256i, _mm256_cmpeq_epi8}; // Import only what's needed
    
    // SIMD optimized reverse complement
    pub unsafe fn reverse_complement_avx2(input: &str) -> String {
        let len = input.len();
        // Create result vector but only use it when needed
        let _result = vec![0u8; len]; // Use _ to indicate intentionally unused for now
        
        // Simplified implementation for now
        // This is a placeholder - real implementation would use AVX2
        super::reverse_complement(input)
    }
    
    // SIMD optimized sequence search
    pub unsafe fn sequence_search_avx2(haystack: &str, needle: &str) -> Vec<i16> {
        let _needle_len = needle.len();
        let _haystack_len = haystack.len();
        let _positions: Vec<i16> = Vec::new();
        
        // Simple implementation for now
        super::regrex_sarch(haystack, needle)
    }
}

// 优化的反向互补函数
pub fn fast_reverse_complement(input: &str) -> String {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            return unsafe { simd_ops::reverse_complement_avx2(input) };
        }
    }
    
    // 回退到标准实现
    reverse_complement(input)
}

// 优化的序列搜索函数
pub fn fast_sequence_search(haystack: &str, needle: &str) -> Vec<i16> {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            return unsafe { simd_ops::sequence_search_avx2(haystack, needle) };
        }
    }
    
    // 回退到标准实现
    regrex_sarch(haystack, needle)
}

fn main() {
    // Initialize logging
    env_logger::init();

    let mut args = Args::parse();  // Add 'mut' keyword here
    
    // Validate that PE1 exists
    if !Path::new(&args.pe1_fastq).exists() {
        error!("PE1 FASTQ file not found: {}", args.pe1_fastq);
        std::process::exit(1);
    }
    
    // Log whether we're in single or paired-end mode
    if args.pe2_fastq.is_empty() {
        info!("Running in single-end mode with PE1 file only");
    } else if !Path::new(&args.pe2_fastq).exists() {
        warn!("PE2 FASTQ file not found: {}. Running in single-end mode.", args.pe2_fastq);
        // Simply modify the pe2_fastq field directly instead of creating a new struct
        args.pe2_fastq = String::new();
    }

    let proceed_result = proceed_artidx_map(&args.five_art_fa, &args.three_art_fa, &args.five_idx_fa, &args.three_idx_fa, args.seed_len);
    if let Ok((five_art_length, five_artificial_idx_map, five_index_map, three_artificial_map, three_index_map, 
        five_artificial_id_map, three_artificial_id_map, five_art_seed_map, three_art_seed_map)) = proceed_result {
        
        let pretrain_result = pretrain(args.train, five_art_length, five_artificial_idx_map, five_index_map.clone(), 
            three_artificial_map.clone(), three_index_map.clone(), five_artificial_id_map.clone(), three_artificial_id_map.clone(), 
            five_art_seed_map.clone(), three_art_seed_map.clone(), 
            args.idx_loc, &args.pe1_fastq, &args.pe2_fastq, args.seed_len,
            args.error_tolerance, args.match_score, args.error_score, args.gap_open_score, args.gap_extend_score);
        
        match pretrain_result {
            Ok(five_artificial_idx_map_ext) => {
                match seg_pe(five_art_length, five_artificial_idx_map_ext, five_index_map, three_artificial_map, three_index_map,  
                    five_artificial_id_map, three_artificial_id_map, five_art_seed_map, three_art_seed_map,
                    args.idx_loc, &args.pe1_fastq, &args.pe2_fastq, args.seed_len,
                    args.merge_pe, args.error_tolerance, args.match_score, args.error_score,
                    args.gap_open_score, args.gap_extend_score, args.num_threads, args.batch, args.trim_name, &args.outdir) {
                    
                    Ok(_) => info!("Processing completed successfully"),
                    Err(err) => {
                        error!("Error during processing: {}", err);
                        std::process::exit(1);
                    }
                }
            },
            Err(err) => {
                error!("Error during pretrain: {}", err);
                std::process::exit(1);
            }
        }
    } else {
        error!("Error during map initialization: {:?}", proceed_result.err());
        std::process::exit(1);
    }
}
