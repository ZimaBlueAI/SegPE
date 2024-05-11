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
use log::{debug,info,error};
use env_logger;
use std::time::Instant;

mod utils;
use utils::*;

#[derive(Parser, Debug)]
#[command(
    author = "BENM",
    version = "version 0.1.14",
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
    #[arg(long)]
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
    gap_open_score: i8, gap_extend_score: i8)  -> (Vec<String>, Vec<String>) {
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
        // 将reads2反向互补
        let read2_rc = reverse_complement(&reads2[i*4+1]);
        let qual2_r: String = reverse(&reads2[i*4+3]);
        // 从read2_rc 5'端第一个碱基开始到第二个碱基（只允许一个错配）取seed_length长度的序列在read1中查找是否有匹配的位置
        let mut score = 0;
        for j in 0..2 {
            let seed = &read2_rc[j..j+seed_length];
            let match_pos = regrex_sarch(&read1, seed);
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
        // let mut cls5 = String::new(); // Change the type of cls5 to String
        if read_len > five_art_length {
            let trimed_seq_bytes = (&reads[i*4+1][0..five_art_length]).as_bytes();
            // if five_artificial_idx_map.contains_key(trimed_seq_bytes) {
            //     left_trim_pos = five_art_length;
            //     cls5 = five_artificial_idx_map[trimed_seq_bytes].to_owned(); // Use to_owned() to convert the value to String
            //     cls5_bool = true;
            // } else 
            if five_artificial_idx_map_ext.contains_key(trimed_seq_bytes) {
                left_trim_pos = five_art_length;
                cls5 = five_artificial_idx_map_ext[trimed_seq_bytes].to_owned(); // Use to_owned() to convert the value to String
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
        if three_artificial_map.len() > 0 {
            for j in 0..1+error_tolerance as usize {
                let read_len = read_len;
                let seed = &reads[i*4+1][read_len-seed_length-j..read_len-j];
                if three_art_seed_map.contains_key(seed) {
                    let vec_map = three_art_seed_map.get(seed); 
                    for map in vec_map.unwrap() {
                        let map_id3 = map.0.clone(); 
                        if three_art_seed_map.contains_key(&map_id3) {
                            let art_seq = String::from_utf8(three_artificial_id_map.get(&map_id3).unwrap().to_owned()).unwrap();
                            let pos = map.1 as usize;
                            let overlapping_len = pos+seed_length+j;                                                              
                            if overlapping_len>=seed_length && pos+overlapping_len < art_seq.len() && read_len-seed_length-j-pos>0 && read_len-seed_length-j-pos<reads[i*4+1].len(){
                                let res = block_align(&art_seq[0..pos+overlapping_len], &reads[i*4+1][read_len-seed_length-j-pos..], match_score, mismatch_score, gap_open_score, gap_extend_score);
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
    if merge_pe_bool {
        let new_pair = merge_pe(pair, seed_length, error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score);
        reads1 = new_pair.0.clone();
        reads2 = new_pair.1.clone();
    } 
    // let mut reads1 = pair.0.clone();
    // let mut reads2 = pair.1.clone();
    let (cls1, _) = classify_reads(five_artificial_idx_map, three_artificial_map,
        five_index_map, three_idx_map, five_artificial_id_map, three_artificial_id_map, five_art_seed_map, three_art_seed_map,
        five_art_length, idx_loc & 0b1 == 0b1, &mut reads1, seed_length,
        error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score, trim_name);
    let (cls2, _) = classify_reads(five_artificial_idx_map, three_artificial_map,
        five_index_map, three_idx_map, five_artificial_id_map, three_artificial_id_map, five_art_seed_map, three_art_seed_map, 
        five_art_length, idx_loc & 0b10 == 0b10, &mut reads2, seed_length,
        error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score, trim_name);

    for i in 0..cls1.len() {
        let mut cls = "";

        if reads2[i*4+1] == "" {
            if idx_loc & 0b1 == 0b1 && cls1[i].0 != "" {
                cls = &cls1[i].0;
            }else if idx_loc & 0b10 == 0b10 && cls1[i].1 != "" {
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
            if reads2[i*4+1]!="" {
                let fq1_name = format!("{}_{}.fastq", prefix1.to_string(), cls);
                let fq2_name  = format!("{}_{}.fastq", prefix2.to_string(), cls);
                let fq1_data = reads1[i*4+0..i*4+4].to_vec();
                let fq2_data = reads2[i*4+0..i*4+4].to_vec();
                // let cls_string = cls.to_string(); // Convert cls to String
                io_tx.send(vec![PEreads{read1_name: fq1_name, read1_vec: fq1_data, read2_name: fq2_name, read2_vec: fq2_data}]).expect("Failed to send data through channel");
            } else {
                let suffix = cls.to_owned() + "_single";
                let fq1_name = format!("{}_{}.fastq", prefix1.to_string(), suffix);
                let fq1_data = reads1[i*4+0..i*4+4].to_vec();
                io_tx.send(vec![PEreads{read1_name: fq1_name, read1_vec: fq1_data, read2_name: "".to_string(), read2_vec: Vec::new()}]).expect("Failed to send data through channel");
            }
        } else {
            if reads2[i*4+1]!="" {
                let suffix =  "unclassified";
                let fq1_name = format!("{}_{}.fastq", prefix1.to_string(), suffix);
                let fq2_name  = format!("{}_{}.fastq", prefix2.to_string(), suffix);
                let fq1_data = reads1[i*4+0..i*4+4].to_vec();
                let fq2_data = reads2[i*4+0..i*4+4].to_vec();
                io_tx.send(vec![PEreads{read1_name: fq1_name, read1_vec: fq1_data, read2_name: fq2_name, read2_vec: fq2_data}]).expect("Failed to send data through channel");
            } else{
                let suffix =  "unclassified";
                let fq1_name = format!("{}_{}.fastq", prefix1.to_string(), suffix);
                let fq1_data = reads1[i*4+0..i*4+4].to_vec();
                io_tx.send(vec![PEreads{read1_name: fq1_name, read1_vec: fq1_data, read2_name: "".to_string(), read2_vec: Vec::new()}]).expect("Failed to send data through channel");
            }
        }
    }

    let duration = start.elapsed();
    info!("Time elapsed in process_reads() is: {:?}", duration);
}


// 创建一个快速读取io_rx，当pe1_fastq/pe2_fastq/single_fastq数据超过batch_size行时，将数据发出给Asyc IO线程，以防止io_rx阻塞
fn cache_channel(io_rx: std::sync::mpsc::Receiver<Vec<PEreads>>, tx: std::sync::mpsc::Sender<BTreeMap<String, Vec<String>>>, batch_size: usize) {
    let mut count1 = 0;
    let mut count2 = 0;

    let mut pe1_batch: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut pe2_batch: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut single_batch: BTreeMap<String, Vec<String>> = BTreeMap::new();

    for pereads in io_rx.iter() {
        for peread in pereads {
            if peread.read2_name != "" {
                if pe1_batch.get(&peread.read1_name).is_some() {
                    pe1_batch
                        .get_mut(&peread.read1_name)
                        .unwrap()
                        .extend(peread.read1_vec.clone());
                    pe2_batch
                        .get_mut(&peread.read2_name)
                        .unwrap()
                        .extend(peread.read2_vec.clone());
                } else {
                    pe1_batch.insert(peread.read1_name.clone(), peread.read1_vec.clone());
                    pe2_batch.insert(peread.read2_name.clone(), peread.read2_vec.clone());
                }
                count2 += 1;
            } else {
                if single_batch.get(&peread.read1_name).is_some() {
                    single_batch
                        .get_mut(&peread.read1_name)
                        .unwrap()
                        .append(&mut peread.read1_vec.clone());
                } else {
                    single_batch.insert(peread.read1_name.clone(), peread.read1_vec.clone());
                }
                count1 += 1;
            }

            if count2 >= batch_size || count2 >= 100 {
                let _ = tx.send(pe1_batch);
                let _ = tx.send(pe2_batch);
                pe1_batch = BTreeMap::new();
                pe2_batch = BTreeMap::new();
                count2 = 0;
            }
            if count1 >= batch_size || count1 >= 100 {
                let _ = tx.send(single_batch);
                count1 = 0;
                single_batch = BTreeMap::new();
            }
        }
    }

    // 发送剩余的批处理数据给Async IO线程
    if pe1_batch.len() > 0 {
        let _ = tx.send(pe1_batch);
    }
    if pe2_batch.len() > 0 {
        let _ = tx.send(pe2_batch);
    }
    if single_batch.len() > 0 {
        let _ = tx.send(single_batch);
    }
}

// Aync IO线程
async fn process_io(rx: std::sync::mpsc::Receiver<BTreeMap<String, Vec<String>>>) {
    loop {
        match rx.recv() {
            Ok(reads_batch) => {
                let fastqwriter = FastqWriter {};
                for (fq_name, fq_data) in reads_batch.iter() {
                    let fq_name = fq_name.to_string(); // Convert fq_name to String
                    let _ = fastqwriter.write_to_fastq(fq_data.clone(), &fq_name).await;
                }
            }
            Err(_) => break, //receiving on a closed channel
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
    let mut pe2_reads: Vec<String>  = read_fastq(&pe2_file, pretrain_reads_number);
    let (_, five_artificial_idx_map_ext1) = classify_reads(&five_artificial_idx_map, &three_artificial_map,
        &five_index_map, &three_index_map,  &five_artificial_id_map, &three_artificial_id_map, &five_art_seed_map, &three_art_seed_map, five_art_length, idx_loc & 0b1 == 0b1, &mut pe1_reads, seed_length,
        error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score, false);    
    let (_,  five_artificial_idx_map_ext2) = classify_reads(&five_artificial_idx_map, &three_artificial_map,
        &five_index_map, &three_index_map, &five_artificial_id_map, &three_artificial_id_map, &five_art_seed_map, &three_art_seed_map, five_art_length, idx_loc & 0b10 == 0b10, &mut pe2_reads, seed_length,
        error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score, false);

    //TODO: combine five_artificial_idx_map_ext1 and five_artificial_idx_map_ext2
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

    info!("pretrain art-idx-map number: {}",five_artificial_idx_map_ext.len());

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

    let out_path = Path::new(outdir);
    if !out_path.exists() {
        match create_dir(&out_path) {
            Ok(_) => info!("Directory created {}", out_path.display()),
            Err(e) => error!("Error creating directory: {}", e),
        }
    }

    // 设定每个单元包含的reads数量
    let reads_per_unit = batch as usize;

    // 打开FASTQ文件
    let file1 = File::open(pe1_file)?;
    let file2 = File::open(pe2_file)?;

    // 获取file1, file2的prefix_names
    let prefix1 = outdir.to_string() + "/" + &get_fastq_prefix(pe1_file);
    let prefix2 = outdir.to_string() + "/" + &get_fastq_prefix(pe2_file);

    let mut channels = Vec::new();
    let (io_tx, io_rx) = std::sync::mpsc::channel();

    let five_artificial_idx_map = Arc::new(five_artificial_idx_map);
    let three_artificial_map = Arc::new(three_artificial_map);
    let five_index_map = Arc::new(five_index_map);
    let three_index_map = Arc::new(three_index_map);
    let five_artificial_id_map = Arc::new(five_artificial_id_map);
    let three_artificial_id_map = Arc::new(three_artificial_id_map);
    let five_art_seed_map = Arc::new(five_art_seed_map);
    let three_art_seed_map = Arc::new(three_art_seed_map);

    let mut handles = Vec::new();
    for _ in 0..num_threads {
        let five_artificial_idx_map: Arc<BTreeMap<Vec<u8>, String>> = Arc::clone(&five_artificial_idx_map);
        let three_artificial_map = Arc::clone(&three_artificial_map);
        let five_index_map = Arc::clone(&five_index_map);
        let three_index_map = Arc::clone(&three_index_map);
        let five_artificial_id_map: Arc<BTreeMap<String, Vec<u8>>> = Arc::clone(&five_artificial_id_map);
        let three_artificial_id_map: Arc<BTreeMap<String, Vec<u8>>>  = Arc::clone(&three_artificial_id_map);
        let five_art_seed_map = Arc::clone(&five_art_seed_map);
        let three_art_seed_map = Arc::clone(&three_art_seed_map);
        let (tx, rx) = std::sync::mpsc::channel();
        let rx = Arc::new(Mutex::new(rx)); // Wrap the receiver in a Mutex
        channels.push((tx, Arc::clone(&rx))); // Clone the receiver instead of moving it

        let prefix1_clone = prefix1.clone(); // Clone prefix1
        let prefix2_clone = prefix2.clone(); // Change the type of prefix2_clone to Arc<String>
        let io_tx_clone = Arc::new(Mutex::new(io_tx.clone())); // Change the type of io_tx_clone to Arc<Mutex<Sender<(String, Vec<String>, String, Vec<String>)>>>
        
        let handle = thread::spawn(move || {
            for pair in rx.lock().unwrap().iter() { // Lock the mutex before accessing the receiver
                process_reads(&five_artificial_idx_map, &three_artificial_map, &five_index_map, &three_index_map,  
                    &five_artificial_id_map, &three_artificial_id_map,
                    &five_art_seed_map, &three_art_seed_map,
                    five_art_length, idx_loc,
                    &prefix1_clone, &prefix2_clone,
                    pair, seed_length, merge_pe_bool, 
                    error_tolerance, match_score, mismatch_score, gap_open_score, gap_extend_score,
                    io_tx_clone.lock().unwrap().clone(), trim_name);
            }
        });
        handles.push(handle);
    }

    // 创建单线程cache channel
    let (tx, rx) = std::sync::mpsc::channel();
    let cache_handle = thread::spawn(move || {
        cache_channel(io_rx, tx, (batch/10) as usize);
    });

    // 创建Asnyc IO线程
    let rt = tokio::runtime::Runtime::new().unwrap();
    let io_handle = rt.spawn(async move {
        process_io(rx).await;
    });

    // 按内在空间容量适当范围读mmap1, mmap2，发送给处理线程
    // Memory-map the files
    let mmap1 = unsafe { MmapOptions::new().map(&file1)? };
    let mmap2 = unsafe { MmapOptions::new().map(&file2)? };
    let reader1: Box<dyn std::io::BufRead> = if pe1_file.ends_with(".gz") || pe1_file.ends_with(".gzip") {
        Box::new(BufReader::with_capacity(8192, GzDecoder::new(&mmap1[..])))
    } else {
        Box::new(BufReader::with_capacity(8192, &mmap1[..]))
    };
    let reader2: Box<dyn std::io::BufRead> = if pe2_file.ends_with(".gz") || pe2_file.ends_with(".gzip") {
        Box::new(BufReader::with_capacity(8192, GzDecoder::new(&mmap2[..])))
    } else {
        Box::new(BufReader::with_capacity(8192, &mmap2[..]))
    };

    let mut t: usize = 0;
    let mut m: usize = 0;
    let mut pe1_buffer: Vec<String> = Vec::new();
    let mut pe2_buffer: Vec<String>  = Vec::new();

    for (line1, line2) in reader1.lines().zip(reader2.lines()) {
        let line1 = line1.unwrap();
        let line2 = line2.unwrap();
        pe1_buffer.push(line1);
        pe2_buffer.push(line2);
        m += 1;
        if m % (reads_per_unit * 4) == 0 {
            channels[t].0.send((pe1_buffer.clone(), pe2_buffer.clone())).unwrap();
            pe1_buffer.clear();
            pe2_buffer.clear();
            t = (t + 1) % num_threads;
        }
    }
    info!("Total handle paired reads: {} x 2", m);

    // 发送最后一个块（如果有的话）
    if !pe1_buffer.is_empty() {
        channels[t].0.send((pe1_buffer.clone(), pe2_buffer.clone())).unwrap();
    }

    // 关闭所有发送端
    for chnnel in channels {
        drop(chnnel); // 关闭发送端
    }

    // 等待所有处理线程完成
    for handle in handles {
        handle.join().unwrap();
    }

    drop(io_tx);

    let _ = tokio::runtime::Runtime::new().unwrap().block_on(io_handle);

    cache_handle.join().unwrap();

    let duration = start.elapsed();
    info!("Time elapsed in seg_pe() is: {:?}", duration);
    Ok(())
}

fn main() {
    // 读参
    env_logger::init();

    let args = Args::parse();
    let proceed_result = proceed_artidx_map(&args.five_art_fa, &args.three_art_fa, &args.five_idx_fa, &args.three_idx_fa, args.seed_len);
    if let Ok((five_art_length, five_artificial_idx_map, five_index_map, three_artificial_map, three_index_map, 
        five_artificial_id_map, three_artificial_id_map, five_art_seed_map, three_art_seed_map)) = proceed_result {
        let pretrain_result = pretrain(args.train, five_art_length, five_artificial_idx_map, five_index_map.clone(), 
        three_artificial_map.clone(), three_index_map.clone(), five_artificial_id_map.clone(), three_artificial_id_map.clone(), 
        five_art_seed_map.clone(), three_art_seed_map.clone(), 
        args.idx_loc, &args.pe1_fastq, &args.pe2_fastq, args.seed_len,
        args.error_tolerance, args.match_score, args.error_score, args.gap_open_score, args.gap_extend_score);
        if let Ok(five_artificial_idx_map_ext) = pretrain_result {
            if let Err(err) = seg_pe(five_art_length, five_artificial_idx_map_ext, five_index_map, three_artificial_map, three_index_map,  
                five_artificial_id_map, three_artificial_id_map, five_art_seed_map, three_art_seed_map,
                args.idx_loc, &args.pe1_fastq, &args.pe2_fastq, args.seed_len,
                args.merge_pe, args.error_tolerance, args.match_score, args.error_score,
                args.gap_open_score, args.gap_extend_score, args.num_threads, args.batch, args.trim_name, &args.outdir) {
                eprintln!("Error: {}", err);
                std::process::exit(1);
            }
        } else if let Err(err) = pretrain_result {
            eprintln!("Error: {}", err);
            std::process::exit(1);
        }
    } else if let Err(err) = proceed_result {
        eprintln!("Error: {}", err);
        std::process::exit(1);
    }
    
}
