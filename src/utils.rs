use std::{
    path::Path,
    collections::BTreeMap,
    io::{BufWriter,Write},
    future
};
use regex::Regex;

use bio::alignment::{
    pairwise::*,
    sparse::*
};
use bio::io::fasta;
use fastq::{each_zipped, parse_path};
use needletail::parse_fastx_file;

use block_aligner::{cigar::*, scan_block::*, scores::*};

#[allow(dead_code)]
pub fn read_fastq(fastq_path: &str, cutoff_reads_number: usize) -> Vec<String>{
    let mut reader = parse_fastx_file(&fastq_path).expect("valid path/file");
    let mut reads_buffer = Vec::new();
    let mut m = 0;
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let seq = seqrec.seq().to_owned().to_vec();
        let qual = seqrec.qual().unwrap_or_default().to_vec();
        reads_buffer.push(format!("@{}", String::from_utf8(seqrec.id().to_vec()).unwrap().as_str()));
        reads_buffer.push(String::from_utf8(seq).unwrap());
        reads_buffer.push(String::from("+"));
        reads_buffer.push(String::from_utf8(qual).unwrap());
        m += 1;
        if cutoff_reads_number>0 && m >= cutoff_reads_number * 4 {
            println!("Number of reads {}: {}", fastq_path, m);
            break;
        }
    }
    reads_buffer
}

#[allow(dead_code)]
pub fn read_pair_fastq(pe1_path: &str, pe2_path: &str) -> Result<(Vec<String>, Vec<String>), std::io::Error>{
    let mut pe1 = Vec::new();
    let mut pe2 = Vec::new();
    let mut m = 0;
    let mut n = 0;
    parse_path(Some(&pe1_path), |parser1| {
        parse_path(Some(&pe2_path), |parser2| {
            each_zipped(parser1, parser2, |rec1, rec2| {
                if let Some(rec) = rec1 {
                    let record = rec.to_owned_record();
                    pe1.push(format!("@{}", String::from_utf8(record.head.to_vec()).unwrap().as_str()));
                    pe1.push(String::from_utf8(record.seq.to_vec()).unwrap());
                    pe1.push(String::from("+"));
                    pe1.push(String::from_utf8(record.qual.to_vec()).unwrap());
                    m += 1;
                }
                if let Some(rec) = rec2 {
                    let record = rec.to_owned_record();
                    pe2.push(format!("@{}",String::from_utf8(record.head.to_vec()).unwrap()));
                    pe2.push(String::from_utf8(record.seq.to_vec()).unwrap());
                    pe2.push(String::from("+"));
                    pe2.push(String::from_utf8(record.qual.to_vec()).unwrap());
                    n += 1;
                }
                (true, true)
            })
            .expect("Invalid record.");
        })
        .expect("Unknown format for file 2.");
    })
    .expect("Unknown format for file 1.");
    Ok((pe1, pe2))
}

#[allow(dead_code)]
pub fn get_fastq_prefix(file_path: &str) -> String {
    Path::new(file_path)
        .file_stem()         // 获取文件名（不包括扩展名）
        .and_then(|stem| stem.to_str())  // 将 OsStr 转换为 &str
        .map(|s| s.to_string())          // 将 &str 转换为 String
        .unwrap_or_else(|| "".to_string()) // 如果上述步骤失败，则返回空字符串
}

#[allow(dead_code)]
pub fn reverse(input: &str) -> String {
    input.chars().rev().collect()
}

#[allow(dead_code)]
pub fn reverse_complement(sequence: &str) -> String {
    let reversed_sequence: String = sequence.chars()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
             _ => 'N',
        })
        .collect();

    return reverse(&reversed_sequence);
}

fn write_io(data: Vec<String>,file_name: &str)->future::Ready<std::io::Result<()>> {
    let fastq_name = file_name.to_string();

    let file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&fastq_name)
        .or_else(|_| std::fs::File::create(&fastq_name))
        .expect("Failed to open or create file");

    let mut writer = BufWriter::new(file);

    for line in data {
        writer.write_all(line.as_bytes()).expect("Failed to write to file");
        writer.write_all(b"\n").expect("Failed to write to file");
    }

    writer.flush().expect("Failed to flush file");

    future::ready(Ok(()))
}

pub trait WriteFastq {
    #[allow(dead_code)]
    async fn write_to_fastq(
        &self,
        data: Vec<String>,
        file_name: &str,
    ) -> future::Ready<std::io::Result<()>>;
}

pub struct FastqWriter;

impl WriteFastq for FastqWriter {
    async fn write_to_fastq(
        &self,
        data: Vec<String>,
        file_name: &str,
    ) -> future::Ready<std::io::Result<()>>  {
        match write_io(data,file_name).await {
            Ok(_) => future::ready(Ok(())),
            Err(e) => future::ready(Err(e)),
        }
    }
}


#[allow(dead_code)]
pub fn find_pattern(seq: &str, pattern: &str) -> Vec<usize> {
    seq.match_indices(pattern).map(|(i, _)| i).collect()
}

#[allow(dead_code)]
pub fn hamming_distance(s1: &str, s2: &str) -> usize {
    s1.chars().zip(s2.chars()).filter(|&(a, b)| a != b).count()
}

#[allow(dead_code)]
pub fn regrex_sarch(seq: &str, pattern: &str) -> Vec<i16> {
    let pattern = Regex::new(pattern).unwrap();

    let match_pos: Vec<i16> = pattern
        .find_iter(seq)
        .map(|mat| mat.start() as i16)
        .collect();

    let mut match_pos = match_pos;
    match_pos.reverse();
    match_pos
}

#[allow(dead_code)]
pub fn sw_align(seq1: &str, seq2: &str, match_score: i8, mismatch_score: i8, 
    gap_open_score: i8, gap_extend_score: i8) -> i32 {

    let x = seq1.as_bytes();
    let y = seq2.as_bytes();
    
    let scoring = Scoring {
        gap_open: gap_open_score as i32,
        gap_extend: gap_extend_score as i32,
        match_fn: |a: u8, b: u8| if a == b { 1i32 } else { -3i32 },
        match_scores: Some((match_score as i32, mismatch_score as i32)),
        xclip_prefix: -10,
        xclip_suffix: MIN_SCORE,
        yclip_prefix: 0,
        yclip_suffix: 0,
    };
    let mut aligner = Aligner::with_capacity_and_scoring(x.len(), y.len(), scoring);
    let alignment = aligner.custom(x, y);
    // 通过alignment.pretty获取比对的错误mismatch, indels(gap_open, gap_extend)值
    // return (alignment.score, alignment.cigar(false));
    return alignment.score
}

#[allow(dead_code)]
pub fn sparse_align(seq1: &str, seq2: &str, kmer: usize) -> i32 {
    let matches = find_kmer_matches(seq1.as_bytes(), seq2.as_bytes(), kmer);
    let sparse_al = lcskpp(&matches, kmer);
    sparse_al.score as i32
}

#[allow(dead_code)]
pub fn block_align(seq1: &str, seq2: &str, match_score: i8, mismatch_score: i8, 
    gap_open_score: i8, gap_extend_score: i8) -> AlignResult {
    let min_block_size = 32;
    let max_block_size = 256;

    // A gap of length n will cost: open + extend * (n - 1)
    let gaps = Gaps { open: gap_open_score, extend: gap_extend_score };
    let nw: NucMatrix = NucMatrix::new_simple(match_score, mismatch_score);

    // Note that PaddedBytes, Block, and Cigar can be initialized with sequence length
    // and block size upper bounds and be reused later for shorter sequences, to avoid
    // repeated allocations.
    let r = PaddedBytes::from_bytes::<NucMatrix>(seq1.as_bytes(), max_block_size);
    let q = PaddedBytes::from_bytes::<NucMatrix>(seq2.as_bytes(), max_block_size);

    // Align with traceback, but no X-drop threshold (global alignment).
    let mut a = Block::<true, false>::new(q.len(), r.len(), max_block_size);
    
    a.align(&q, &r, &nw, gaps, min_block_size..=max_block_size, 0);
    let res = a.res();
    let mut cigar = Cigar::new(res.query_idx, res.reference_idx);
    a.trace().cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);
    res

}

fn ascii_to_decimals(ascii_string: &str) -> Vec<u8> {
    ascii_string.chars().map(|c| c as u8 - 33).collect()
}

#[allow(dead_code)]
pub fn merge_pe_by_qual(seq1: &str, seq2: &str, qual1: &str, qual2: &str, pos1: usize) -> (String, String) {
    let overlap_len = qual1[pos1..].len();
    let read1_overalp_qual = ascii_to_decimals(&qual1[pos1..]);
    let read2_overalp_qual = ascii_to_decimals(&qual2[..overlap_len]);
    // let read1_overalp_qual = qual1[pos1..].chars().map(|c| c.to_digit(10).unwrap() as i32).sum::<i32>();
    // let read2_overalp_qual = qual2[..overlap_len].chars().map(|c| c.to_digit(10).unwrap() as i32).sum::<i32>();
    let mut merge_seq = String::new();
    let mut merge_qual = String::new();
    if read1_overalp_qual > read2_overalp_qual {
        merge_seq.push_str(seq1);
        merge_seq.push_str(&seq2[overlap_len..]);
        merge_qual.push_str(qual1);
        merge_qual.push_str(&qual2[overlap_len..]);
    } else {
        merge_seq.push_str(&seq1[..pos1]);
        merge_seq.push_str(seq2);
        merge_qual.push_str(&qual1[..pos1]);
        merge_qual.push_str(qual2);
    }
    (merge_seq, merge_qual)
}

#[allow(dead_code)]
pub fn readfa_to_btreemap(fa_file: &str) -> BTreeMap<Vec<u8>, String> {
    let mut fa_map: BTreeMap<Vec<u8>, String> = BTreeMap::new();
    let fa_path = Path::new(&fa_file);
    if fa_path.exists() {
        let ref_reader = fasta::Reader::from_file(&fa_file).unwrap();
        
        fa_map.extend(ref_reader.records().map(|result| {
            let result_data = result.unwrap();
            let id = result_data.id().into();
            let seq = result_data.seq().to_owned();
            (seq, id)
        }));
    }
    fa_map
}

#[allow(dead_code)]
pub fn readfa_to_id_btreemap(fa_file: &str) -> BTreeMap<String, Vec<u8>> {
    let mut fa_map: BTreeMap<String, Vec<u8>> = BTreeMap::new();
    let fa_path = Path::new(&fa_file);
    if fa_path.exists() {
        let ref_reader = fasta::Reader::from_file(&fa_file).unwrap();
        
        fa_map.extend(ref_reader.records().map(|result| {
            let result_data = result.unwrap();
            let id = result_data.id().into();
            let seq = result_data.seq().to_owned();
            (id, seq)
        }));
    }
    fa_map
}

// 集中管理所有修剪和质控参数
#[derive(Clone, Debug)]
pub struct TrimConfig {
    pub qual_trim: u8,
    pub n_trim: bool,
    pub length_offset: usize,
    pub ascii_offset: u8,
    pub poly_trim: Option<usize>,
    pub trim_name: bool,
}

#[derive(Clone, Debug)]
pub struct AlignConfig {
    pub seed_len: usize,
    pub error_tolerance: u8,
    pub match_score: i8,
    pub error_score: i8,
    pub gap_open_score: i8,
    pub gap_extend_score: i8,
}

#[derive(Clone, Debug)]
pub struct InputConfig {
    pub five_art_fa: String,
    pub three_art_fa: String,
    pub five_idx_fa: String,
    pub three_idx_fa: String,
    pub idx_loc: u8,
    pub pe1_fastq: String,
    pub pe2_fastq: String,
    pub outdir: String,
    pub merge_pe: bool,
    pub num_threads: usize,
    pub batch: u32,
    pub train: usize,
}

pub trait LowQualityControlTrait {
    /// 对序列和质量值进行两端修剪，返回修剪后的序列、质量、左/右修剪长度
    /// 新增poly_trim参数：如Some(6)表示去除两端连续大于等于6bp的A/C/G/T
    fn trim_low_quality(
        &self,
        seq: &str,
        qual: &str,
        trim_config: &TrimConfig,
    ) -> Option<(String, String, usize, usize)>;
}

pub struct LowQualityControlTrimmer;

impl LowQualityControlTrait for LowQualityControlTrimmer {
    fn trim_low_quality(
        &self,
        seq: &str,
        qual: &str,
        trim_config: &TrimConfig,
    ) -> Option<(String, String, usize, usize)> {
        let qual_threshold = if trim_config.qual_trim > 0 { Some(trim_config.qual_trim) } else { None };
        let n_trim = trim_config.n_trim;
        let min_length = trim_config.length_offset;
        let ascii_offset = trim_config.ascii_offset;
        let poly_trim = trim_config.poly_trim;
        let mut left = 0;
        let mut right = seq.len();
        let seq_bytes = seq.as_bytes();
        let qual_bytes = qual.as_bytes();

        // 1. 5'端低质量修剪
        if let Some(thresh) = qual_threshold {
            let (mut l_q, mut r_q, mut l_i, mut r_i) = (0, 0, 0, 0);
            // 3'端
            while right > left && right >= seq.len() / 2 {
                let q = qual_bytes[right - 1] - ascii_offset;
                r_q += q;
                r_i += 1;
                if r_q < thresh * r_i {
                    right -= 1;
                } else {
                    break;
                }
            }
            // 5'端
            while (right - left) >= min_length - 1 {
                let q = qual_bytes[left] - ascii_offset;
                l_q += q;
                l_i += 1;
                if l_q < thresh * l_i {
                    left += 1;
                } else {
                    break;
                }
            }
            // 5'端修剪后，补充判断是否需要3'端修剪
            if (right - left) >= min_length && right <= seq.len() / 2 {
                while (right - left) >= min_length - 1 {
                    let q = qual_bytes[right - 1] - ascii_offset;
                    r_q += q;
                    r_i += 1;
                    if r_q < thresh * r_i {
                        right -= 1;
                    } else {
                        break;
                    }
                }
            }
        }

        // 2. N/非ACGT修剪
        if n_trim {
            // 3'端
            let mut r_n = false;
            let mut tmp_right = right;
            while tmp_right > left && tmp_right >= seq.len() / 2 {
                let c = seq_bytes[tmp_right - 1] as char;
                if "ACGTacgt".contains(c) && r_n {
                    break;
                } else if !"ACGTacgt".contains(c) {
                    tmp_right -= 1;
                    r_n = true;
                } else {
                    tmp_right -= 1;
                }
            }
            if r_n {
                right = tmp_right;
            }
            // 5'端
            let mut l_n = false;
            let mut tmp_left = left;
            while tmp_left < right && right - tmp_left >= min_length - 1 {
                let c = seq_bytes[tmp_left] as char;
                if "ACGTacgt".contains(c) && l_n {
                    break;
                } else if !"ACGTacgt".contains(c) {
                    tmp_left += 1;
                    l_n = true;
                } else {
                    tmp_left += 1;
                }
            }
            if l_n {
                left = tmp_left;
            }
            // 5'端修剪后，补充判断是否需要3'端修剪
            if (right - left) >= min_length && right <= seq.len() {
                let mut r_n = false;
                let mut tmp_right = right;
                while (tmp_right - left) >= min_length - 1 {
                    let c = seq_bytes[tmp_right - 1] as char;
                    if "ACGTacgt".contains(c) && r_n {
                        break;
                    } else if !"ACGTacgt".contains(c) {
                        tmp_right -= 1;
                        r_n = true;
                    } else {
                        tmp_right -= 1;
                    }
                }
                if r_n {
                    right = tmp_right;
                }
            }
        }

        // 3. poly-A/C/G/T末端修剪（可选）
        if let Some(poly_n) = poly_trim {
            // 3'端
            if right - left >= min_length {
                let base = seq_bytes[right - 1] as char;
                let mut count = 1;
                let mut temp_right = right;
                while temp_right > left && (temp_right - left) >= min_length {
                    let c = seq_bytes[temp_right - 1] as char;
                    if c == base && "ACGTacgt".contains(c) {
                        count += 1;
                        temp_right -= 1;
                    } else {
                        break;
                    }
                }
                if count >= poly_n {
                    right = temp_right;
                }
            }
            // 5'端
            if right - left >= min_length {
                let base = seq_bytes[left] as char;
                let mut count = 1;
                let mut temp_left: usize = left;
                while (right - temp_left) >= min_length {
                    let c = seq_bytes[temp_left] as char;
                    if c == base && "ACGTacgt".contains(c) {
                        count += 1;
                        temp_left += 1;
                    } else {
                        break;
                    }
                }
                if count >= poly_n {
                    left = temp_left;
                }
            }
        }

        // 4. 长度判断
        if right > left && (right - left) >= min_length {
            Some((
                seq[left..right].to_string(),
                qual[left..right].to_string(),
                left,
                seq.len() - right,
            ))
        } else {
            None
        }
    }
}
