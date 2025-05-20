# SegPE

A simple program for classifing and separating paired-ends (PE) and single-ends (SE) FASTQ after removing adapters via adapter-index seq

Major algorithms:

- Exact match search: This is a more direct method for finding exact matches of artificial sequences.
Regular expression matching and Hamming distance:
This is suitable for detecting and locating index sequences, especially when mismatches of a certain length are taken into account.
- Process and classify PE/SE sequences: After removing the adapter and index sequences, classify the PE/SE sequences and create new PE/SE FASTQ files.
- Remove low-quality reads: This is a common step in bioinformatics to ensure that the data used for analysis is of high quality.
- Use Multi-threads, SIMD and AsyncIO to handle large amounts of data.

### Setup

```

# Install rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# or update
rustup update

cd SegPE
cargo build --release

## SIMD  edit Cargo.toml first
# Intel x86
cargo build --features simd_avx2 --release
# Arm 
cargo build --features simd_neon --release
```

## Example:

```
time RUST_LOG=INFO ./target/release/segpe --five-art-fa data/5_art.fa --three-art-fa data/3_art.fa --five-idx-fa data/idx.fa --pe1-fastq data/PE1.fastq.gz --pe2-fastq data/PE2.fastq.gz -n 1000 -b 100 -o data/output
```

### Version

0.1.16, build-240520

### Usage

[**飞书文档**](https://zimablueai.feishu.cn/wiki/MFyEw1nmAi6W6BkHJOBcBCP3nLc?from=from_copylink)

[**Notion中文文档**](https://past-midnight-b4b.notion.site/SegPE-299c19073d5c4f1b95452cbc04f7e650)

[**Notion Usage**](https://past-midnight-b4b.notion.site/SegPE-Usage-Documentation-675058e883314ea38ef6e9e90193a476)

### License

^Licensed under either of [**Apache License, Version 2.0**](./LICENSE) or [**MIT license**](./LICENSE) at your option.^
Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in this crate by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without any additional terms or conditions.
