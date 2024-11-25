// 移除不需要的导入
use anyhow::{Context, Result};
use clap::Parser;
use edlib_rs::edlibrs::{edlibAlignRs, EdlibAlignConfigRs, EdlibAlignModeRs, EDLIB_STATUS_OK, EdlibAlignTaskRs};
use bio::alphabets::dna::revcomp;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;
use serde::Serialize;
use std::path::PathBuf;
use std::env::args;
use rayon::prelude::*;
use std::sync::mpsc;
use std::thread;

// Args 结构保持不变
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// 输入的fastq.gz文件（单端测序）或者第一端序列文件（双端测序）
    #[arg(short, long)]
    input: String,

    /// 第二端序列文件（双端测序，可选）
    #[arg(short = '2', long)]
    input2: Option<String>,

    /// 引物序列文件(TSV格式：name\tsequence)
    #[arg(short, long)]
    primers: String,

    /// 输出目录
    #[arg(short = 'O', long, default_value = "output")]
    outdir: String,

    /// 样本名称
    #[arg(short = 'S', long)]
    sample: String,

    /// 最大允许错配数
    #[arg(short = 'e', long, default_value = "3")]
    max_errors: i32,

    /// 判定为二聚体的最小距离
    #[arg(short = 'd', long, default_value = "100")]
    min_distance: usize,

    /// 最大输出序列数（0表示输出所有序列）
    #[arg(short = 'n', long, default_value = "10000")]
    max_output: usize,

    /// 双端序列最小重叠长度
    #[arg(short = 'o', long, default_value = "10")]
    min_overlap: usize,

    /// 双端序列重叠区域最大错配率
    #[arg(short = 'm', long, default_value = "0.1")]
    max_mismatch_rate: f64,
}

// 其他结构体定义保持不变
#[derive(Debug, Clone)]
struct PrimerMatch {
    position: Option<usize>,
    errors: Option<usize>,
    found: bool,
}

#[derive(Debug)]
struct ReadAnalysis {
    read_id: String,
    length: usize,
    strand: char,
    f_primer: String,
    r_primer: String,
    f_match: PrimerMatch,
    r_match: PrimerMatch,
    distance: Option<usize>,
    is_dimer: bool,
}

#[derive(Debug, Clone)]
struct AlignmentResult {
    edit_distance: i32,
    position: usize,
}

#[derive(Debug, Clone)]
struct FastqRecord {
    id: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl FastqRecord {
    fn new() -> Self {
        FastqRecord {
            id: String::new(),
            seq: Vec::new(),
            qual: Vec::new(),
        }
    }
}


// FASTQ 解析器
struct FastqParser<R: BufRead> {
    reader: R,
    buffer: String,
}

impl<R: BufRead> FastqParser<R> {
    fn new(reader: R) -> Self {
        FastqParser {
            reader,
            buffer: String::with_capacity(1024),
        }
    }

    fn next_record(&mut self, record: &mut FastqRecord) -> Result<bool> {
        // 清空缓冲区和记录
        record.id.clear();
        record.seq.clear();
        record.qual.clear();
        
        // 读取 ID 行
        self.buffer.clear();
        match self.reader.read_line(&mut self.buffer) {
            Ok(0) => return Ok(false), // 文件结束
            Ok(_) => {
                if !self.buffer.starts_with('@') {
                    return Err(anyhow::anyhow!("FASTQ格式错误：ID行必须以@开头 - {}", self.buffer));
                }
                record.id = self.buffer[1..]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();
            },
            Err(e) => return Err(anyhow::anyhow!("读取ID行时发生错误: {}", e)),
        }

        // 读取序列行
        self.buffer.clear();
        match self.reader.read_line(&mut self.buffer) {
            Ok(0) => return Err(anyhow::anyhow!("FASTQ文件不完整：序列行意外结束")),
            Ok(_) => {
                let seq = self.buffer.trim();
                if seq.is_empty() {
                    return Err(anyhow::anyhow!("FASTQ格式错误：空序列行"));
                }
                record.seq = seq.as_bytes().to_vec();
            },
            Err(e) => return Err(anyhow::anyhow!("读取序列行时发生错误: {}", e)),
        }

        // 读取 + 行
        self.buffer.clear();
        match self.reader.read_line(&mut self.buffer) {
            Ok(0) => return Err(anyhow::anyhow!("FASTQ文件不完整：缺少+行")),
            Ok(_) => {
                if !self.buffer.starts_with('+') {
                    return Err(anyhow::anyhow!("FASTQ格式错误：第三行必须以+开头 - {}", self.buffer));
                }
            },
            Err(e) => return Err(anyhow::anyhow!("读取+行时发生错误: {}", e)),
        }

        // 读取质量行
        self.buffer.clear();
        match self.reader.read_line(&mut self.buffer) {
            Ok(0) => return Err(anyhow::anyhow!("FASTQ文件不完整：缺少质量行")),
            Ok(_) => {
                let qual = self.buffer.trim();
                if qual.is_empty() {
                    return Err(anyhow::anyhow!("FASTQ格式错误：空质量行"));
                }
                record.qual = qual.as_bytes().to_vec();
            },
            Err(e) => return Err(anyhow::anyhow!("读取质量行时发生错误: {}", e)),
        }

        // 检查序列长度和质量值长度是否匹配
        if record.seq.len() != record.qual.len() {
            return Err(anyhow::anyhow!(
                "FASTQ格式错误：序列长度({})与质量值长度({})不匹配",
                record.seq.len(),
                record.qual.len()
            ));
        }

        Ok(true)
    }
}

// 配对的FastqRecord结构
#[derive(Debug, Clone)]
struct PairedFastqRecord {
    r1: FastqRecord,
    merged: Option<FastqRecord>,
}


impl PairedFastqRecord {
    fn new(r1: FastqRecord) -> Self {
        PairedFastqRecord {
            r1,
            merged: None,
        }
    }
}

// 新增：提取序列ID的核心部分
fn get_sequence_id(full_id: &str) -> String {
    // 移除可能的方向标识（如 /1 或 /2）
    full_id.split_whitespace()
        .next()
        .unwrap_or(full_id)
        .trim_end_matches("/1")
        .trim_end_matches("/2")
        .to_string()
}


fn merge_paired_reads(
    r1: &FastqRecord,
    r2: &FastqRecord,
    min_overlap: usize,
    max_mismatch_rate: f64,
) -> Option<FastqRecord> {
    // 基本验证
    if r1.seq.is_empty() || r2.seq.is_empty() {
        return None;
    }

    // 将R2序列反向互补
    let r2_rc = revcomp(&r2.seq);
    let r2_rc_qual: Vec<u8> = r2.qual.iter().rev().copied().collect();

    // 寻找最佳重叠
    let max_overlap = r1.seq.len().min(r2_rc.len());
    let mut best_overlap_len = 0;
    let mut found_overlap = false;

    // 只有当最大可能重叠长度大于最小重叠要求时才尝试寻找重叠
    if max_overlap >= min_overlap {
        // 从最长的可能重叠开始检查
        for overlap_len in (min_overlap..=max_overlap).rev() {
            if r1.seq.len() < overlap_len || r2_rc.len() < overlap_len {
                continue;
            }

            let r1_part = &r1.seq[r1.seq.len() - overlap_len..];
            let r2_part = &r2_rc[..overlap_len];

            let mismatches = r1_part.iter()
                .zip(r2_part.iter())
                .filter(|(b1, b2)| b1 != b2)
                .count();

            let mismatch_rate = mismatches as f64 / overlap_len as f64;

            if mismatch_rate <= max_mismatch_rate {
                best_overlap_len = overlap_len;
                found_overlap = true;
                break;
            }
        }
    }

    // 构建合并序列
    let (merged_id, merged_seq, merged_qual) = if found_overlap {
        // 使用重叠区域合并，选择质量值较高的碱基
        let mut merged_seq = Vec::with_capacity(r1.seq.len() + r2_rc.len() - best_overlap_len);
        let mut merged_qual = Vec::with_capacity(r1.qual.len() + r2_rc_qual.len() - best_overlap_len);

        // 添加R1非重叠部分
        merged_seq.extend_from_slice(&r1.seq[..r1.seq.len() - best_overlap_len]);
        merged_qual.extend_from_slice(&r1.qual[..r1.qual.len() - best_overlap_len]);

        // 处理重叠区域，选择质量值较高的碱基
        let r1_overlap_seq = &r1.seq[r1.seq.len() - best_overlap_len..];
        let r1_overlap_qual = &r1.qual[r1.qual.len() - best_overlap_len..];
        let r2_overlap_seq = &r2_rc[..best_overlap_len];
        let r2_overlap_qual = &r2_rc_qual[..best_overlap_len];

        for i in 0..best_overlap_len {
            if r1_overlap_qual[i] >= r2_overlap_qual[i] {
                merged_seq.push(r1_overlap_seq[i]);
                merged_qual.push(r1_overlap_qual[i]);
            } else {
                merged_seq.push(r2_overlap_seq[i]);
                merged_qual.push(r2_overlap_qual[i]);
            }
        }

        // 添加R2剩余部分
        merged_seq.extend_from_slice(&r2_rc[best_overlap_len..]);
        merged_qual.extend_from_slice(&r2_rc_qual[best_overlap_len..]);

        (
            format!("{}_merged_overlap_{}", get_sequence_id(&r1.id), best_overlap_len),
            merged_seq,
            merged_qual
        )
    } else {
        // 直接连接序列
        (
            format!("{}_merged_concat", get_sequence_id(&r1.id)),
            [&r1.seq[..], &r2_rc[..]].concat(),
            [&r1.qual[..], &r2_rc_qual[..]].concat()
        )
    };

    // 验证合并结果
    if merged_seq.is_empty() || merged_qual.is_empty() || merged_seq.len() != merged_qual.len() {
        return None;
    }

    Some(FastqRecord {
        id: merged_id,
        seq: merged_seq,
        qual: merged_qual,
    })
}

// 统计相关的结构体
#[derive(Debug, Default)]
struct Statistics {
    total_reads: usize,
    both_primers_found: usize,
    plus_strand: usize,
    minus_strand: usize,
    dimer_count: usize,
    primer_pairs: HashMap<(String, String), usize>,
}

#[derive(Serialize)]
struct StatisticsOutput {
    sample_name: String,
    total_reads: usize,
    both_primers_found: usize,
    success_rate: f64,
    plus_strand: usize,
    minus_strand: usize,
    dimer_count: usize,
    dimer_rate: f64,
    primer_pairs: Vec<PrimerPairStat>,
}

#[derive(Serialize)]
struct PrimerPairStat {
    forward_primer: String,
    reverse_primer: String,
    count: usize,
    percentage: f64,
}

// 结果处理器
struct AnalysisWriter {
    writer: flate2::write::GzEncoder<std::io::BufWriter<File>>,
    sample_name: String,
    output_dir: PathBuf,
    count: usize,
    max_output: usize,
    stats: Statistics,
}

impl AnalysisWriter {
    fn new(output_file: &str, sample_name: &str, max_output: usize) -> Result<Self> {
        let output_path = PathBuf::from(output_file);
        let output_dir = output_path.parent()
            .ok_or_else(|| anyhow::anyhow!("无法获取输出目录"))?
            .to_path_buf();

        let file = File::create(output_file)?;
        let buf_writer = std::io::BufWriter::with_capacity(64 * 1024, file);
        let mut writer = flate2::write::GzEncoder::new(buf_writer, flate2::Compression::default());
        
        // 移除表头中的 F_Alignment 和 R_Alignment 列
        writeln!(
            writer,
            "Read_ID\tLength\tStrand\tF_Primer\tR_Primer\tF_Found\tF_Pos\tF_Errors\t\
             R_Found\tR_Pos\tR_Errors\tDistance\tIs_Dimer"
        )?;

        Ok(AnalysisWriter {
            writer,
            sample_name: sample_name.to_string(),
            output_dir,
            count: 0,
            max_output,
            stats: Statistics::default(),
        })
    }

    fn process(&mut self, analysis: &ReadAnalysis) -> Result<()> {
        // 更新统计信息部分保持不变
        self.stats.total_reads += 1;
        
        if analysis.f_match.found && analysis.r_match.found {
            self.stats.both_primers_found += 1;
        }
        
        match analysis.strand {
            '+' => self.stats.plus_strand += 1,
            '-' => self.stats.minus_strand += 1,
            _ => {},
        }
        
        if analysis.is_dimer {
            self.stats.dimer_count += 1;
        }
        
        let pair = (analysis.f_primer.clone(), analysis.r_primer.clone());
        *self.stats.primer_pairs.entry(pair).or_insert(0) += 1;

        // 如果达到最大输出数量，只收集统计信息不写入文件
        if self.max_output > 0 && self.count >= self.max_output {
            return Ok(());
        }
        
        // 写入分析结果，移除 f_alignment 和 r_alignment
        writeln!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            analysis.read_id,
            analysis.length,
            analysis.strand,
            analysis.f_primer,
            analysis.r_primer,
            analysis.f_match.found,
            analysis.f_match.position.map_or("-".to_string(), |p| p.to_string()),
            analysis.f_match.errors.map_or("-".to_string(), |e| e.to_string()),
            analysis.r_match.found,
            analysis.r_match.position.map_or("-".to_string(), |p| p.to_string()),
            analysis.r_match.errors.map_or("-".to_string(), |e| e.to_string()),
            analysis.distance.map_or("-".to_string(), |d| d.to_string()),
            analysis.is_dimer,
        )?;
        
        self.count += 1;
        Ok(())
    }

    fn get_statistics(&self) -> StatisticsOutput {
        StatisticsOutput {
            sample_name: self.sample_name.clone(),
            total_reads: self.stats.total_reads,
            both_primers_found: self.stats.both_primers_found,
            success_rate: if self.stats.total_reads > 0 {
                (self.stats.both_primers_found as f64 / self.stats.total_reads as f64) * 100.0
            } else {
                0.0
            },
            plus_strand: self.stats.plus_strand,
            minus_strand: self.stats.minus_strand,
            dimer_count: self.stats.dimer_count,
            dimer_rate: if self.stats.total_reads > 0 {
                (self.stats.dimer_count as f64 / self.stats.total_reads as f64) * 100.0
            } else {
                0.0
            },
            primer_pairs: self.stats.primer_pairs
                .iter()
                .map(|((f, r), count)| PrimerPairStat {
                    forward_primer: f.clone(),
                    reverse_primer: r.clone(),
                    count: *count,
                    percentage: if self.stats.total_reads > 0 {
                        (*count as f64 / self.stats.total_reads as f64) * 100.0
                    } else {
                        0.0
                    },
                })
                .collect(),
        }
    }

    fn finalize(&mut self) -> Result<()> {
        self.writer.try_finish()?;
        Ok(())
    }

    fn save_statistics(&self, stats: &StatisticsOutput) -> Result<()> {
        let stats_path = self.output_dir.join(format!("{}_statistics.json", self.sample_name));
        let stats_json = serde_json::to_string_pretty(stats)?;
        std::fs::write(&stats_path, stats_json)?;
        Ok(())
    }
}


fn load_primers(primer_file: &str) -> Result<HashMap<String, (String, String)>> {
    let file = File::open(primer_file)?;
    let reader = BufReader::new(file);
    let mut primers = HashMap::new();

    for line in reader.lines() {
        let line = line.context("无法读取引物文件行，可能存在编码问题")?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            let name = parts[0].trim_start_matches('\u{feff}').trim();
            let seq = parts[1].trim();
            
            if !seq.chars().all(|c| matches!(c.to_ascii_uppercase(), 'A' | 'T' | 'G' | 'C' | 'N')) {
                eprintln!("警告: 跳过包含无效字符的序列: {} - {}", name, seq);
                continue;
            }
            
            let seq = seq.to_uppercase();
            let rc_seq = String::from_utf8_lossy(&revcomp(seq.as_bytes())).into_owned();
            primers.insert(name.to_string(), (seq.to_uppercase(), rc_seq.to_uppercase()));
        } else {
            eprintln!("警告: 跳过格式不正确的行: {}", line);
        }
    }

    if primers.is_empty() {
        anyhow::bail!("未能加载任何有效的引物序列");
    }

    println!("成功加载 {} 个引物序列", primers.len());
    Ok(primers)
}



fn align_sequence(query: &[u8], target: &[u8], max_errors: i32) -> Option<AlignmentResult> {
    let additional_equalities: Vec<edlib_rs::edlibrs::EdlibEqualityPairRs> = Vec::new();
    let config = EdlibAlignConfigRs {
        k: max_errors,
        mode: EdlibAlignModeRs::EDLIB_MODE_HW,
        task: EdlibAlignTaskRs::EDLIB_TASK_PATH,
        additionalequalities: &additional_equalities,
    };

    let result = edlibAlignRs(query, target, &config);
    if result.status != EDLIB_STATUS_OK || result.editDistance > max_errors {
        return None;
    }

    // 获取最佳位置
    let end_locations = result.endLocations?;
    let start_locations = result.startLocations?;
    
    if end_locations.is_empty() || start_locations.is_empty() {
        return None;
    }

    let start_pos = start_locations[0] as usize;
    
    Some(AlignmentResult {
        edit_distance: result.editDistance,
        position: start_pos,
    })
}

fn create_primer_match(result: Option<AlignmentResult>) -> PrimerMatch {
    match result {
        Some(r) => PrimerMatch {
            position: Some(r.position),
            errors: Some(r.edit_distance as usize),
            found: true,
        },
        None => PrimerMatch {
            position: None,
            errors: None,
            found: false,
        },
    }
}


fn analyze_read(
    record: &FastqRecord,
    primers: &HashMap<String, (String, String)>,
    max_errors: i32,
    min_distance: usize,
) -> Option<ReadAnalysis> {
    let seq = &record.seq;
    let mut best_result = None;
    let mut best_score = i32::MAX;

    // 对引物对进行排序以确保一致的顺序
    let mut primer_pairs: Vec<_> = primers.iter().collect();
    primer_pairs.sort_by(|a, b| a.0.cmp(b.0));

    for (i, (name1, (seq1, seq1_rc))) in primer_pairs.iter().enumerate() {
        for (name2, (seq2, seq2_rc)) in primer_pairs[i+1..].iter() {
            // 正向链检查
            if let (Some(f), Some(r)) = (
                align_sequence(seq1.as_bytes(), seq, max_errors),
                align_sequence(seq2_rc.as_bytes(), seq, max_errors)
            ) {
                if f.position < r.position {
                    let score = f.edit_distance + r.edit_distance;
                    if score < best_score {
                        best_score = score;
                        best_result = Some((
                            name1.to_string(),  // 使用 to_string() 而不是 clone()
                            name2.to_string(),  // 使用 to_string() 而不是 clone()
                            '+',
                            create_primer_match(Some(f)),
                            create_primer_match(Some(r)),
                        ));
                    }
                }
            }

            // 反向链检查
            if let (Some(f), Some(r)) = (
                align_sequence(seq2.as_bytes(), seq, max_errors),
                align_sequence(seq1_rc.as_bytes(), seq, max_errors)
            ) {
                if f.position < r.position {
                    let score = f.edit_distance + r.edit_distance;
                    if score < best_score {
                        best_score = score;
                        best_result = Some((
                            name2.to_string(),  // 使用 to_string() 而不是 clone()
                            name1.to_string(),  // 使用 to_string() 而不是 clone()
                            '-',
                            create_primer_match(Some(f)),
                            create_primer_match(Some(r)),
                        ));
                    }
                }
            }
        }
    }


    let default_match = PrimerMatch {
        position: None,
        errors: None,
        found: false,
    };

    match best_result {
        Some((f_name, r_name, strand, f_match, r_match)) => {
            let distance = match (f_match.position, r_match.position) {
                (Some(f_pos), Some(r_pos)) => {
                    if r_pos > f_pos {
                        Some(r_pos - f_pos)
                    } else {
                        None
                    }
                },
                _ => None,
            };

            let is_dimer = match distance {
                Some(d) => d < min_distance,
                None => false,
            };

            Some(ReadAnalysis {
                read_id: record.id.clone(),
                length: record.seq.len(),
                strand,
                f_primer: f_name,
                r_primer: r_name,
                f_match,
                r_match,
                distance,
                is_dimer,
            })
        },
        None => Some(ReadAnalysis {
            read_id: record.id.clone(),
            length: record.seq.len(),
            strand: '?',
            f_primer: String::from("-"),
            r_primer: String::from("-"),
            f_match: default_match.clone(),
            r_match: default_match,
            distance: None,
            is_dimer: false,
        })
    }
}


// 添加一个用于并行处理的批次结构
#[derive(Default)]
struct ReadBatch {
    records: Vec<PairedFastqRecord>,
}

impl ReadBatch {
    fn new() -> Self {
        ReadBatch {
            records: Vec::with_capacity(1000), // 每批处理1000条记录
        }
    }

    fn is_full(&self) -> bool {
        self.records.len() >= 1000
    }
}


fn process_reads<R1: BufRead, R2: BufRead>(
    reader1: R1,
    reader2: Option<R2>,
    primers: &HashMap<String, (String, String)>,
    max_errors: i32,
    min_distance: usize,
    min_overlap: usize,
    max_mismatch_rate: f64,
    outdir: &str,
    sample: &str,
    max_output: usize,
) -> Result<()> {
    let mut parser1 = FastqParser::new(reader1);
    let mut parser2 = reader2.map(|r| FastqParser::new(r));
    
    let mut record_count = 0;
    let mut record1 = FastqRecord::new();
    let mut record2 = FastqRecord::new();

    // 创建结果写入器
    let result_file = PathBuf::from(outdir)
        .join(format!("{}_primer_analysis.txt.gz", sample));
    let mut writer = AnalysisWriter::new(
        result_file.to_str().unwrap(),
        sample,
        max_output
    )?;

    // 创建通道用于传输分析结果
    let (tx, rx) = mpsc::channel();

    // 启动写入线程
    let writer_thread = thread::spawn(move || {
        let mut writer = writer;
        for analysis in rx {
            if let Err(e) = writer.process(&analysis) {
                eprintln!("写入结果时发生错误: {}", e);
            }
        }
        if let Err(e) = writer.finalize() {
            eprintln!("完成写入时发生错误: {}", e);
        }
        writer
    });

    // 用于批处理的缓冲区
    let mut current_batch = ReadBatch::new();

    // 处理批次的闭包
    let process_batch = |batch: ReadBatch, 
                        primers: &HashMap<String, (String, String)>,
                        tx: &mpsc::Sender<ReadAnalysis>| {
        batch.records.par_iter().for_each(|record| {
            if let Some(analysis) = if let Some(ref merged) = record.merged {
                analyze_read(merged, primers, max_errors, min_distance)
            } else {
                analyze_read(&record.r1, primers, max_errors, min_distance)
            } {
                if let Err(e) = tx.send(analysis) {
                    eprintln!("发送分析结果时发生错误: {}", e);
                }
            }
        });
    };

    loop {
        match parser1.next_record(&mut record1) {
            Ok(true) => {
                let paired_record = if let Some(ref mut parser2) = parser2 {
                    match parser2.next_record(&mut record2) {
                        Ok(true) => {
                            let r1_id = get_sequence_id(&record1.id);
                            let r2_id = get_sequence_id(&record2.id);
                            
                            if r1_id != r2_id {
                                continue;
                            }

                            let mut paired = PairedFastqRecord::new(record1.clone());
                            paired.merged = merge_paired_reads(
                                &record1,
                                &record2,
                                min_overlap,
                                max_mismatch_rate
                            );
                            paired
                        },
                        Ok(false) => break,
                        Err(_) => continue,
                    }
                } else {
                    PairedFastqRecord::new(record1.clone())
                };
                
                record_count += 1;
                if record_count % 100_000 == 0 {
                    println!("已处理 {} 条序列", record_count);
                }

                // 将记录添加到当前批次
                current_batch.records.push(paired_record);

                // 如果批次已满，进行并行处理
                if current_batch.is_full() {
                    let batch = std::mem::replace(&mut current_batch, ReadBatch::new());
                    process_batch(batch, primers, &tx);
                }
            },
            Ok(false) => break,
            Err(_) => continue,
        }
    }

    // 处理最后一个批次
    if !current_batch.records.is_empty() {
        process_batch(current_batch, primers, &tx);
    }

    // 关闭发送端，让接收线程知道没有更多数据
    drop(tx);

    // 等待写入线程完成并获取writer
    let mut writer = writer_thread.join().unwrap();

    // 保存统计信息
    let statistics = writer.get_statistics();
    writer.save_statistics(&statistics)?;

    // 打印统计信息
    println!("总处理序列数: {}", record_count);
    println!("\n统计信息:");
    println!("样本名称: {}", statistics.sample_name);
    println!("总读数: {}", statistics.total_reads);
    println!("成功找到两个引物的读数: {}", statistics.both_primers_found);
    println!("成功率: {:.2}%", statistics.success_rate);
    println!("正链数量: {}", statistics.plus_strand);
    println!("负链数量: {}", statistics.minus_strand);
    println!("二聚体数量: {} ({:.2}%)", 
             statistics.dimer_count, 
             statistics.dimer_rate);

    println!("\n引物对使用统计:");
    for pair in &statistics.primer_pairs {
        println!("{} - {}: {} ({:.2}%)",
                pair.forward_primer,
                pair.reverse_primer,
                pair.count,
                pair.percentage);
    }

    Ok(())
}

// main 函数需要相应修改
fn main() -> Result<()> {
    let start_time = Instant::now();
    let args = Args::parse();
    
    // 创建输出目录
    std::fs::create_dir_all(&args.outdir)?;
    
    println!("正在加载引物文件...");
    let primers = load_primers(&args.primers)
        .context("加载引物文件失败")?;
    println!("成功加载 {} 个引物", primers.len());

    println!("正在读取FASTQ文件...");
    let file1 = File::open(&args.input)
        .context("无法打开R1文件")?;
        
    let reader1: Box<dyn BufRead + Send> = if args.input.ends_with(".gz") {
        Box::new(BufReader::with_capacity(
            8 * 1024 * 1024,
            flate2::read::MultiGzDecoder::new(file1)
        ))
    } else {
        Box::new(BufReader::with_capacity(8 * 1024 * 1024, file1))
    };

    let reader2 = if let Some(input2) = args.input2.as_ref() {
        println!("检测到双端测序数据，正在读取R2文件...");
        let file2 = File::open(input2)
            .context("无法打开R2文件")?;
        
        let reader: Box<dyn BufRead + Send> = if input2.ends_with(".gz") {
            Box::new(BufReader::with_capacity(
                8 * 1024 * 1024,
                flate2::read::MultiGzDecoder::new(file2)
            ))
        } else {
            Box::new(BufReader::with_capacity(8 * 1024 * 1024, file2))
        };
        Some(reader)
    } else {
        None
    };

    println!("开始分析序列...");
    process_reads(
        reader1,
        reader2,
        &primers,
        args.max_errors,
        args.min_distance,
        args.min_overlap,
        args.max_mismatch_rate,
        &args.outdir,
        &args.sample,
        args.max_output,  // 传递 max_output 参数
    ).context("序列分析过程中发生错误")?;
    
    let total_time = start_time.elapsed();
    println!("分析完成! 总运行时间: {:.2}s", total_time.as_secs_f64());

    Ok(())
}