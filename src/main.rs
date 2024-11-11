use anyhow::{Context, Result};
use chrono::Local;
use clap::Parser;
use edlib_rs::edlibrs::{edlibAlignRs, EdlibAlignConfigRs, EdlibAlignModeRs, EDLIB_STATUS_OK, EdlibAlignTaskRs};
use flate2::read::GzDecoder;
use itertools::Itertools;
use bio::alphabets::dna::revcomp;
use bio::io::fastq;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;
use rayon::prelude::*;
use bio::io::fastq::FastqRead;
use serde::Serialize;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// 输入的fastq.gz文件
    #[arg(short, long)]
    input: String,

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
}

#[derive(Debug, Clone)]
struct PrimerMatch {
    position: Option<usize>,
    errors: Option<usize>,
    alignment: String,
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
    query_aligned: String,
    target_aligned: String,
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


fn align_sequence(primer: &str, seq: &str, max_errors: i32) -> Option<AlignmentResult> {
    let mut config = EdlibAlignConfigRs::default();
    config.mode = EdlibAlignModeRs::EDLIB_MODE_HW;
    config.task = EdlibAlignTaskRs::EDLIB_TASK_PATH;
    config.k = max_errors;
    
    let result = edlibAlignRs(primer.as_bytes(), seq.as_bytes(), &config);

    if result.status != EDLIB_STATUS_OK || result.editDistance < 0 || result.editDistance > max_errors {
        return None;
    }
    
    let (start_pos, end_pos) = result.startLocations
        .zip(result.endLocations)
        .and_then(|(starts, ends)| {
            if starts.is_empty() || ends.is_empty() {
                None
            } else {
                Some((starts[0] as usize, ends[0] as usize))
            }
        })?;
    
    if start_pos >= seq.len() || end_pos >= seq.len() || end_pos < start_pos {
        return None;
    }
    
    let alignment = result.alignment.filter(|a| !a.is_empty())?;
    let aligned_region = &seq[start_pos..=end_pos];
    
    let (query_aligned, target_aligned) = format_alignment(primer, aligned_region, &alignment);
    
    Some(AlignmentResult {
        edit_distance: result.editDistance,
        position: start_pos,
        query_aligned,
        target_aligned,
    })
}



fn format_alignment(query: &str, target: &str, path: &[u8]) -> (String, String) {
    let capacity = path.len();
    let mut query_aligned = String::with_capacity(capacity);
    let mut target_aligned = String::with_capacity(capacity);
    let mut q_idx = 0;
    let mut t_idx = 0;

    let query_chars: Vec<char> = query.chars().collect();
    let target_chars: Vec<char> = target.chars().collect();

    for &op in path {
        match op {
            0 | 3 => { // 匹配
                if let (Some(&q), Some(&t)) = (query_chars.get(q_idx), target_chars.get(t_idx)) {
                    query_aligned.push(q);
                    target_aligned.push(t);
                    q_idx += 1;
                    t_idx += 1;
                }
            },
            1 => { // 插入
                if let Some(&t) = target_chars.get(t_idx) {
                    query_aligned.push('-');
                    target_aligned.push(t);
                    t_idx += 1;
                }
            },
            2 => { // 删除
                if let Some(&q) = query_chars.get(q_idx) {
                    query_aligned.push(q);
                    target_aligned.push('-');
                    q_idx += 1;
                }
            },
            _ => continue
        }
    }

    query_aligned.extend(std::iter::repeat('-').take(target_aligned.len().saturating_sub(query_aligned.len())));
    target_aligned.extend(std::iter::repeat('-').take(query_aligned.len().saturating_sub(target_aligned.len())));

    (query_aligned, target_aligned)
}


fn create_primer_match(result: Option<AlignmentResult>) -> PrimerMatch {
    match result {
        Some(r) => {
            let capacity = r.query_aligned.len() * 3 + 2; // 为三行加两个分隔符预分配空间
            let mut alignment_str = String::with_capacity(capacity);
            let match_line: String = r.query_aligned.chars()
                .zip(r.target_aligned.chars())
                .map(|(q, t)| {
                    if q == t { '|' }
                    else if q == '-' || t == '-' { ' ' }
                    else { '*' }
                })
                .collect();
            
            alignment_str.push_str(&r.query_aligned);
            alignment_str.push('|');
            alignment_str.push_str(&match_line);
            alignment_str.push('|');
            alignment_str.push_str(&r.target_aligned);

            PrimerMatch {
                position: Some(r.position),
                errors: Some(r.edit_distance as usize),
                alignment: alignment_str,
                found: true,
            }
        },
        None => PrimerMatch {
            position: None,
            errors: None,
            alignment: String::from("-"),
            found: false,
        },
    }
}

fn get_best_primer_pair(
    seq: &str,
    primers: &HashMap<String, (String, String)>,
    max_errors: i32,
) -> Option<(String, String, char, PrimerMatch, PrimerMatch, String, String)> {
    let mut best_result = None;
    let mut best_score = i32::MAX;

    for ((name1, (seq1, seq1_rc)), (name2, (seq2, seq2_rc))) in primers.iter().tuple_combinations() {
        // 检查正向链
        if let (Some(f), Some(r)) = (align_sequence(seq1, seq, max_errors), 
                                    align_sequence(seq2_rc, seq, max_errors)) {
            if f.position < r.position {
                let score = f.edit_distance + r.edit_distance;
                if score < best_score {
                    best_score = score;
                    best_result = Some((
                        name1.clone(),
                        name2.clone(),
                        '+',
                        create_primer_match(Some(f)),
                        create_primer_match(Some(r)),
                        seq1.clone(),
                        seq2_rc.clone()
                    ));
                }
            }
        }

        // 检查负向链
        if let (Some(f), Some(r)) = (align_sequence(seq2, seq, max_errors),
                                    align_sequence(seq1_rc, seq, max_errors)) {
            if f.position < r.position {
                let score = f.edit_distance + r.edit_distance;
                if score < best_score {
                    best_score = score;
                    best_result = Some((
                        name2.clone(),
                        name1.clone(),
                        '-',
                        create_primer_match(Some(f)),
                        create_primer_match(Some(r)),
                        seq2.clone(),
                        seq1_rc.clone()
                    ));
                }
            }
        }
    }

    best_result
}


fn analyze_read(
    record: &fastq::Record,
    primers: &HashMap<String, (String, String)>,
    max_errors: i32,
    min_distance: usize,
) -> Option<ReadAnalysis> {
    let seq = std::str::from_utf8(record.seq()).ok()?;
    let best_match = get_best_primer_pair(seq, primers, max_errors);

    let default_match = PrimerMatch {
        position: None,
        errors: None,
        alignment: String::from("-"),
        found: false,
    };

    match best_match {
        Some((f_name, r_name, strand, f_match, r_match, _f_seq, _r_seq)) => {
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
                read_id: record.id().to_string(),
                length: record.seq().len(),
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
            read_id: record.id().to_string(),
            length: record.seq().len(),
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

fn write_analysis_results(
    output_file: &str,
    analyses: &[ReadAnalysis],
) -> Result<()> {
    let mut file = File::create(output_file)?;
    
    writeln!(
        file,
        "Read_ID\tLength\tStrand\tF_Primer\tR_Primer\tF_Found\tF_Pos\tF_Errors\t\
         R_Found\tR_Pos\tR_Errors\tDistance\tIs_Dimer\tF_Alignment\tR_Alignment"
    )?;

    for analysis in analyses {
        let f_alignment = analysis.f_match.alignment.replace('\n', "|");
        let r_alignment = analysis.r_match.alignment.replace('\n', "|");

        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
            f_alignment,
            r_alignment,
        )?;
    }

    Ok(())
}

#[derive(Serialize)]
struct Statistics {
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

fn collect_statistics(analyses: &[ReadAnalysis], sample_name: &str) -> Statistics {
    let total = analyses.len();
    let both_found = analyses.iter()
        .filter(|a| a.f_match.found && a.r_match.found)
        .count();
    let plus_strand = analyses.iter()
        .filter(|a| a.strand == '+')
        .count();
    let minus_strand = analyses.iter()
        .filter(|a| a.strand == '-')
        .count();
    let dimer_count = analyses.iter()
        .filter(|a| a.is_dimer)
        .count();

    let mut primer_pairs: HashMap<(String, String), usize> = HashMap::new();
    for analysis in analyses {
        let pair = (analysis.f_primer.clone(), analysis.r_primer.clone());
        *primer_pairs.entry(pair).or_insert(0) += 1;
    }

    let mut pairs: Vec<PrimerPairStat> = primer_pairs.iter()
        .map(|((f, r), count)| PrimerPairStat {
            forward_primer: f.clone(),
            reverse_primer: r.clone(),
            count: *count,
            percentage: (*count as f64 / total as f64) * 100.0,
        })
        .collect();
    
    pairs.sort_by(|a, b| b.count.cmp(&a.count));

    Statistics {
        sample_name: sample_name.to_string(),
        total_reads: total,
        both_primers_found: both_found,
        success_rate: (both_found as f64 / total as f64) * 100.0,
        plus_strand,
        minus_strand,
        dimer_count,
        dimer_rate: (dimer_count as f64 / total as f64) * 100.0,
        primer_pairs: pairs,
    }
}

fn process_reads_parallel<R: BufRead + Send>(
    reader: fastq::Reader<R>,
    primers: &HashMap<String, (String, String)>,
    max_errors: i32,
    min_distance: usize,
) -> Vec<ReadAnalysis> {
    let batch_size = 100_000;
    let mut all_analyses = Vec::new();
    let mut records = Vec::with_capacity(batch_size);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_cpus::get())
        .build()
        .unwrap();

    // 使用互斥锁来安全地共享 reader
    let reader = std::sync::Mutex::new(reader);
    let mut record_count = 0;
    let mut record = fastq::Record::new();

    loop {
        // 清空批次缓冲区
        records.clear();
        
        // 读取一个批次的数据
        {
            let mut reader = reader.lock().unwrap();
            loop {
                match reader.read(&mut record) {
                    Ok(()) => {
                        if record.is_empty() {
                            break;  // 文件结束
                        }
                        records.push(record.clone());
                        record_count += 1;
                        if records.len() >= batch_size {
                            break;
                        }
                    },
                    Err(e) => {
                        eprintln!("读取记录时发生错误: {}", e);
                        break;
                    }
                }
            }
        }

        // 如果没有读取到任何记录，说明处理完成
        if records.is_empty() {
            break;
        }

        // 每处理100万条序列打印一次进度
        if record_count % 1_000_000 == 0 {
            println!("已处理 {} 条序列...", record_count);
        }

        // 并行处理当前批次
        let analyses: Vec<_> = pool.install(|| {
            records.par_iter()
                .with_max_len(1000)
                .filter_map(|record| analyze_read(record, primers, max_errors, min_distance))
                .collect()
        });

        // 添加到结果集
        all_analyses.extend(analyses);
    }

    println!("共处理 {} 条序列", record_count);
    all_analyses
}

fn main() -> Result<()> {
    let start_time = Instant::now();
    let args = Args::parse();
    
    // 创建输出目录
    std::fs::create_dir_all(&args.outdir)?;
    
    println!("正在加载引物文件...");
    let primers = load_primers(&args.primers)
        .context("Failed to load primers")?;
    println!("成功加载 {} 个引物", primers.len());
    
    // 构建输出文件路径
    let result_file = PathBuf::from(&args.outdir)
        .join(format!("{}_primer_analysis.txt", args.sample));
    let stats_file = PathBuf::from(&args.outdir)
        .join(format!("{}_statistics.json", args.sample));

    println!("正在读取FASTQ文件: {}", args.input);
    let file = File::open(&args.input)
        .context("无法打开输入文件")?;
    let gz_decoder = GzDecoder::new(file);
    let buf_reader = BufReader::new(gz_decoder);
    let reader = fastq::Reader::new(buf_reader);
    
    println!("开始并行分析序列...");
    let analyses = process_reads_parallel(reader, &primers, args.max_errors, args.min_distance);
    
    let total_time = start_time.elapsed();
    println!("分析完成! 总运行时间: {:.2}s", total_time.as_secs_f64());
    
    // 写入分析结果
    println!("正在写入结果到文件: {}", result_file.display());
    write_analysis_results(result_file.to_str().unwrap(), &analyses)
        .context("写入结果失败")?;

    // 收集并写入统计结果
    let statistics = collect_statistics(&analyses, &args.sample);
    println!("正在写入统计结果到文件: {}", stats_file.display());
    let stats_json = serde_json::to_string_pretty(&statistics)?;
    std::fs::write(&stats_file, stats_json)?;

    // 打印主要统计信息
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
