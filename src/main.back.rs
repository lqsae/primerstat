use anyhow::{Context, Result};
use bio::alignment::pairwise::{Aligner, MatchParams};
use bio::alignment::AlignmentOperation;
use bio::alphabets::dna::revcomp;
use bio::io::fastq;
use chrono::Local;
use clap::Parser;
use flate2::read::GzDecoder;
use itertools::Itertools;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;
use rayon::prelude::*;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// 输入的fastq.gz文件
    #[arg(short, long)]
    input: String,

    /// 引物序列文件(TSV格式：name\tsequence)
    #[arg(short, long)]
    primers: String,

    /// 输出文件名
    #[arg(short, long)]
    output: Option<String>,

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

#[derive(Debug)]
struct AlignmentResult {
    edit_distance: i32,
    position: usize,
    query_aligned: String,
    target_aligned: String,
}

fn format_alignment_path(
    path: &[(usize, usize, AlignmentOperation)],
    query: &[u8],
    target: &[u8],
) -> (String, String) {
    let mut query_aligned = Vec::new();
    let mut target_aligned = Vec::new();

    for &(i, j, op) in path {
        if i >= query.len() || j >= target.len() {
            continue;
        }

        match op {
            AlignmentOperation::Match => {
                query_aligned.push(query[i]);
                target_aligned.push(target[j]);
            }
            AlignmentOperation::Subst => {
                query_aligned.push(query[i]);
                target_aligned.push(target[j]);
            }
            AlignmentOperation::Del => {
                query_aligned.push(query[i]);
                target_aligned.push(b'-');
            }
            AlignmentOperation::Ins => {
                query_aligned.push(b'-');
                target_aligned.push(target[j]);
            }
            AlignmentOperation::Xclip(_) | AlignmentOperation::Yclip(_) => {}
        }
    }

    (
        String::from_utf8_lossy(&query_aligned).into_owned(),
        String::from_utf8_lossy(&target_aligned).into_owned(),
    )
}


fn load_primers(primer_file: &str) -> Result<HashMap<String, (String, String)>> {
    let file = File::open(primer_file)?;
    let reader = BufReader::new(file);
    let mut primers = HashMap::new();

    for line in reader.lines() {
        let line = line.context("无法读取引物文件行，可能存在编码问题")?;
        // 跳过空行和注释行
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            // 清理可能存在的 BOM 和空白字符
            let name = parts[0].trim_start_matches('\u{feff}').trim();
            let seq = parts[1].trim();
            
            // 验证序列是否只包含有效的 DNA 字符
            if !seq.chars().all(|c| matches!(c, 'A' | 'T' | 'G' | 'C' | 'N')) {
                eprintln!("警告: 跳过包含无效字符的序列: {} - {}", name, seq);
                continue;
            }
            
            let seq = seq.to_uppercase();
            let rc_seq = String::from_utf8_lossy(&revcomp(seq.as_bytes())).into_owned();
            primers.insert(name.to_string(), (seq, rc_seq));
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
    thread_local! {
        static ALIGNER: std::cell::RefCell<Aligner<MatchParams>> = 
            std::cell::RefCell::new(Aligner::new(-1, -1, MatchParams::new(1, -1)));
    }
    
    ALIGNER.with(|aligner| {
        let mut aligner = aligner.borrow_mut();
        let alignment = aligner.semiglobal(primer.as_bytes(), seq.as_bytes());
        
        let edit_distance = (primer.len() as i32) - alignment.score;
        if edit_distance > max_errors {
            return None;
        }

        let (query_aligned, target_aligned) = format_alignment_path(
            &alignment.path().as_slice(),
            primer.as_bytes(),
            seq.as_bytes(),
        );
        
        Some(AlignmentResult {
            edit_distance,
            position: alignment.ystart,
            query_aligned,
            target_aligned,
        })
    })
}

fn create_primer_match(result: Option<AlignmentResult>) -> PrimerMatch {
    match result {
        Some(r) => {
            // 创建对齐的可视化表示
            let mut alignment_str = String::new();
            let mut match_line = String::new();
            
            // 遍历比对结果，创建匹配线
            for (q, t) in r.query_aligned.chars().zip(r.target_aligned.chars()) {
                if q == t {
                    match_line.push('|');  // 完全匹配
                } else if q == '-' || t == '-' {
                    match_line.push(' ');  // 插入或删除
                } else {
                    match_line.push('*');  // 错配
                }
            }
            
            // 组合三行对齐结果
            alignment_str.push_str(&r.query_aligned);
            alignment_str.push('|');  // 使用'|'作为分隔符
            alignment_str.push_str(&match_line);
            alignment_str.push('|');  // 使用'|'作为分隔符
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
    let mut best_score = i32::MAX;
    let mut best_result = None;

    // 直接遍历所有引物对
    for ((name1, (seq1, seq1_rc)), (name2, (seq2, seq2_rc))) in primers.iter().tuple_combinations() {
        let f_result_plus = align_sequence(seq1, seq, max_errors);
        let r_result_plus = align_sequence(seq2_rc, seq, max_errors);
        let f_result_minus = align_sequence(seq2, seq, max_errors);
        let r_result_minus = align_sequence(seq1_rc, seq, max_errors);

        let plus_score = match (f_result_plus.as_ref(), r_result_plus.as_ref()) {
            (Some(f), Some(r)) => f.edit_distance + r.edit_distance,
            _ => i32::MAX,
        };

        let minus_score = match (f_result_minus.as_ref(), r_result_minus.as_ref()) {
            (Some(f), Some(r)) => f.edit_distance + r.edit_distance,
            _ => i32::MAX,
        };

        if plus_score < best_score {
            best_score = plus_score;
            best_result = Some((
                name1.clone(),
                name2.clone(),
                '+',
                create_primer_match(f_result_plus),
                create_primer_match(r_result_plus),
                seq1.clone(),
                seq2_rc.clone(),
            ));
        }
        if minus_score < best_score {
            best_score = minus_score;
            best_result = Some((
                name2.clone(),
                name1.clone(),
                '-',
                create_primer_match(f_result_minus),
                create_primer_match(r_result_minus),
                seq2.clone(),
                seq1_rc.clone(),
            ));
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
        None => {
            Some(ReadAnalysis {
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

fn print_statistics(analyses: &[ReadAnalysis]) {
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

    // 计引物对使用情况
    let mut primer_pairs: HashMap<(String, String), usize> = HashMap::new();
    for analysis in analyses {
        let pair = (analysis.f_primer.clone(), analysis.r_primer.clone());
        *primer_pairs.entry(pair).or_insert(0) += 1;
    }

    println!("\n统计信息:");
    println!("总读数: {}", total);
    println!("成功找到两个引物的读数: {}", both_found);
    println!("成功率: {:.2}%", (both_found as f64 / total as f64) * 100.0);
    println!("正链数量: {}", plus_strand);
    println!("负链数量: {}", minus_strand);
    println!("二聚体数量: {} ({:.2}%)", 
             dimer_count, 
             (dimer_count as f64 / total as f64) * 100.0);

    println!("\n引物对使用统计:");
    let mut pairs: Vec<_> = primer_pairs.iter().collect();
    pairs.sort_by_key(|&(_, count)| std::cmp::Reverse(*count));
    
    for ((f, r), count) in pairs {
        println!("{} - {}: {} ({:.2}%)",
                f, r, count,
                (*count as f64 / total as f64) * 100.0);
    }
}

fn process_reads_parallel<R: BufRead>(
    reader: fastq::Reader<R>,
    primers: &HashMap<String, (String, String)>,
    max_errors: i32,
    min_distance: usize,
) -> Vec<ReadAnalysis> {
    let batch_size = 10000;
    let mut all_analyses = Vec::new();
    let mut batch = Vec::with_capacity(batch_size);
    
    for record in reader.records() {
        if let Ok(record) = record {
            batch.push(record);
            
            if batch.len() >= batch_size {
                let analyses: Vec<_> = batch.par_iter()
                    .filter_map(|record| analyze_read(record, primers, max_errors, min_distance))
                    .collect();
                all_analyses.extend(analyses);
                batch.clear();
            }
        }
    }
    
    if !batch.is_empty() {
        let analyses: Vec<_> = batch.par_iter()
            .filter_map(|record| analyze_read(record, primers, max_errors, min_distance))
            .collect();
        all_analyses.extend(analyses);
    }
    
    all_analyses
}

fn main() -> Result<()> {
    let start_time = Instant::now();
    let args = Args::parse();
    
    println!("正在加载引物文件...");
    let primers = load_primers(&args.primers)
        .context("Failed to load primers")?;
    println!("成功加载 {} 个引物", primers.len());
    
    let output_file = args.output.unwrap_or_else(|| {
        format!(
            "primer_analysis_{}.txt",
            Local::now().format("%Y%m%d_%H%M%S")
        )
    });

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
    
    println!("正在写入结果到文件: {}", output_file);
    write_analysis_results(&output_file, &analyses)
        .context("写入结果失败")?;

    println!("分析结果已保存到文件: {}", output_file);
    
    print_statistics(&analyses);

    Ok(())
}