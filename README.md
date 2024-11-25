# PrimerStat

PrimerStat 是一个高性能的 FASTQ 测序数据引物分析工具。它可以并行处理大规模测序数据，识别序列中的引物对，分析它们的位置、方向和匹配质量，并生成详细的统计报告。

## 功能特点

- 支持单端和双端测序数据分析
- 支持压缩的 FASTQ 文件（.fastq.gz）
- 多线程并行处理，提高分析速度
- 自动识别正向和反向链上的引物
- 检测引物二聚体
- 生成详细的分析报告和 JSON 格式的统计结果
- 支持自定义错配容忍度
- 提供引物对使用频率统计
- 双端序列智能合并功能

## 安装要求

- Rust 1.70 或更高版本
- 系统依赖：
  - zlib（用于处理 gzip 压缩文件）
  - clang/LLVM（用于 edlib-rs 编译）

## 安装

1. 克隆仓库
   ```bash
   git clone https://github.com/qsliu2017/primerstat.git
   ```

2. 进入项目目录
   ```bash
   cd primerstat
   ```

3. 编译项目
   ```bash
   cargo build --release
   ```

## 使用方法

### 基本用法
```bash
primerstat -i input.fastq.gz -p primers.tsv -S sample_name -O output_dir
```

### 完整参数说明
```
必需参数:
-i, --input <FILE>             输入的fastq.gz文件（单端测序）或者第一端序列文件（双端测序）
-p, --primers <FILE>           引物序列文件(TSV格式)
-S, --sample <NAME>            样本名称

可选参数:
-2, --input2 <FILE>            第二端序列文件（双端测序，可选）
-O, --outdir <DIR>             输出目录 [default: output]
-e, --max-errors <NUM>         最大允许错配数 [default: 3]
-d, --min-distance <NUM>       判定为二聚体的最小距离 [default: 100]
-n, --max-output <NUM>         详细结果文件最大输出序列数 [default: 10000]
-o, --min-overlap <NUM>        双端序列最小重叠长度 [default: 10]
-m, --max-mismatch-rate <NUM>  双端序列重叠区域最大错配率 [default: 0.1]
-h, --help                     显示帮助信息
-V, --version                  显示版本信息
```

### 引物文件格式

引物文件应为 TSV（制表符分隔）格式，包含两列：
```
Primer1 ATCGATCG
Primer2 GCTAGCTA
```
- 第一列：引物名称
- 第二列：引物序列（仅支持 A、T、G、C、N）

## 输出文件

### 1. 分析结果文件：`{sample}_primer_analysis.txt.gz`

压缩格式的制表符分隔文件，包含以下列：

| 字段名 | 说明 | 示例值 |
|--------|------|--------|
| Read_ID | 序列标识符 | @SRR1234567.1 |
| Length | 序列长度（bp） | 150 |
| Strand | 链方向（+/-/?） | + |
| F_Primer | 正向引物名称 | Primer1 |
| R_Primer | 反向引物名称 | Primer2 |
| F_Found | 是否找到正向引物 | true |
| F_Pos | 正向引物起始位置 | 0 |
| F_Errors | 正向引物错配数 | 1 |
| R_Found | 是否找到反向引物 | true |
| R_Pos | 反向引物起始位置 | 130 |
| R_Errors | 反向引物错配数 | 0 |
| Distance | 引物间距离 | 112 |
| Is_Dimer | 是否为二聚体 | false |

### 2. 统计结果文件：`{sample}_statistics.json`

包含以下主要信息：
- 总体统计
  - 样本名称
  - 总读数
  - 成功匹配率
  - 正/负链比例
  - 二聚体比例
- 引物对使用统计
  - 每对引物的使用次数
  - 使用频率百分比

## 双端测序数据处理

当提供双端测序数据时，程序会：

1. 自动识别配对的序列
2. 尝试合并双端读段：
   - 寻找最佳重叠区域
   - 根据质量值选择重叠区域的碱基
3. 合并策略：
   - 找到有效重叠：生成合并序列（ID标记为"merged_overlap_X"）
   - 无有效重叠：直接连接序列（ID标记为"merged_concat"）

### 双端数据参数优化建议

- 调整 min_overlap：
  - 增大可提高合并可靠性
  - 建议范围：10-30bp
- 调整 max_mismatch_rate：
  - 降低可提高准确性
  - 建议范围：0.05-0.15

## 使用示例

1. 单端数据分析：
```bash
primerstat -i input.fastq.gz -p primers.tsv -S sample01 -O results
```

2. 双端数据分析（严格参数）：
```bash
primerstat -i read1.fastq.gz -2 read2.fastq.gz -p primers.tsv -S sample01 \
          -O results -o 15 -m 0.05
```

3. 自定义参数分析：
```bash
primerstat -i input.fastq.gz -p primers.tsv -S sample01 -O results \
          -e 2 -d 150 -n 20000
```

## 性能优化

- 批处理序列读取
- Rayon 并行计算
- edlib-rs 高效比对
- 内存预分配
- 优化的字符串处理

## 注意事项

1. 输入要求：
   - FASTQ 文件需为 gzip 格式
   - 引物序列限 ATGCN
   - 引物文件需为 TSV 格式

2. 资源使用：
   - 推荐多核系统
   - 预留足够磁盘空间
   - 注意内存占用

## 作者

刘青山

## 许可证

MIT

## 更新日志

- v1.0.0 (2024-03-15)
  - 初始版本发布
  - 支持单端和双端测序数据
  - 实现基本分析功能