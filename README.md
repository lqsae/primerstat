# PrimerStat

PrimerStat 是一个高性能的 FASTQ 测序数据引物分析工具。它可以并行处理大规模测序数据，识别序列中的引物对，分析它们的位置、方向和匹配质量，并生成详细的统计报告。

## 功能特点

- 支持压缩的 FASTQ 文件（.fastq.gz）
- 多线程并行处理，提高分析速度
- 自动识别正向和反向链上的引物
- 检测引物二聚体
- 生成详细的分析报告和 JSON 格式的统计结果
- 支持自定义错配容忍度
- 提供引物对使用频率统计

## 安装要求

- Rust 1.70 或更高版本
- 系统依赖：
  - zlib（用于处理 gzip 压缩文件）
  - clang/LLVM（用于 edlib-rs 编译）

## 安装

1. 克隆仓库
   ```bash
   git clone https://github.com/yourusername/primerstat.git
   ```

2. 进入项目目录
   ```bash
   cd primerstat
   ```

3. 编译项目
   ```bash
   cargo build --release
   ```

## 使用

```
primerstat -i input.fastq.gz -p primers.tsv -S sample_name -O output_dir
```
### 完整参数说明   
```
Options:
-i, --input <FILE> 输入的 fastq.gz 文件
-p, --primers <FILE> 引物序列文件(TSV格式)
-O, --outdir <DIR> 输出目录 [default: output]
-S, --sample <NAME> 样本名称
-e, --max-errors <NUM> 最大允许错配数 [default: 3]
-d, --min-distance <NUM> 判定为二聚体的最小距离 [default: 100]
-h, --help 显示帮助信息
-V, --version 显示版本信息
```
### 引物文件格式

引物文件应为 TSV（制表符分隔）格式，包含两列：
Primer1 ATCGATCG
Primer2 GCTAGCTA
- 第一列：引物名称
- 第二列：引物序列

## 输出文件
PrimerStat 生成两种输出文件：分析结果文件和统计结果文件。它们分别提供了详细的序列级别信息和总体统计数据。

###1. `{sample}_primer_analysis.txt`：详细的分析结果，包含：
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
| F_Alignment | 正向引物比对可视化 | ATCG\|||||*|\|ATTG |
| R_Alignment | 反向引物比对可视化 | GCTA\|\|\|\|\|\|GCTA |

#### 比对可视化说明
- `|` 表示完全匹配
- `*` 表示错配
- ` ` (空格) 表示插入/删除
- `-` 表示空位
###2. `{sample}_statistics.json`：统计信息，包含：
   - 总读数
   - 成功匹配率
   - 正/负链比例
   - 二聚体比例
   - 引物对使用频率

   | 字段名 | 说明 | 单位/格式 |
   |--------|------|-----------|
   | sample_name | 样本名称 | 字符串 |
   | total_reads | 总读数 | 整数 |
   | both_primers_found | 找到两个引物的读数 | 整数 |
   | success_rate | 成功率 | 百分比 |
   | plus_strand | 正链数量 | 整数 |
   | minus_strand | 负链数量 | 整数 |
   | dimer_count | 二聚体数量 | 整数 |
   | dimer_rate | 二聚体比例 | 百分比 |
   | primer_pairs | 引物对使用统计 | 数组 |

#### primer_pairs 字段说明
- forward_primer: 正向引物名称
- reverse_primer: 反向引物名称
- count: 该引物对出现次数
- percentage: 该引物对占总读数的百分比

#### 数据解释说明

1. **链方向判定**
   - `+`: 正向链（正向引物在5'端）
   - `-`: 负向链（反向引物在5'端）
   - `?`: 未能判定方向

2. **二聚体判定**
   - 当两个引物之间的距离小于设定的最小距离（默认100bp）时判定为二聚体
   - Is_Dimer 字段为 true 表示检测到二聚体

3. **错配计算**
   - F_Errors 和 R_Errors 包含替换、插入和删除的总和
   - 当错配数超过设定阈值时，视为未找到引物

4. **特殊值说明**
   - 位置值为 "-" 表示未找到对应引物
   - 距离值为 "-" 表示无法计算（未找到一个或两个引物）

## 性能优化

- 使用批处理方式读取序列
- 采用 rayon 实现并行计算
- 使用 edlib-rs 进行高效序列比对
- 预分配内存减少重新分配
- 优化的字符串处理

## 注意事项

1. 输入文件必须是 gzip 压缩的 FASTQ 格式
2. 引物序列只能包含 ATCGN 碱基
3. 建议根据系统内存大小调整批处理大小
4. 处理大文件时请确保有足够的磁盘空间

## 示例

假设输入文件为 `input.fastq.gz`，引物文件为 `primers.tsv`，样本名称为 `sample_name`，输出目录为 `output`。运行命令：

```bash
primerstat -i input.fastq.gz -p primers.tsv -S sample_name -O output
```
## 许可证

[您的许可证类型]

## 作者

[作者信息]

## 贡献

欢迎提交问题和改进建议！

1. Fork 该项目
2. 创建您的特性分支
3. 提交您的更改
4. 推送到分支
5. 创建新的 Pull Request
