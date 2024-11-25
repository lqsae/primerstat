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
   git clone http://172.17.25.11:8081/cg_researchprojects/primerstat.git
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
-i, --input <FILE>             输入的fastq.gz文件（单端测序）或者第一端序列文件（双端测序）
-2, --input2 <FILE>            第二端序列文件（双端测序，可选）
-p, --primers <FILE>           引物序列文件(TSV格式)
-O, --outdir <DIR>             输出目录 [default: output]
-S, --sample <NAME>            样本名称
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
Primer1 ATCGATCG
Primer2 GCTAGCTA
- 第一列：引物名称
- 第二列：引物序列

## 输出文件
PrimerStat 生成两种输出文件：分析结果文件和统计结果文件。它们分别提供了详细的序列级别信息和总体统计数据。

### 1. `{sample}_primer_analysis.txt`：
- 当序列数超过 max_output_records 时，只输出前 N 条记录
- 文件末尾会添加截断说明
- 包含以下列：
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


#### 比对可视化说明
- `|` 表示完全匹配
- `*` 表示错配
- ` ` (空格) 表示插入/删除
- `-` 表示空位
### 2. `{sample}_statistics.json`：统计信息，包含：
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

1. 输入文件要求：
   - FASTQ 文件必须是 gzip 压缩格式
   - 引物序列只能包含 A、T、G、C、N 碱基
   - 引物文件必须是 TSV 格式

2. 输出控制：
   - 详细分析结果默认限制为 10000 条记录
   - 可通过 -n 参数调整输出记录数
   - 分析结果自动压缩为 gz 格式
   - 统计结果包含所有记录的信息

3. 资源使用：
   - 建议在多核系统上运行
   - 确保足够的磁盘空间
   - 监控内存使用情况


## 示例

假设输入文件为 `input.fastq.gz`，引物文件为 `primers.tsv`，样本名称为 `sample_name`，输出目录为 `output`。运行命令：

1. 基本使用：
```bash
primerstat -i input.fastq.gz -p primers.tsv -S sample01 -O results
```

2. 双端测序数据分析：
```bash
primerstat -i read1.fastq.gz -2 read2.fastq.gz -p primers.tsv -S sample01 -O results \
          -o 15 -m 0.05
```

3. 调整错配和二聚体参数：
```bash
primerstat -i input.fastq.gz -p primers.tsv -S sample01 -O results -e 2 -d 150 -n 20000
```

## 双端测序数据处理

当提供双端测序数据时，程序会：

1. 自动识别配对的序列
2. 尝试合并双端读段：
   - 寻找最佳重叠区域
   - 根据质量值选择重叠区域的碱基
   - 生成合并后的序列用于引物分析
3. 如果无法找到有效重叠：
   - 直接连接两端序列
   - 在结果中标记为 "merged_concat"
4. 合并成功的序列：
   - 在序列ID中标记重叠长度
   - 使用合并后的序列进行引物分析

### 双端数据合并参数

- `min_overlap`：要求的最小重叠长度
  - 默认值：10bp
  - 较大的值可以提高合并可靠性
  - 过大的值可能导致部分序列无法合并

- `max_mismatch_rate`：重叠区域允许的最大错配率
  - 默认值：0.1（10%）
  - 较小的值可以提高合并准确性
  - 过小的值可能导致部分序列无法合并


## 作者

刘青山

## 贡献

欢迎提交问题和改进建议！

1. Fork 该项目
2. 创建您的特性分支
3. 提交您的更改
4. 推送到分支
5. 创建新的 Pull Request

# 更新日志

- 2024-11-09 初始版本   
- 2024-11-15 添加了新的 max_output_records 参数说明，增加输出控制，优化了输出文件格式
- 2024-11-20 添加双端测序支持，新增序列合并功能和相关参数