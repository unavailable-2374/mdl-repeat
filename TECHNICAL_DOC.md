# mdl-repeat 技术手册

> 一份面向**使用者**的完整参考手册:讲清这个工具是什么、算法怎么跑、每个阶段的逻辑、全部参数的含义与默认值、输入输出格式,以及调参与边界。
>
> 本手册与 `src/` 当前代码、以及项目内权威的 `CLAUDE.md` 对齐。凡本手册与旧版 `README.md` 冲突之处,以本手册与代码为准。基准测试数字、历史实验与被回退的尝试记录在 `FINAL_REPORT.md`,本手册不重复。

---

## 目录

1. [工具是什么](#1-工具是什么)
2. [构建、运行与测试](#2-构建运行与测试)
3. [快速上手](#3-快速上手)
4. [理论基础:MDL 目标函数](#4-理论基础mdl-目标函数)
5. [流水线总览](#5-流水线总览)
6. [发现引擎](#6-发现引擎)
7. [K-mer 表与位置索引](#7-k-mer-表与位置索引)
8. [精修阶段](#8-精修阶段)
9. [MDL 选择](#9-mdl-选择)
10. [大基因组扩展](#10-大基因组扩展)
11. [Recall rescue(召回补救)](#11-recall-rescue召回补救)
12. [外部工具现状](#12-外部工具现状)
13. [输入与输出格式](#13-输入与输出格式)
14. [参数完整参考](#14-参数完整参考)
15. [调参指南与常见场景](#15-调参指南与常见场景)(含空输出 / 边界场景)
16. [核心数据结构](#16-核心数据结构)
17. [不变量与坑(改代码前必读)](#17-不变量与坑改代码前必读)
18. [源码地图](#18-源码地图)
19. [参考文献](#19-参考文献)

---

## 1. 工具是什么

**mdl-repeat 是一个 de novo 的重复序列 / 转座元件(TE)家族库构建器**,输入是一条基因组 FASTA、不需要任何先验重复库,输出是一组共识序列(每个重复家族一条)外加每个家族在基因组上的实例位置。

核心思想是把两件事结合起来:

- **RepeatScout 风格的 seed-and-extend 发现**——从高频 l-mer 种子出发,把一个家族的所有拷贝同时对齐,逐列长出共识序列。
- **最小描述长度(Minimum Description Length, MDL)模型选择**——一个家族只有在「把它的拷贝编码成相对共识的若干 edit」比「把这些拷贝逐字面编码」更省 bit 时才被保留。

**关键点:是 MDL 而不是人工阈值决定保留哪些家族、保留多少个。** 工具里没有一个需要你手调的 recall / precision 阈值——编码方案本身就惩罚了短家族 / 低拷贝家族,奖励了长家族 / 高拷贝家族。

### 1.1 适用范围(scope)

**在范围内:** 高拷贝、低分歧的散在重复——LTR 反转录转座子、LINE、SINE、MITE、DNA 转座子(含 Helitron,但见下方说明)。一个家族大致需要 **≥3 个拷贝**、且共识有一定长度,才能越过 MDL 门槛。

这里「低分歧」指拷贝-vs-共识的**替换率大致 ≤ 0.30**(默认 `-max-divergence`),覆盖近期到中等活跃的家族;古老、高度分歧(如 >40%)的家族被刻意排除——它们更适合 RECON 这类自比对工具。老 TE 富集的基因组可把 `-max-divergence` 调到 0.40,但要预期更多噪声实例(见 [§15.1](#151-我想提高灵敏度漏家族太多))。

> **divergence 的定义:** 全工具中的 divergence = 对齐区内的**替换率(仅错配,不计 indel)**(`align.c` DESIGN NOTE §C3)。这与多数 TE 工具「拷贝 vs 共识」用替换率的惯例一致;与计入 indel 的工具(如 RepeatMasker `.align` 的 `div`)相比,本工具的值会略低(通常 <3% 相对差,偏向接受)。
>
> **Helitron 等的特殊性:** seed-and-extend 对「内部区域保守、适合 k-mer 打种子」的元件(LTR/LINE)效果最好。Helitron(滚环复制、无 TIR、末端 motif 短而内部多变)以及有内部缺失的长元件,常被发现成**多个短共识片段**——这是 [§8.3](#83-fragment-assemblyrefine_assemble_fragments) 片段拼接存在的生物学动因;对老化 / 降解的 Helitron 家族报出碎片化共识是预期行为,不是 bug,且未对 Helitron 专门做基准。

**设计上不在范围内(这是设计取舍,不是 bug):**

- 1–2 拷贝的碎片(那是 RECON 的活)。
- 串联重复 / 卫星阵列 / 着丝粒重复(那是 TRF 的活,本工具用 tandem-distance 与周期性过滤主动排除它们)。
- 低复杂度序列(用 Shannon 熵过滤掉)。

> **怎么评估结果好坏:** 主指标是**家族级别 recall**(cd-hit + BLAST,RepeatModeler2 / EDTA 风格)。per-instance 的碱基级 recall 是个噪声较大的代理指标,会低估库的质量,不要拿它当主数字。
>
> **没有 curated 参考库的非模式基因组怎么办:** 家族级 recall 需要一个 ground-truth 库(如 TAIR10),非模式物种没有。实务替代:与在同一基因组上跑的 RepeatModeler2 库比对(高拷贝家族应大量重叠)、用 TE 蛋白库 BLASTX 共识库看已知 TE 结构域是否被覆盖、或人工检视最大的若干家族。

### 1.2 实现概况

约 12,000 行 C11,**构建期不依赖任何外部库**(只用 `-lm -pthread`)。唯一可选的**运行期**外部工具是 `seqkit`(仅用于可选 QC,见 [§12](#12-外部工具现状))。

---

## 2. 构建、运行与测试

### 2.1 构建

```bash
make                 # 产出 bin/mdl-repeat  (gcc -O3 -march=native -flto)
make PORTABLE=1      # 去掉 -march=native(可移植二进制,跨机分发用)
make clean
```

系统要求:GCC(C11),`libm` 与 pthread。无其它构建期依赖。

### 2.2 运行

```bash
# 最小运行
bin/mdl-repeat -sequence genome.fa -output families.fa -v

# 完整输出 + 多线程
bin/mdl-repeat -sequence genome.fa -output families.fa \
    -instances inst.bed -stats stats.tsv -threads 4 -vv

# 不带任何参数 → 打印权威、与代码一致的 flag 列表
bin/mdl-repeat
```

> **`-sequence` 与 `-output` 是仅有的两个必填参数。** 其余全部可选,有合理默认值。

### 2.3 测试

```bash
bash tests/run_tests.sh                       # 端到端合成基因组测试
gcc -O2 -std=c11 -I src -o /tmp/t tests/test_mdl.c       src/mdl.c    -lm && /tmp/t
gcc -O2 -std=c11 -I src -o /tmp/s tests/test_sweepline.c src/refine.c src/mdl.c -lm && /tmp/s
```

CI(`.github/workflows/ci.yml`)在每次 push 上跑:native 构建、portable 构建 + 测试套件、以及一个小用例的 valgrind。

---

## 3. 快速上手

最常见的用法只需要给输入和输出:

```bash
bin/mdl-repeat -sequence genome.fa -output library.fa -threads 8 -v
```

跑完你会得到 `library.fa`——每条记录是一个被 MDL 接受的重复家族共识。想同时拿到实例坐标和逐家族统计:

```bash
bin/mdl-repeat -sequence genome.fa -output library.fa \
    -instances library.bed -stats library.tsv -threads 8 -vv
```

几条经验法则:

- **`-threads`** 给到物理核数。它加速 l-mer 计数、精修、以及大基因组的分块发现;发现主循环本身在小基因组上是串行的。
- **大基因组(>200 Mb / >1 Gb)无需特殊操作**,采样与分块会自动触发(见 [§10](#10-大基因组扩展));要可复现请固定 `-seed`。
- **想提高灵敏度**(尤其高分歧 / 低拷贝家族遗漏多时)→ 试 `-recall-rescue`(见 [§11](#11-recall-rescue召回补救))。
- **`-vv`** 打印每个阶段的家族数变化,排查「家族在哪一步被丢掉」非常有用;`-trace-dir <dir>` 还会 dump 每个精修阶段后的家族集。

---

## 4. 理论基础:MDL 目标函数

MDL 把「描述这条基因组」的总代价最小化,代价分两部分——描述库本身,加上「给定库之后再描述基因组」:

```
DL(Genome) = DL(Library) + DL(Genome | Library)
```

所有代价单位都是 **bit**。下面是代码里真正实现的几块(`mdl.c`)。

### 4.1 Rissanen 通用整数码 `L_int(n)`(`mdl.c:18`)

把一个正整数编码成 `log2*(n) + log2(c0)` bit,其中 `log2*` 是迭代对数之和(反复 `x = log2(x)` 累加,直到 `x ≤ 0`),常数 `LOG2_C0 = 1.5179605508986484`(= log2(2.865…),保证满足 Kraft 不等式)。

> **关键设计:** `n ≤ 0` 故意返回 `INFINITY`,这样畸形输入会让对应家族**因代价无穷大而被拒**,而不是被算成「免费」。

参考值:`L_int(1)≈1.52`,`L_int(2)≈2.52`,`L_int(3)≈3.77`,`L_int(10)≈7.36`,`L_int(100)≈12.34`。

### 4.2 库代价(每条共识)`mdl_model_cost`(`mdl.c:130`)

```
model_cost = L_int(len) + 2·len
```

即:一个长度头(用 `L_int`)+ 每碱基 2 bit 的字面共识(均匀碱基分布,刻意保守)。

### 4.3 每实例代价 `mdl_instance_cost_full`(`mdl.c:73`)

```
C_instance = L_int(a) + L_int(m+1) + m·log2(3) + [位置项]
```

- `a` = 对齐长度,`m` = edit 数(`+1` 避免 `L_int(0)`)。
- `m·log2(3)`:每个被编辑位点指定 3 种替换碱基之一。
- **位置项**取决于 `-mdl-mode`:
  - `none` —— 无位置项。
  - `exact`(默认)—— `log2 C(a, m)`,用 `lgamma` 算精确二项系数,即「在 `a` 个位置中放 `m` 个 edit」的精确方案数。
  - `upper` —— `m·log2(a)`,是 `exact` 那个精确二项数 `log2 C(a,m)` 的**上界**(对典型 `m ≪ a` 通常更大),即更保守的编码(见 [§15.4](#154-mdl-编码模式怎么选))。

### 4.4 每家族打分 `mdl_score_family`(`mdl.c:143`)

```
total_savings = Σ_i (2·a_i − C_instance,i)
mdl_score     = total_savings − model_cost
```

一个家族「值得保留」当且仅当:它的拷贝压缩出来的收益(正 savings)超过存储其共识的代价。edit 数 `m_i` 来自**发现 / 对齐阶段**的 DP(`1 / −1 / −5` 打分),**不是**精修 DP;若 `num_edits` 为负则回退到 `divergence · a_i`。

### 4.5 ⚠️ 载重级注意:代价与 R 无关

**当前的每实例公式刻意丢掉了「家族数开销」**(type bit、`log2(R)` 家族编号、strand bit)。`consensus_length` 与 `num_families` 虽然作为参数传入,但被 `(void)` 忽略(`mdl.c:73`)。

**后果:一个家族的分数不随被接受家族数 R 变化。** 这正是为什么旧版里那个「迭代 R 到收敛」的循环和「prune 后恢复 pass」都被移除了——在 R 无关的代价下它们是无效的死代码。

> 注意:`R`(`mdl_select_library` 里的 `R_estimate`,带 `≥2` 下限)**仍作为 API 参数贯穿**,但只是为「将来重新引入 R 相关项」保留的脚手架,当前对分数零影响——删它会破坏接口,而不是简化逻辑。

> **改代码警告:** 不要在没有先重新引入 R 相关代价项的前提下,加入任何「假设分数会随家族数变化」的控制流。

---

## 5. 流水线总览

下面是 `main.c` 里的**确切当前执行顺序**:

```
genome_load                                              (1)
  └─ [若 raw_length > sample-size] genome_sample_windows (1b,超大基因组)
default_k  → l-mer 长度
discover_families  或  discover_chunked                  (2)  按 raw_length > chunk-size 分支
  └─ [若 -recall-rescue] recall_rescue_run               (2b,可选补救 pass)
kmer_count → kmer_trim → kmer_build_positions            (3)  k-mer 表 + 位置索引
compact(丢弃 <2 实例或共识过短的家族)                    (4)
refine_merge_families        (80-80-80 + union-find)     (5)
refine_split_families        (Otsu 双峰)                 (6)
refine_assemble_fragments    (空间共现)                   (6b)
mdl_select_library           (单遍 MDL 选择)              (7)
refine_prune_families        (排他覆盖规则)               (8)
  └─ [若 coalesce-factor > 0] refine_coalesce_tandem…   (8b,仅影响报告)
output_fasta / output_bed / output_stats                 (9)  若采样过则重映射坐标
  └─ [若 -external-qc <file>] external_qc_run_seqkit…    (9b,可选 QC)
```

几个对理解逻辑很关键的顺序点:

- **MDL 选择在 merge/split/assembly 之后跑**——它打分的是精修后的家族集,不是原始发现输出。
- **Prune 在 MDL 选择之后跑**,用的是它**自己的排他覆盖规则**(不是 MDL 分数)来丢边缘家族。
- **Tandem coalesce 与外部 QC 都在决策之后**——它们只改变报告出来的实例 / 多写一个 QC 文件,**不回馈**到选择里。
- 若 compaction 后 0 家族存活,提前空输出退出。
- `-trace-dir` 会 dump 阶段 01–08 后的家族集供调试。

---

## 6. 发现引擎

文件:`discover.c`、`align.c`、`kmer.c`、`discover_mask.c`。在基因组上做 seed-and-extend。基因组加载时**前端 pad 了 `PADLENGTH = 11000`**,使得向两侧延伸永远不会越出数组边界。

### 6.1 L-mer 频率表

用**对称 / canonical 哈希**统计每个 l-mer(一个正向 l-mer 与它的反向互补共用一个桶;`hash_function`,`discover.c:40`)。插入时有三道过滤:

- **Tandem-distance 过滤**——一次命中只有在距离上一次同链命中 ≥ `TANDEMDIST`(500 bp)时才计入频率,这样串联阵列不会虚高计数。
- **熵过滤 + 周期性过滤**——Shannon 熵拒绝(`MAXENTROPY −0.70`)以及周期扫描(`PERIODIC_MATCH_PCT 85`),丢掉低复杂度与短周期种子。
- **自动 l-mer 长度** `l = ceil(1 + log4(N))`(`discover.c:1496`)。哈希表大小是**动态**的:`max(16000057, 4N/l)` 取奇数(16M 只是地板值,不是固定质数)。31-base 的打包上限只存在于并行计数路径(`kmer.c`)。

### 6.2 种子选择(`find_besttmp`,`discover.c:594`)

贪心地取还没被 mask 的最高频 l-mer;有个 locality 快捷路径,会先复用上一次的 best hash,失败再回退到全表扫描。

### 6.3 N 序列同时带状延伸(`extend_right` / `extend_left`)

这是核心算法:把一个种子的**全部 N 个出现位置同时对齐**,一列一列地向两个方向长出共识。对每个候选碱基 ∈ {A,C,G,T},引擎以该碱基为约束重跑带状 DP,把所有出现位置上的最佳分数加总;让总和最大的那个碱基成为这一列的共识。

- 带半宽 `MAXOFFSET = 5`;打分 `match 1 / mismatch −1 / gap −5`,每个出现位置有个 `CAPPENALTY −20` 的下限。
- 当某列不能把总分提升 ≥ `MINIMPROVEMENT`(3),且这种「静默」持续了一个**自适应**窗口 `max(WHEN_TO_STOP=100, extended/10)` 时停止。
- 只有共识 ≥ `GOODLENGTH`(30 bp)才保留。

### 6.4 Masking(`mask_headptr`,`discover_mask.c`)

新家族发现后,用 1-vs-1 带状 DP 把它的出现位置重新定位回基因组,并把这些位置标记为「已占用」(`MAX_FAMILY_CLAIMS = 1`;被 tandem 剪掉的位置写 `CLAIM_PERMANENT = 255` 哨兵),这样下一轮种子迭代不会重新发现同一个家族。然后回到 [§6.2](#62-种子选择find_besttmpdiscoverc594)。

### 6.5 对齐期实例招募器(`align.c`)

这是与发现独立的一套机制:用带位置索引的 k-mer 表来招募 / 重对齐实例——`seed_genome_scan` 收集共识 k-mer 命中,`cluster_seed_hits` 把它们按 locus 聚成锚点(聚类窗口 `1.5·cons_len`),然后 `align_banded` 逐个精对齐。链向由 `cons_is_rc == genome_is_rc` 决定。

> **实例数上限(会影响 `copies=` 与 savings,务必知道):** 招募带两道上限——原始 seed hits 上限 `-max-seed-hits`(默认 50000),聚类后**每家族实例上限** `-max-instances`(默认 10000)。因 MDL savings 随拷贝数增长,拷贝数超过上限的**超高拷贝家族**(如 >10k 拷贝的 SINE/MITE,正是 [§1.1](#11-适用范围scope) 声称的 in-scope)会被截断,导致报告的 `copies=` 与 MDL savings 被**低估**(共识本身不受影响,它在发现阶段已建好)。**当任何家族触顶时程序会在 stderr 打 WARNING**,提示调大 `-max-instances`(`-max-seed-hits` 会自动 ≥ 它)。注意 `-maxoccurrences`(默认 10000)是**发现期**的每种子出现数上限,与招募期的 `-max-instances` 是两个独立的量。

### 6.6 并行

- l-mer 计数用 striped-lock 哈希(`NUM_STRIPES = 4096`),每线程独立池最后合并(`kmer.c`)。
- 位置索引并行计数、但按基因组顺序填充以保证确定性,受 `KMER_MAX_POSITIONS = 50000` 截断上限约束。
- **发现主循环本身是串行的;大基因组靠跨 chunk 并行**(见 [§10](#10-大基因组扩展))。

---

## 7. K-mer 表与位置索引

发现完成后,为精修流水线单独建一个 canonical k-mer 表(`kmer.c`)。

- **Canonical k-mer:** 打包进 64-bit 整数(每碱基 2 bit,最大 k = 31),取 `min(packed, revcomp(packed))` 为 canonical 形式;Fibonacci 哈希分布。
- **Trim:** `kmer_trim` 移除频率低于 `MINTHRESH`(默认 2)的 k-mer,减小内存、加速查找。
- **位置索引:** `kmer_build_positions` 为每个 k-mer 建位置数组,O(1) 查到某 k-mer 在基因组的所有出现位置;并行计数 + 基因组顺序填充保证确定性,单 k-mer 位置数截断在 `KMER_MAX_POSITIONS = 50000`。正负号编码链向(正=正链,负=反向互补)。

---

## 8. 精修阶段

文件:`refine.c`。在发现与 MDL 选择之间有四个变换,外加一个选择后 prune 和一个仅影响报告的 coalesce。**每个会改变家族内容的变换都有一道局部 MDL 接受检查**——保证该变换在其自身评分下不劣化(split 用的是放宽门槛)。注意这是**逐变换的局部门**,发生在 `mdl_select_library` 与 25%-prune 之前,**并不等于**对最终被选库总 DL 的单调改善保证。

### 8.1 Merge —— 80-80-80(`refine_merge_families`)

- **预筛:** k-mer Jaccard(`REFINE_SCREEN_K 8`,`MIN_JACCARD 0.15`),再双链 semi-global DP。
- **合并条件:** identity ≥ **0.80**、coverage ≥ **0.80**、aligned ≥ **80 bp**(放宽的 0.70/0.70 档只有在「实例重叠」确认时才允许)。
- **护栏:** length-ratio ≥ 0.7;**nested-element 否决**(当一条共识 ≥3× 另一条、且短的那条有 ≥3 拷贝时,要求 ≥50% 包含);**MDL 否决**(`estimate_merge_score ≤ 0 → 跳过`)。
- 传递性合并用 union-find 解决。共识-vs-共识 DP 打分(`2 / −3 / −2`)只产出 identity/coverage,**绝不**喂给 MDL。

### 8.2 Split —— Otsu 双峰(`refine_split_families`)

- 每个家族建 100-bin 分歧度直方图,找 Otsu 阈值,当双峰性(组间 / 总方差)≥ **0.20** 时拆分(0.20–0.40 的边界带由谷深检查兜底)。
- 需要 ≥3 实例、每个子簇 ≥3、最小分歧度间隔 0.03。
- 只有当两个子家族的 **MDL** 分数越过一个放宽门槛(原家族为正时,子家族 ≥ 0)才接受。

> **限制:** Otsu 每个家族只产一个阈值、即**最多拆成两组**。有 3+ 个真实亚家族的家族(如灵长类 AluJ/S/Y、或多次转座爆发的植物 LTR)最多被部分解析;某个被拆出的子组若仍双峰,当前实现**不会再继续拆**。

### 8.3 Fragment assembly(`refine_assemble_fragments`)

用的是**空间共现 sweep-line**,不是序列相似度——找那些实例在基因组上反复彼此靠近(邻近距离 `D = median_cons·4`,clamp 到 500–30000 bp)的家族。

- **护栏:** ≥3 次共现、同向 ≥0.80、size-ratio ≥0.10、**nesting guard**(≥50% 包含则跳过)、以及 gap 几何 sanity(median gap、MAD)。
- 构建拼接共识,只有当其 MDL 分数超过**各部分之和**时才接受。

> **嵌合风险(已知限制):** nesting guard 只拦**包含**(一个元件的实例整段落在另一个里)。它**拦不住**两个不相关 TE 家族在 TE 富集区(着丝粒周边 / 嵌套插入热点)反复**相邻**共现的情况——在 `D` 最大可达 30 kb 时,一个 LTR 片段与邻近 LINE 片段若 size-ratio 与同向性都过关,可能被错误拼成一条嵌合共识。MDL 门与「≥3 次共现」只是经验性缓解、不能根除。检查输出时留意异常长(>10 kb)的共识,它可能是相邻但无关 TE 的嵌合拼接。

### 8.4 Prune —— 排他覆盖(`refine_prune_families`,sweep-line)

- 把被接受的家族按「最弱优先」排序;对每个实例算其**排他碱基**(没有被任何其它家族覆盖的位置)。
- 一个实例只有在其长度 ≥ **25%** 是排他时才「算数」;一个家族在**没有任何实例**越过 25% 时被 prune(`CAND_ACCEPT_PRUNED`)。
- 操作判据是这个 25%-排他规则,**不是**直接的「分数 vs 代价」比较。

### 8.5 Tandem coalesce(`refine_coalesce_tandem_instances`,选择后,仅报告)

把同链、间隔满足 `−10 ≤ gap ≤ coalesce_factor·consensus_length`(地板 50)的连续实例合并;默认 `coalesce-factor = 20.0`。会重算分歧度,**不改变**库里有哪些家族。

### 8.6 ⚠️ 跨记录安全(正确性关键)

**每一处跨实例的位置比较都必须先 guard `Instance.seq_index`**,使不同 FASTA 记录的位置永不被当作相邻。merge 的包含 / 重叠 helper、assembly sweep、coalesce 都带了这些 guard。

> **任何新写的、跨家族比较 `.position` 的代码,必须先检查 `seq_index`。** 这曾是一类系统性的多记录 bug。

---

## 9. MDL 选择

`mdl_select_library`(`mdl.c:281`)。**单遍,没有 R 收敛循环**(`R_estimate = num_families`,只用一次)。家族按分数降序、贪心地在一个**两分支门**下接纳;只有**排他(不重叠)的** savings 才累加进报告总量,从而保持两部分编码的界。

一个家族被接纳,当**满足其一**:

1. **排他分支** —— 它有未被覆盖的(排他)实例碱基,且 `exclusive_savings − model_cost > 0`;**或**
2. **standalone 回退分支** —— `standalone_score > 0` **且** 共识长度 ≥ 50 **且** 实例 ≥ 3。这一支用来挽救那些「因为拷贝与已接纳家族重叠、会被纯 unique-coverage 贪心误杀」的真实家族;这样被接纳的家族标记 `CAND_ACCEPT_STANDALONE` / `CAND_QF_STANDALONE_FALLBACK`。

**预筛拒绝:** `mdl_score ≤ 0` 或 `< 2` 实例。

**覆盖跟踪是区间 sweep-line**(`MdlInterval`,二分 + 双指针归并),**不是**每碱基 bitmap——内存是 O(#instances),这是多 Gb 基因组可行的关键。

**压缩比**(`mdl.c:713`):`dl_total / (2N)`,其中 `dl_total = 2N − total_savings + dl_library`,clamp 到 ≥ 0。这个 clamp 是这里唯一强制的硬不变量。

---

## 10. 大基因组扩展

两个机制按 `raw_length`(未 pad 的真实碱基数)自动触发(`main.c`):

- **基因组采样**(> `sample-size`,默认 1 Gb)。把基因组缩成有代表性的样本:对 `window-size` 大小的 tile 做 Fisher–Yates(`-seed` 可复现),选出若干 `window-size` 大小的 tile(默认 1 Mb,可调)。发现在样本上跑,实例坐标在输出前**重映射回原基因组**。可用 `-sample-output` 把采样基因组写出来复现。
- **分块发现**(> `chunk-size`,默认 200 Mb)。把长序列切段(段间有边界重叠),用 LPT 把段均衡分到各线程 bin,并行发现,各 chunk 坐标各自重映射。每个 chunk 自算 l-mer 长度以提升灵敏度。在 `chunk-size·1.8` 以上会递归二分。

所有位置 / 长度都是 64-bit(`gpos_t` / `glen_t = int64_t`),避免溢出。

两者自然复合:采样(若触发)先把 10+ Gb 压到 ~1 Gb → 分块(若触发)再把 ~1 Gb 切成 ~200 Mb 并行 chunk → 精修在整个采样基因组上跑(不是逐 chunk)→ 若采样过,输出坐标重映射回原基因组。小基因组(< chunk-size)不切分,发现单线程跑。

---

## 11. Recall rescue(召回补救)

`rescue.c`,opt-in 的 `-recall-rescue`(流水线步骤 2b)。这是一个用运行时间换灵敏度的可选第二遍发现:用**更短的 l-mer**(`l − rescue-l-delta`,地板 8)重跑**进程内**发现引擎,抓主遍漏掉的家族,再把非重复的家族追加回候选列表、标记 `CAND_QF_RESCUE_DISCOVERY`。

- **Targeted 模式(默认):** 从**未覆盖的基因组 gap**(现有实例区间合并后,长度 ≥ `rescue-min-gap` 且两侧 ±L 的 gap)构建补救段,只在那里重发现,再映射回去。
- **Full-genome 模式(`-rescue-full-genome`):** 在整条基因组上重发现(大则分块)。
- **重复门:** 当一个补救家族相对某已有家族 length-ratio ≥ 0.80、identity ≥ 0.80、containment ≥ 0.80(正向 / RC,滑窗 identity)时,作为 duplicate 丢弃。`-rescue-audit <file>` 记录逐 target / 逐 candidate 的决策。整个过程不涉及任何外部工具。

---

## 12. 外部工具现状

**默认情况下,二进制不 fork 任何子进程。** 这一点要特别说明,因为更早的散文档暗示了 BLAST 招募:

- **rmblastn / BLAST 短元件招募代码存在,但未接线(UNWIRED)。** `align_blast_recruit_short_families` 与 `find_rmblastn`(优先环境变量 `$RMBLASTN_BIN`,否则 `which rmblastn` 从 PATH 解析;**无硬编码路径**)在 `align.c` 里,但**没有任何活的调用者**——只有 `main.c` 里一条 revert 注释。RepeatMasker 任何地方都没被调用。
- **流水线唯一能跑的外部工具是 `seqkit stats`**,且仅当给了 `-external-qc <file>`(它会把 `-external-tools` 自动提升为 `auto`)。它对**已写好的** FASTA 做非破坏性 QC、写一个 TSV;绝不回馈到发现 / 精修。`seqkit` 经 PATH 或 `-seqkit <path>` 定位。只有在 `-external-tools require` 下缺工具才硬失败;否则缺工具只是软警告。
- **`tool_runner.c`** 是一个通用、零依赖的 `fork`/`execv` 启动器(带超时,默认 300 s);目前只有 `external_qc.c` 调它。

---

## 13. 输入与输出格式

### 13.1 输入

- **`-sequence <file>`**:基因组 FASTA(可多记录)。
- 可选 **`-freq <file>`**:预算的 l-mer 频率表(`build_lmer_table` 格式),用 `-freq-output <file>` 写出复用。

### 13.2 FASTA 库(`-output`)

每个被接受的家族一条共识。Header schema(v6.1+,**向后兼容**——只匹配 `>R=N length=L copies=C mdl=M` 的旧解析器仍可工作,新字段是追加的):

```
>R=42 length=312 copies=15 mdl=1204.3 div=0.082 topo=linear accept=exclusive tier=... flags=0x00000000 qflags=...
ACGTACGT...
```

字段含义:

| 字段 | 含义 |
|---|---|
| `R=` | 家族编号 |
| `length=` | 共识长度 |
| `copies=` | 实例数 |
| `mdl=` | MDL 分数(节省的 bit) |
| `div=` | 平均分歧度 |
| `topo=` | 拓扑:`linear` / `cyclic` / `complex` |
| `accept=` | 接纳状态:`exclusive` / `standalone` / …(见 [§16](#16-核心数据结构)) |
| `tier=` | 质量层级 |
| `flags=` / `qflags=` | quality-flag bitset(十六进制 + 可读名) |

> **共识链向:** 输出共识按家族**被发现时的方向**给出,**不**强制规范到生物学正义链,也不做 TE 分类(`topo=` 只是 linear/cyclic/complex 拓扑,不是 TE class)。下游分类 / 注释不应假定共识就是 sense 链。

### 13.3 BED6 实例(`-instances`)

每个实例一行;若用了采样,坐标会**重映射回原基因组**:

```
chr1    10000   10312   R=42    850   +
```

- 第 4 列:家族编号 `R=<id>`。
- 第 5 列:`1000·(1−divergence)`(0–1000 刻度)。
- 第 6 列:strand(`+` / `−`)。

### 13.4 TSV 统计(`-stats`)

逐家族一行,表头(tab 分隔):

```
family_id  consensus_length  num_instances  divergence_mean  mdl_score  model_cost  topology
  standalone_score  exclusive_score  exclusive_bases  exclusive_instances
  acceptance  quality_tier  quality_flags  quality_notes
```

这是诊断「家族为什么被接纳 / 拒绝」最有用的文件——`acceptance`、`standalone_score`、`exclusive_score`、`exclusive_bases` 直接对应 [§9](#9-mdl-选择) 的两分支门。

### 13.5 诊断输出

- `-trace-dir <dir>`:每个精修阶段后 dump FASTA+BED+TSV。
- `-split-audit <file>`:split 阶段逐决策审计 TSV。
- `-rescue-audit <file>`:recall-rescue 逐 target / candidate 审计 TSV。
- `-external-qc <file>`:对最终 FASTA 跑 `seqkit stats`,写 QC TSV。

---

## 14. 参数完整参考

> 权威、与代码一致的列表始终是 `bin/mdl-repeat`(无参数运行)。下表按代码默认值整理。

### 14.1 必填

| 参数 | 说明 |
|---|---|
| `-sequence <file>` | 输入 FASTA 基因组 |
| `-output <file>` | 输出重复库(FASTA) |

### 14.2 发现(discovery)

| 参数 | 默认 | 说明 |
|---|---|---|
| `-freq <file>` | — | 预算 l-mer 频率表 |
| `-freq-output <file>` | — | 写出 l-mer 频率表供复用 |
| `-l #` | auto `ceil(1+log4(N))` | l-mer(种子)长度 |
| `-L #` | 10000 | 每侧最大延伸距离(bp) |
| `-minthresh #` | 2 | 种子所需最小 l-mer 频率 |
| `-goodlength #` | 30 | 共识长度预过滤下限 |
| `-maxgap #` | 5 | DP 带最大 offset(`MAXOFFSET`) |
| `-match #` | 1 | DP match 分 |
| `-mismatch #` | −1 | DP mismatch 分 |
| `-gap #` | −5 | DP gap 罚分 |
| `-cappenalty #` | −20 | 退出对齐的罚分上限 |
| `-minimprovement #` | 3 | 每步最小总分提升 |
| `-stopafter #` | 100 | 无进展 N 列后停(自适应 `·/10`) |
| `-maxentropy #` | −0.70 | Shannon 熵过滤阈值 |
| `-tandemdist #` | 500 | 计数的同链 l-mer 最小间距 |
| `-maxoccurrences #` | 10000 | **发现期**每种子最大出现数(注:招募期每家族实例上限见 `-max-instances`,§14.4) |
| `-maxrepeats #` | 100000 | 最大发现家族数 |

### 14.3 Recall rescue

| 参数 | 默认 | 说明 |
|---|---|---|
| `-recall-rescue` | off | 跑有界的二次发现(更短 l-mer 种子) |
| `-rescue-full-genome` | off | 在整条基因组而非 gap 上做补救 |
| `-rescue-l-delta #` | 1 | 补救 l-mer 缩短量(最小 l=8) |
| `-rescue-maxrepeats #` | 2000 | 最多追加的补救家族数 |
| `-rescue-min-gap #` | 200 | targeted 补救的最小未覆盖 gap 长度 |

### 14.4 精修(refinement)

| 参数 | 默认 | 说明 |
|---|---|---|
| `-threads #` | 1 | 精修与分块发现的线程数 |
| `-mdl-mode <mode>` | exact | MDL 位置编码:`none` / `exact` / `upper` |
| `-max-divergence #` | 0.30 | 实例接受的最大替换率(0.0–1.0) |
| `-refine-gap #` | −5 | 精修 gap 罚分(高 indel 物种建议 −3) |
| `-refine-maxoffset #` | 12(最大 32) | 精修 DP 带半宽 |
| `-max-dp-cells #` | 10000000 | 共识 merge 的最大 DP cell 数(~40 MB) |
| `-coalesce-factor #` | 20.0(0=关) | tandem 实例 coalesce 的 gap 容忍(碱基) |
| `-max-instances #` | 10000 | 招募期每家族实例上限;超高拷贝家族调大它(否则 `copies=`/savings 被截断,触顶会 warn) |
| `-max-seed-hits #` | 50000 | 招募期原始 seed-hit 上限(自动 ≥ `-max-instances`) |

> 关于 merge/split 内部固定常数(identity/coverage/aligned = 0.80/0.80/80;split 双峰 0.20 等):它们是编译期常量,不通过 CLI 暴露,详见 [§8](#8-精修阶段)。

### 14.5 大基因组

| 参数 | 默认 | 说明 |
|---|---|---|
| `-chunk-size #` | 200 | 分块发现的 chunk 大小(Mb,最小 10) |
| `-sample-size #` | 1000 | 触发采样的阈值(Mb,最小 100) |
| `-window-size #` | 1000 | 采样 tile 大小(kb,100–10000) |
| `-seed #` | 42 | 采样随机种子(可复现) |
| `-sample-output <file>` | — | 写出采样基因组 FASTA |

### 14.6 输出与诊断

| 参数 | 默认 | 说明 |
|---|---|---|
| `-instances <file>` | — | 输出实例 BED |
| `-stats <file>` | — | 输出家族统计 TSV |
| `-trace-dir <dir>` | — | 每个精修阶段后 dump(诊断) |
| `-split-audit <file>` | — | split 决策审计 TSV |
| `-rescue-audit <file>` | — | recall-rescue 审计 TSV |
| `-external-tools <mode>` | off | 外部工具策略:`off` / `auto` / `require` |
| `-external-qc <file>` | — | 对最终 FASTA 写 seqkit stats TSV |
| `-seqkit <path>` | — | seqkit 可执行文件路径 |
| `-v` / `-vv` | off | 详细级别(1 或 2) |

---

## 15. 调参指南与常见场景

下面是基于各参数实际作用的实务建议(不是 benchmark 数字,benchmark 见 `FINAL_REPORT.md`)。

### 15.1 我想提高灵敏度(漏家族太多)

按优先级:

1. **先开 `-recall-rescue`。** 这是为这个目的专门设计的低风险手段——更短的 l-mer 二次发现,带 duplicate 门防止冗余。高分歧 / 低拷贝家族遗漏时首选。
2. **调小 `-l`** 或 **调大 `-rescue-l-delta`**(在 rescue 内)。更短种子 = 更灵敏,但代价是更多噪声种子与更慢。
3. **放宽 `-max-divergence`**(默认 0.30)。如果你的目标家族分歧度更高,实例会在 0.30 处被截掉;调到比如 0.40 能招回更分歧的拷贝。
4. **高 indel 物种用 `-refine-gap -3`**(默认 −5)。usage 里直接给了这条建议。

> 别指望靠下游 polish 工具补灵敏度。项目记录显示:把库链过 Refiner_mdl 之类的外部 polish,经验上贡献 ≈ 0 甚至变差;召回缺口要在 mdl-repeat **内部**解决。

### 15.2 我想要更干净 / 更小的库

- MDL 本身已经是 precision 的主把关者,通常不需要额外动作。
- 极致追求干净时,可在**下游**只应用最保守的 hard filter;不要启用激进的多阶段 polish。
- `-stats` 里的 `acceptance` 列能告诉你哪些家族是靠 `standalone` 回退进来的(相对更边缘)。

### 15.3 大基因组 / 内存或时间紧张

- **给足 `-threads`。** l-mer 计数、精修、分块发现都并行。
- **采样阈值** `-sample-size`(默认 1 Gb):若想在更小基因组上也强制采样以提速,调小它(最小 100 Mb);想关掉采样、全量发现,调到比基因组还大。
- **`-chunk-size`**(默认 200 Mb):更小的 chunk → 更短的 per-chunk l-mer → 对低拷贝家族更灵敏,但调度开销上升。
- **务必固定 `-seed`** 以保证采样可复现。
- 想留存采样输入做复现 → `-sample-output`(与 `-seed` 配合即可重放整个采样)。
- **确定性的范围(实测):** **固定 `-threads`(配合固定 `-seed`)下完全可复现**——同参数重复运行产出逐字节一致的家族集。但**不同 `-threads` 之间家族集合可能有细微差异**(实测一个 200 kb 多记录基因组:`-threads 1` 得 6 个家族、`-threads 4` 得 5 个,相差一个边缘家族)。根因是**并行精修 / 招募路径**(merge/split worker 与实例招募的处理顺序随线程数变化,偶尔翻转某个边缘家族的留/弃),并非只有大基因组的跨-chunk 并行发现。这对库质量影响很小,但若需严格复现,请固定 `-threads`。

### 15.4 MDL 编码模式怎么选

- **`exact`(默认)**:精确二项位置编码 `log2 C(a,m)`,推荐保持默认。
- **`upper`**:位置项用 `m·log2(a)`,是 `exact` 的**上界**——实例代价更高 → savings 更小 → 家族分数更低,因此**更保守、更不灵敏**(注意:不是更灵敏)。
- **`none`**:完全去掉位置项 → 实例代价最低 → 最宽松 / 最灵敏,但失去位置编码的理论依据,一般仅作向后兼容。

### 15.5 排查「家族在哪一步消失」

- 加 `-vv` 看每阶段家族数。
- 加 `-trace-dir <dir>`,逐阶段对比 dump 的家族集。
- 看 `-stats` 的 `acceptance` / `exclusive_*` / `standalone_score`,判断是被 MDL 选择拒、还是被 prune 的 25%-排他规则拒。
- split 异常 → `-split-audit`;rescue 行为 → `-rescue-audit`。

### 15.6 空输出与边界场景(不是 bug 的常见原因)

- **极小基因组 → 合法空输出。** 自动 l-mer 长度 `l = ceil(1+log4(N))` 对亚 kb 输入只有 4–6;叠加 `GOODLENGTH=30` 与「≥3 拷贝」的 MDL 门,几 kb 的玩具基因组**完全可能 0 家族、空 FASTA**(流水线会提前空输出退出)。这是预期,不是工具坏了——要有意义的结果,基因组需有足够长度与真实高拷贝重复。手动测试可考虑用合成高拷贝序列(见 `tests/`)。
- **N 密集 / gap 多的装配。** N 编码为 99 并被熵 / mask 逻辑抑制,大量 N 区会静默压低可用种子;这是设计行为,不会报错。
- **低复杂度 / satellite 基因组。** 熵、周期性、tandem-distance 三道过滤会主动排除这些(它们属 [§1.1](#11-适用范围scope) 的 out-of-scope),所以以串联 / 卫星为主的序列产出很少家族属正常。
- **单记录 vs 多记录。** [§8.6](#86-️跨记录安全正确性关键) 的 `seq_index` guard 只在多 FASTA 记录间起作用;单记录输入天然不触及这类跨记录融合问题。BED 实例按记录(染色体)给坐标。

---

## 16. 核心数据结构

文件:`candidates.h`、`types.h`、`genome.h`。

- **`gpos_t` / `glen_t = int64_t`** —— 每个基因组位置与长度。> 2 Gb 基因组必须 64-bit。
- **DNA 编码:** A/C/G/T = 0/1/2/3,N = 99;互补 = `3 − c`(N-safe)。每条加载的基因组前端 pad `PADLENGTH = 11000`。
- **`Genome`:** 数值化 `sequence`(已 pad)、`length`(含 pad)、`raw_length`(真实碱基数——驱动采样 / 分块阈值)、`boundaries`(逐记录)、`sequence_ids`。
- **`Instance`:** `position`(pad 坐标)、`aligned_length`、`cons_start/end`、`num_edits`、`divergence`(0–1)、`score`、`strand`(±1)、**`seq_index`**(FASTA 记录号——任何跨实例比较的正确性关键)。
- **`CandidateFamily`:** `consensus`(数值)、`consensus_length`、`topology`(linear/complex/cyclic)、`estimated_copies`、`instances[]`、`discovery_flags`,以及一个 `CandidateMdlState mdl`(model cost、standalone/exclusive savings & scores、exclusive bases/instances、accept state、quality tier)。迁移期保留 `mdl_score` / `model_cost` 别名并同步。
- **`CandidateAcceptState`:** `UNSCORED=0, REJECTED, EXCLUSIVE, STANDALONE, PRUNED`(0 = 尚未选择,所以 `memset` 初始化合法)。
- **`CAND_QF_*`** quality-flag bitset —— 尤其 `CAND_QF_RESCUE_DISCOVERY` 与 `CAND_QF_STANDALONE_FALLBACK`。

---

## 17. 不变量与坑(改代码前必读)

- **每一处跨实例位置比较都要 guard `seq_index`。** 漏一处就会静默把不同 FASTA 记录的无关 locus 融合。
- **位置一律 64-bit** —— 永远别把基因组偏移存进 `int`。
- **MDL 代价在当前公式里与 R 无关。** 不要加入「假设分数随家族数变化」的控制流,除非先恢复一个 R 相关代价项(见 [§4.5](#45-️载重级注意代价与-r-无关))。
- **`PADLENGTH`(11000)必须 ≥ 最大单侧延伸距离**,否则延伸 / masking / rescue 的边界数学会坏。
- **MDL 不变量:** 压缩比 ∈ [0,1];`dl_total` clamp ≥ 0;只有排他 savings 累加进报告总量。改选择时务必保持「排他-OR-standalone」的接纳语义。
- **精修变换是 MDL-gated 的(逐变换局部门)。** merge/split/assembly 各有一道局部 MDL 接受检查(split 为放宽门槛);它保证单个变换在自身评分下不劣化,但这是局部性质、**不是**对最终库总 DL 的单调保证——改代码时保持这道局部门。
- **主指标是家族级 recall**(cd-hit + BLAST,RepeatModeler2/EDTA 风格)。per-instance bp recall 是噪声代理,会低估库质量,别拿它当主数字。

---

## 18. 源码地图

```
src/
  main.c              流水线驱动、CLI 解析、采样、分块发现、usage()
  types.h             gpos_t/glen_t (int64_t)、DNA 编码、PADLENGTH、常量
  genome.c/.h         FASTA 加载、前端 pad、逐记录边界、采样
  kmer.c/.h           Canonical k-mer 计数(striped-lock 并行)、位置索引
  discover.c/.h       seed-and-extend 引擎;l-mer 表;N 序列带状延伸
  discover_internal.h discover.c 与 discover_mask.c 间的共享类型
  discover_mask.c     家族发现后的 masking 子系统(1-vs-1 带状 DP)
  align.c/.h          多 k-mer 种子、实例招募/重对齐;(未接线的 rmblastn)
  candidates.c/.h     CandidateFamily / Instance / CandidateList + MDL 状态
  refine.c/.h         merge / split / fragment-assembly / prune / tandem-coalesce
  mdl.c/.h            L_int (Rissanen)、逐家族打分、单遍选择
  rescue.c/.h         可选的 recall-rescue 二次发现 pass
  external_qc.c/.h    opt-in 的 seqkit-stats QC 策略层
  tool_runner.c/.h    通用 fork/execv 子进程启动器(带超时)
  output.c/.h         FASTA / BED6 / TSV 写出器
  cmd_line_opts.c/.h  通用 CLI 取值解析
```

---

## 19. 参考文献

1. Price, A.L., Jones, N.C., & Pevzner, P.A. (2005). De novo identification of repeat families in large genomes. *Bioinformatics*, 21(suppl_1), i351–i358.(RepeatScout,seed-and-extend 与 N 序列同时延伸的来源)
2. Rissanen, J. (1978). Modeling by shortest data description. *Automatica*, 14(5), 465–471.(MDL 与通用整数码)
3. Otsu, N. (1979). A threshold selection method from gray-level histograms. *IEEE Trans. Systems, Man, and Cybernetics*, 9(1), 62–66.(split 阶段双峰阈值)
4. Grumbach, S., & Tahi, F. (1994). A new challenge for compression algorithms: genetic sequences. *Information Processing & Management*, 30(6), 875–886.
