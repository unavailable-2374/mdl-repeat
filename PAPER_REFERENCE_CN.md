# mdl-repeat：论文写作完整技术参考文档

**用途**：本文档提供 mdl-repeat 所有算法细节、数学公式、数据结构、阈值和复杂度分析，作为 Methods 论文写作的参考素材。每个公式和常量均可追溯到源码的精确位置。

---

## 目录

1. [问题定义与 MDL 框架](#1-问题定义与-mdl-框架)
2. [发现引擎：N 序列同步带状 DP](#2-发现引擎)
3. [K-mer 哈希表与位置索引](#3-k-mer-哈希表与位置索引)
4. [种子锚定带状比对（精炼阶段）](#4-种子锚定带状比对)
5. [精炼流水线](#5-精炼流水线)
6. [MDL 打分与文库筛选](#6-mdl-打分与文库筛选)
7. [大基因组支持](#7-大基因组支持)
8. [数据结构](#8-数据结构)
9. [完整参数参考](#9-完整参数参考)
10. [复杂度分析](#10-复杂度分析)

---

## 1. 问题定义与 MDL 框架

### 1.1 两部分编码（Two-Part Code）

给定长度为 $N$ 的基因组 $G$，重复序列文库 $L$ 是 $R$ 条共识序列（consensus）的集合。最小描述长度（MDL）原理选择使总描述长度最小化的文库：

$$DL(G) = DL(L) + DL(G \mid L)$$

其中 $DL(L)$ 是描述文库本身的代价，$DL(G \mid L)$ 是在已知文库条件下描述基因组的代价（Rissanen, 1978）。不需要任何外部阈值（最小拷贝数、最小共识长度等）。

### 1.2 Rissanen 通用整数编码

编码正整数 $n$：

$$L_{\text{int}}(n) = \log_2(c_0) + \log_2(n) + \log_2(\log_2(n)) + \log_2(\log_2(\log_2(n))) + \cdots$$

仅累加非负项。归一化常数 $c_0 \approx 2.865064$ 保证满足 Kraft 不等式（$\log_2(c_0) = 1.5179605508986484$）。

**实现**（`mdl.c:15-34`）：迭代循环累加 `log2(val)` 直至 `val <= 0.0`，然后加上 `LOG2_C0`。

参考值：

| $n$ | $L_{\text{int}}(n)$（比特） |
|-----|---------------------------|
| 1   | 1.518                     |
| 2   | 2.518                     |
| 3   | 3.768                     |
| 10  | 7.364                     |
| 100 | 12.344                    |
| 1000| 17.830                    |

### 1.3 文库编码代价

$$DL(L) = L_{\text{int}}(R) + \sum_{r=1}^{R} \left[ L_{\text{int}}(\text{len}_r) + 2 \cdot \text{len}_r \right]$$

每个家族的模型代价：

$$\text{model\_cost}_r = L_{\text{int}}(\text{len}_r) + 2 \cdot \text{len}_r$$

$2 \cdot \text{len}_r$ 项以 2 bits/base 编码共识序列（均匀碱基分布，有意保守）。（`mdl.c:123-126`）

### 1.4 单实例编码代价

一个比对长度为 $a_i$、编辑数为 $m_i$ 的重复实例：

$$C_{\text{instance}}(a_i, m_i) = L_{\text{int}}(a_i) + L_{\text{int}}(m_i + 1) + m_i \cdot \log_2 3 + P(a_i, m_i)$$

其中：
- $L_{\text{int}}(a_i)$：比对区域长度
- $L_{\text{int}}(m_i + 1)$：编辑数（+1 避免 $L_{\text{int}}(0)$）
- $m_i \cdot \log_2 3$：编辑类型指定（替换为 3 种可能碱基之一）
- $P(a_i, m_i)$：编辑位置编码

**位置编码模式**（通过 `-mdl-mode` 选择）：

| 模式 | 公式 | 说明 |
|------|------|------|
| `exact`（默认） | $\log_2 \binom{a_i}{m_i} = \frac{\ln\Gamma(a_i+1) - \ln\Gamma(m_i+1) - \ln\Gamma(a_i - m_i + 1)}{\ln 2}$ | 通过 lgamma 精确计算二项式系数 |
| `upper` | $m_i \cdot \log_2(a_i)$ | 上界（每个编辑位置独立编码） |
| `none` | $0$ | 省略（向后兼容） |

（`mdl.c:66-107`）

**截断处理**：$a_i$ 截断至 $[1, \infty)$；$m_i$ 截断至 $[0, a_i]$。（`mdl.c:76-78`）

### 1.5 MDL 筛选准则

每个家族：

$$\text{mdl\_score}_r = \sum_{i=1}^{n_r} \left( 2 \cdot a_i - C_{\text{instance}}(a_i, m_i) \right) - \text{model\_cost}_r$$

一个家族被**接受**当且仅当 $\text{mdl\_score}_r > 0$ 且 $n_r \geq 2$。（`mdl.c:128-156`）

### 1.6 设计决策：省略单实例开销

理论上完整的编码应包含单实例开销：类型位（1 bit）、家族标识符（$\lceil\log_2 R\rceil$ bits）、链位（1 bit）、共识起始指针（$\lceil\log_2(\text{len}_r)\rceil$ bits）——每个实例共约 19–24 额外比特。在真实基因组上测试显示 3–6% 的灵敏度损失。当前公式省略了这些项，优先保证灵敏度而非理论完备性。`mdl_instance_cost_full()` 函数接受 `consensus_length` 和 `num_families` 参数以保留 API 扩展性，但目前不使用它们。（`mdl.c:54-64`）

**结果**：迭代 R 收敛循环（§6.1）实际上处于无效状态（inert），因为单实例代价不依赖于 R。

---

## 2. 发现引擎

### 2.1 概述

发现引擎（`discover.c`，约 2000 行）实现了 RepeatScout（Price et al. 2005）的种子-扩展范式，采用线程安全的 `DiscoverContext` 封装。所有可变状态存储在堆分配的结构体中，支持并行分块发现。

### 2.2 基因组加倍

发现前，基因组被加倍以捕获两条链：

```
[PADLENGTH 个 N] [正向基因组] [PADLENGTH 个 N] [反向互补基因组]
```

- `PADLENGTH = 11000`（`types.h:22`）
- 正向拷贝：位置 $[\text{PADLENGTH}, N + \text{PADLENGTH})$
- RC 拷贝：位置 $[N + 2 \cdot \text{PADLENGTH}, 2N + 2 \cdot \text{PADLENGTH})$
- 总长度：$2N + 2 \cdot \text{PADLENGTH}$

基因组加倍约使种子出现次数翻倍，为同步 DP 扩展提供更强的统计证据。（`discover.c:1575-1597`）

### 2.3 L-mer 频率计数

**L-mer 长度**：自动计算为 $l = \lceil 1 + \log_4 N \rceil$，最大 31。（`discover.c:1758-1759`，`main.c:36-39`）

**哈希函数**（`discover.c:126-141`）：

$$h_{\text{fwd}} = \left(\sum_{x=0}^{l-1} 4 \cdot h + (\text{base}_x \bmod 4)\right) \bmod H$$
$$h_{\text{rc}} = \left(\sum_{x=0}^{l-1} 4 \cdot h + (3 - \text{base}_{l-1-x}) \bmod 4\right) \bmod H$$
$$h = \max(h_{\text{fwd}}, h_{\text{rc}})$$

其中 $H = 16{,}000{,}057$（素数）。对称哈希确保正向和 RC 的 l-mer 映射到同一桶。

**TANDEMDIST 过滤**（`discover.c:332-378`）：某条链上的 l-mer 出现只有在与该链上前一个同向出现的距离超过 `TANDEMDIST`（默认 500 bp）时才被计数。通过 `lastplusocc`（正向）和 `lastminusocc`（反向互补）分链追踪，防止串联重复导致计数膨胀。

**两遍位置去重**（`discover.c:459-506`）：计数后，第二遍使用带符号的位置标记移除同链上 `TANDEMDIST` 以内的串联重复。

### 2.4 N 序列同步带状 DP 扩展

给定一个具有 $N$ 个基因组出现位置的种子 l-mer，**所有 $N$ 个序列同时扩展**，使用带状动态规划。

**DP 状态**：`score[2][N][2 \cdot \text{MAXOFFSET} + 1]`
- 维度 1：当前/前一行（滚动缓冲）
- 维度 2：每个出现位置一个条目（$N \leq$ `MAXN` = 10000）
- 维度 3：从 $-\text{MAXOFFSET}$ 到 $+\text{MAXOFFSET}$ 的带偏移量

**右扩展递推关系**（`discover.c:594-655`）：

对共识位置 $y$、序列 $n$、偏移 $o$ 和候选碱基 $a$：

$$S_y(n, o) = \max \begin{cases}
S_{y-1}(n, o+1) + \text{GAP} & \text{（序列中的空位，} o < \text{MAXOFFSET}） \\[4pt]
S_{y-1}(n, o) + \delta(a, g_n(o,y)) & \text{（对角线匹配/错配）} \\[4pt]
\max_{o' \in [-M, o)} \left[ S_{y-1}(n, o') + (o - o') \cdot \text{GAP} + \delta'(a, g_n, o', o, y) \right] & \text{（共识中的多空位）}
\end{cases}$$

其中：
- $\delta(a, b) = \text{MATCH}$（若 $a = b$），否则 $\text{MISMATCH}$（默认：+1 / -1）
- $g_n(o, y)$：由 `pos[n]`、偏移 $o$ 和共识位置 $y$ 计算的基因组碱基
- $\delta'$：若 $a$ 匹配空位范围 $[o', o]$ 中的**任意**碱基则返回 MATCH，否则 MISMATCH
- $M = \text{MAXOFFSET}$（发现阶段默认为 5，`discover.c:1607`）

**左扩展**（`discover.c:662-723`）：与右扩展对称但方向相反。使用 `(w+1) % 2` 作为行索引。不对称性：左扩展**不**保存 `score_of_besty` 和 `savebestscore` 检查点。

**链处理**：对于反向互补出现（`rev[n] = 1`），基因组位置映射为 `pos[n] - (offset + y - L - l) - 1`，并查找互补碱基。

**边界检查**（`discover.c:601-608`）：每次 DP 单元计算前，验证基因组位置是否在通过 `get_boundaries()` 计算的序列边界 $[\text{bStart}, \text{bEnd})$ 内。

### 2.5 共识碱基选择

在每个扩展位置 $y$，对每个候选碱基 $a \in \{A, C, G, T\}$：

$$\text{total}(a) = \sum_{n=0}^{N-1} \max\left(0, \max_{o} S_y(n, o) - \text{CAPPENALTY}\right)$$

其中单序列最优分数的下限为 `bestbestscore[n] + CAPPENALTY`（默认 CAPPENALTY = -20）。

**共识碱基选择为** $\arg\max_a \text{total}(a)$ —— 纯 argmax 取累加最优分数，不是加权多数投票。（`discover.c:862-886`）

### 2.6 终止条件

三个交互参数控制扩展终止：

**MINIMPROVEMENT**（`discover.c:760-773`）：仅当满足以下条件时更新检查点：
$$\text{totalbestscore}_y \geq \text{besttotalbestscore} + (y - y^*) \cdot \text{MINIMPROVEMENT}$$
其中 $y^*$ 是上一个检查点位置。默认 MINIMPROVEMENT = 3。

**WHEN_TO_STOP**（`discover.c:889`）：若 $y - y^* \geq \text{WHEN\_TO\_STOP}$ 则停止扩展。默认 = 100。

**CAPPENALTY**（`discover.c:867`）：单序列分数下限：$\max(0, \text{bestbestscore}[n] + \text{CAPPENALTY})$。默认 = -20。

### 2.7 Shannon 熵过滤

扩展后，低复杂度序列被拒绝：

$$H = \sum_{b \in \{A,C,G,T\}} \frac{c_b}{l} \cdot \ln\left(\frac{c_b}{l}\right)$$

其中 $c_b$ 是 l-mer 中碱基 $b$ 的计数。**若 $H > \text{MAXENTROPY}$（默认 -0.70，自然对数，负值）则拒绝**。（`discover.c:188-202`）

### 2.8 掩蔽（Masking）

发现一个家族后，通过 1-vs-1 带状 DP 比对（而非 N 序列同步 DP）掩蔽其出现位置。掩蔽 DP（`compute_maskscore_right/left`，`discover.c:985-1066`）使用与 §2.4 相同的递推关系，但只对单条序列与共识进行比对。扩展持续至 `WHEN_TO_STOP` 个连续无改进位置。被掩蔽的位置设为 DNA_N。（`discover.c:1073-1216`）

### 2.9 从扩展中收集实例

对参与扩展的每个出现 $n$（`discover.c:1501-1565`）：

**正向链**（$\text{rev}[n] = 0$）：
- $\text{genome\_start} = \text{pos}[n] + w^* - L$
- $\text{genome\_end} = \text{pos}[n] + y^* - L + 1$

**反向链**（$\text{rev}[n] = 1$）：
- $\text{genome\_start} = \text{pos}[n] - (y^* - L - l) - 1$
- $\text{genome\_end} = \text{pos}[n] - (w^* - L - l)$

其中 $y^*$ 和 $w^*$ 分别是最优右扩展和左扩展位置。编辑数通过与共识序列逐碱基直接比较计算。仅报告映射到加倍基因组正向拷贝的实例（位置 < $N + \text{PADLENGTH}$）。

### 2.10 DiscoverContext：线程安全封装

所有可变状态存储在堆分配的 `DiscoverContext` 结构体中（`discover.c:62-110`），包含：

| 区域 | 字段 | 用途 |
|------|------|------|
| 参数 | `l, L, MAXOFFSET, MAXN, MAXR, MATCH, MISMATCH, GAP, CAPPENALTY, MINIMPROVEMENT, WHEN_TO_STOP, MAXENTROPY, GOODLENGTH, MINTHRESH, TANDEMDIST` | 算法参数 |
| 基因组数据 | `sequence, sequence_owned, removed, length, orig_length, disc_boundaries, disc_num_sequences` | 加倍基因组和掩蔽数组 |
| 工作空间 | `master, masters[], masterstart[], masterend[], pos[], rev[], upperBoundI[], N` | 扩展缓冲区 |
| DP 数组 | `score[2][MAXN][2M+1], score_of_besty[MAXN][2M+1], maskscore[2][2M+1], bestbestscore[], savebestscore[]` | DP 状态 |
| 循环状态 | `besty, bestw, R, totalbestscore, besttotalbestscore` | 扩展进度 |

公共 API：`CandidateList *discover_families(const Genome *genome, const DiscoverParams *params)` —— 只读输入、无全局状态、多次调用可安全并行运行。（`discover.c:1724-1971`）

---

## 3. K-mer 哈希表与位置索引

### 3.1 规范 K-mer（Canonical K-mers）

K-mer 压缩为 64 位整数（每碱基 2 bit，最大 $k = 31$）：

$$\text{packed}(s) = \sum_{i=0}^{k-1} s_i \cdot 4^{k-1-i}$$

其中 $s_i \in \{0,1,2,3\}$ 为数值编码。（`kmer.c:91-100`）

规范形式：$\text{canonical}(s) = \min(\text{packed}(s), \text{packed}(\overline{s}^R))$，其中 $\overline{s}^R$ 为反向互补。（`kmer.c:112-116`）

### 3.2 哈希表

**哈希函数**：Fibonacci 哈希，分布优良：

$$h(k) = (k \times 11400714819323198485) \bmod T$$

其中 $11400714819323198485 = \lfloor 2^{64} / \varphi \rfloor$（$\varphi$ = 黄金比例）。（`kmer.c:120-125`）

**表大小**：$T = \text{next\_prime}(\max(16{,}000{,}057, N/4))$。（`kmer.c:232-234`）

**冲突解决**：链表法（chaining）。

### 3.3 条带锁并行计数（Striped Lock Parallel Counting）

- $\text{NUM\_STRIPES} = 4096$ 个互斥锁（`kmer.c:131`）
- 锁分配：$\text{stripe} = h(k) \bmod \text{NUM\_STRIPES}$（`kmer.c:162-163`）
- 线程本地内存池：`POOL_BLOCK_SIZE = 4096` 个 KmerEntry/块（`kmer.c:10`）
- 计数时的 TANDEMDIST 过滤：通过 `last_plus_occ` / `last_minus_occ` 分链追踪（`kmer.c:174-185`）

### 3.4 位置索引：两阶段并行构建

**阶段 1 — 计数**（`kmer.c:480-499`）：工作线程扫描分配的基因组区块，原子递增位置计数。

**阶段 2 — 分配**（`kmer.c:560-572`）：对每个 k-mer，分配大小为 $\min(\text{num\_positions}, \text{KMER\_MAX\_POSITIONS})$ 的位置数组，其中 $\text{KMER\_MAX\_POSITIONS} = 50{,}000$（`kmer.h:24`）。

**阶段 3 — 填充**（`kmer.c:574-597`）：按基因组坐标顺序顺序填充以保证确定性。截断时始终保留坐标最小的位置。（`kmer.c:662-665`）

**位置编码**：正值表示正向链，负值表示反向互补链（符号编码链方向）。

### 3.5 修剪（Trimming）

`kmer_trim` 移除频率低于 `MINTHRESH`（默认 2）的 k-mer。（`kmer.c:399-428`）

---

## 4. 种子锚定带状比对（精炼阶段）

精炼流水线使用不同于发现引擎的比对策略：基于精确 k-mer 匹配锚定的多 k-mer 种子带状比对。

### 4.1 共识 K-mer 集合

对长度为 $L_c$ 的家族共识序列，提取 $L_c - k + 1$ 个 k-mer（上限 `MAX_CONS_KMERS` = 10000）。在哈希集合中去重规范 k-mer。（`align.c:27-93`）

### 4.2 种子基因组扫描

对每个唯一的共识 k-mer，从预构建的 KmerTable 位置索引中查找所有基因组位置。复杂度：$O(\sum_j f_j)$，其中 $f_j$ 是共识 k-mer $j$ 的频率 —— 避免 $O(N)$ 基因组扫描。（`align.c:103-149`）

返回至多 `MAX_SEED_HITS = 50000` 个种子命中，每个记录 `(genome_pos, cons_pos, strand)`。

### 4.3 种子命中聚类

命中按 `(strand, genome_pos)` 排序，使用距离阈值聚类：

$$d_{\text{merge}} = \lfloor 1.5 \cdot L_c \rfloor$$

每个聚类内，**锚点**是共识位置最接近 $L_c / 2$（中点）的命中。（`align.c:171-210`）

### 4.4 带状 DP 比对

对每个锚点，带状 DP 比对双向扩展。

**带宽**：$W = 2 \cdot \text{MAXOFFSET} + 1$（默认 MAXOFFSET = 12，故 $W = 25$）。（`align.c:266`）

**DP 递推关系**（`align.c:297-330`）：对共识位置 $i$、带偏移 $o$：

$$D_i(o) = \max \begin{cases}
D_{i-1}(o) + \delta(\text{cons}_i, \text{genome}(i, o)) & \text{（对角线）} \\
D_{i-1}(o+1) + g & \text{（基因组中的空位）} \\
D_{i-1}(o-1) + g & \text{（共识中的空位）} \\
S^* + \text{CAPPENALTY} & \text{（分数下限）}
\end{cases}$$

其中：
- $\delta(a,b) = +1$ 若匹配，$-1$ 若错配（ALIGN_MATCH / ALIGN_MISMATCH）
- $g = -5$（g_align_gap，可通过 `-refine-gap` 配置）
- $S^*$：该方向上目前的最优分数
- CAPPENALTY = -20

**基因组位置映射**（`align.c:235-242`）：
- 正向：$\text{gp} = \text{anchor\_genome} + (i - \text{anchor\_cons}) + o$
- 反向：$\text{gp} = \text{anchor\_genome} - (i - \text{anchor\_cons}) + o$

**停滞检测**：连续 `ALIGN_WHEN_TO_STOP = 100` 个无改进位置后停止扩展。（`align.c:339-340`）

**散度计算**（`align.c:423-445`）：仅计替换，不计空位：

$$\text{divergence} = \frac{\text{edits}}{\text{compared\_positions}}$$

**拒绝条件**：散度 $> 0.30$（`g_align_max_divergence`）则拒绝该实例。（`align.c:445`）

### 4.5 实例去重

新实例与已有较短实例重叠超过 50% 则被拒绝：

$$\frac{\text{overlap}}{\min(\text{len}_{\text{new}}, \text{len}_{\text{existing}})} > 0.5 \implies \text{拒绝}$$

（`align.c:525-545`）

### 4.6 共识重建：分数加权多数投票

对每个共识位置 $p$，碱基由加权投票选出：

$$\text{consensus}[p] = \arg\max_{b \in \{A,C,G,T\}} \sum_{i : p \in [\text{cons\_start}_i, \text{cons\_end}_i)} w_i \cdot \mathbb{1}[\text{genome}_i(p) = b]$$

其中 $w_i = \max(1, \text{score}_i)$ 为比对分数。（`align.c:579-637`）

### 4.7 共识扩展

双向扩展，通过多数投票在侧翼上下文中进行：

**动态支撑阈值**（`align.c:707-708`）：
- 扩展 > 1000 bp：要求 $\geq 3$ 个支撑实例
- 扩展 > 5000 bp：要求 $\geq 5$ 个支撑实例
- 始终要求 $\geq 2$ 个实例且最优碱基 $\geq 50\%$ 一致

**扩展上限**（基于实例数，`align.c:791-796`）：
- 2–3 个实例：最多 500 bp
- 4–9 个实例：最多 3000 bp
- $\geq 10$ 个实例：最多 10000 bp

**来源资格**：仅距共识边缘 `EXTENSION_SLACK = 15` bp 以内的实例可贡献扩展。（`align.c:677-679`）

### 4.8 迭代精炼

循环：`收集实例 → 扩展 → 重建`。最多 `ALIGN_MAX_ITERATIONS = 10` 次迭代。收敛条件：无扩展 且 (<1% 共识变化或零变化)。（`align.c:808-829`）

### 4.9 并行化

通过 `__atomic_fetch_add` 在共享家族计数器上实现工作窃取（work-stealing）。每个线程获取下一个未处理的家族。（`align.c:835-952`）

---

## 5. 精炼流水线

### 5.1 合并：80-80-80 规则

**阶段 1：K-mer profile 筛选**（`refine.c:26-57`）
- 8-mer 位集 profile：`PROFILE_BITS = 65536`（$4^8$），存储为 1024 个 `uint64_t` 字
- 通过 popcount 计算 Jaccard 指数：$J = |A \cap B| / |A \cup B|$
- $J < 0.15$ 的配对跳过（`REFINE_MIN_JACCARD`）

**阶段 2：半全局比对**（`refine.c:79-183`）

DP 打分：匹配 $+2$，错配 $-3$，空位 $-2$（`REFINE_MATCH/MISMATCH/GAP`）。

仅对查询序列（较短序列）允许首尾自由空位（free end-gaps）。

DP 单元限制：$(\text{qlen}+1) \times (\text{tlen}+1) \leq 10^7$（`g_refine_max_dp_cells`）。（`refine.c:87-92`）

测试正向和反向互补两种方向。

**合并标准**（`refine.c:543-569`）：

| 标准 | 一致性 | 覆盖度 | 比对碱基数 | 附加条件 |
|------|--------|--------|-----------|----------|
| **严格** | $\geq 80\%$ | $\geq 80\%$ | $\geq 80$ bp | — |
| **宽松** | $\geq 70\%$ | $\geq 70\%$ | — | $\geq 50\%$ 实例重叠 |
| **DP 限制回退** | — | — | — | Jaccard $\geq 0.80$ + $\geq 50\%$ 实例重叠 |

**阶段 3：传递合并**——通过带路径压缩的并查集（union-find）实现（`refine.c:294-311`）。代表家族：拥有最多实例的家族。

**阶段 4：合并后重精炼**——通过 `align_refine_family()` 对吸收了其他家族的代表家族进行重精炼。（`refine.c:683-705`）

### 5.2 拆分：基于散度分布的 Otsu 方法

**直方图**：100 个 bin，覆盖散度范围 $[d_{\min}, d_{\max}]$。（`refine.c:753-762`）

**Otsu 阈值**（`refine.c:787-802`）：对每个候选阈值 $t$，计算类间方差：

$$\sigma^2_B(t) = \frac{w_0 \cdot w_1 \cdot (\mu_0 - \mu_1)^2}{n^2}$$

其中 $w_0, w_1$ 是类大小，$\mu_0, \mu_1$ 是类均值。选择 $t^* = \arg\max_t \sigma^2_B(t)$。

**双峰性分数**：$B = \sigma^2_B(t^*) / \sigma^2_{\text{total}}$。（`refine.c:804`）

**接受标准**（`refine.c:944-1061`）：

| 检查 | 阈值 | 来源 |
|------|------|------|
| 最少实例数 | $\geq 10$ | `REFINE_MIN_SPLIT_INSTANCES` |
| 双峰性分数 | $B \geq 0.20$ | `REFINE_BIMODALITY_THRESH` |
| 谷深度检查（若 $0.20 \leq B < 0.40$） | 谷高度 $< \text{lower\_mode\_height}/2$ | `refine.c:807-825` |
| 最小聚类大小 | 每个子家族 $\geq 3$ | `REFINE_MIN_CLUSTER_SIZE` |
| 平均散度差 | $\bar{d}_{\text{hi}} - \bar{d}_{\text{lo}} \geq 0.05$ | `REFINE_MIN_DIV_GAP` |
| MDL 验证 | $\text{mdl}_{\text{lo}} + \text{mdl}_{\text{hi}} > \text{mdl}_{\text{orig}}$ | `refine.c:1053-1061` |

**共识重建**：从重新分配的实例中通过加权多数投票重建（`refine.c:835-888`）。

### 5.3 片段组装：空间共现（Spatial Co-occurrence）

转座子可能被发现为多个片段。不重叠的片段不共享任何 k-mer，因此使用基因组中的**空间共现**替代 k-mer Jaccard。

**距离阈值**（`refine.c:1714-1721`）：

$$D = \min(\max(2 \times \text{median\_consensus\_length}, 500), 5000)$$

**扫描线算法**（`refine.c:1749-1811`）：
- 为所有实例构建按位置排序的数组 `(start, end, family_id, instance_idx, strand)`
- 二次抽样：实例数 >10K 的家族缩减至 10K
- 对距离 $D$ 以内来自不同家族的每对实例：递增共现计数器并追踪方向一致性

**配对过滤**（`refine.c:1830-1879`）：

| 检查 | 阈值 |
|------|------|
| 共现次数 | $\geq 3$ |
| 嵌套保护 | 任一家族均不应 $\geq 50\%$ 包含在另一个中 |
| 大小比例 | $\min(\text{len}_A, \text{len}_B) / \max(\text{len}_A, \text{len}_B) \geq 0.10$ |
| 方向一致性 | $\geq 80\%$ 同向 |

**间隔分析**（`refine.c:1936-1996`）：
- 计算共现实例对之间的间隔分布
- 中位间隔和中位绝对偏差（MAD）
- MAD > 100（间隔不一致）或中位间隔 $\notin [-20, \min(\text{median} + 2\text{MAD}, 500)]$ 则拒绝
- 允许至 -20 bp 的负间隔（靶位点重复 TSD 的重叠）

**组装**（`refine.c:1998-2062`）：
- 正间隔：`consensus_A + N 填充 + consensus_B`
- 重叠（$\leq 0$）：`consensus_A[0..len_A-overlap) + consensus_B`

**MDL 验证**：仅当 $\text{mdl}_{\text{assembled}} > \text{mdl}_A + \text{mdl}_B$ 时接受。（`refine.c:2064-2109`）

**链长度上限**：并查集最大链长度 10，防止失控合并。（`refine.c:1903-1925`）

### 5.4 修剪：独占覆盖度（Exclusive Coverage）

移除独占基因组贡献不足以证明其模型代价的边缘家族。

**算法**（`refine.c:1487-1605`）：
1. 按 MDL 分数升序排列已接受家族（最弱的优先）
2. 构建基因组覆盖数组：`cov[pos]` = 覆盖位置 `pos` 的已接受家族数
3. 对每个家族（最弱优先）：
   - 对每个实例，计算独占碱基数（$\text{cov}[\text{pos}] = 1$）
   - **实例级过滤**：独占碱基 $< \text{aligned\_length} / 4$（< 25%）则跳过该实例（`refine.c:1565`）
   - 计算独占节省：$\sum_i (2 \cdot a_i^{\text{excl}} - C_{\text{instance}}(a_i^{\text{excl}}, m_i^{\text{excl}}))$
   - 若无实例通过 25% 过滤，或独占节省 $< \text{model\_cost}$，则**修剪**
   - 递减被修剪家族的覆盖计数

---

## 6. MDL 打分与文库筛选

### 6.1 迭代 R 收敛

（`mdl.c:172-201`）

1. 初始化 $R = \max(n_{\text{families}}, 2)$
2. 重复（最多 3 次迭代）：
   a. 以当前 $R$ 对所有家族打分
   b. 计算被接受的家族数：分数 $> 0$ 且实例 $\geq 2$
   c. $R_{\text{new}} = \max(\text{accepted}, 2)$
   d. 若 $R_{\text{new}} = R$：收敛。否则 $R \leftarrow R_{\text{new}}$。

**注意**：由于单实例代价目前不依赖于 $R$（§1.6），此循环在 1 次迭代后即收敛。

### 6.2 贪心筛选

（`mdl.c:203-223`）

按 MDL 分数降序排列家族。依次接受 $\text{mdl\_score} > 0$ 且 $n_{\text{instances}} \geq 2$ 的家族。

### 6.3 覆盖位图

（`mdl.c:225-260`）

$\lceil N/8 \rceil$ 字节的位图。对每个被接受家族的实例，标记基因组位置：

```c
BIT_SET(pos): covered[pos >> 3] |= (1 << (pos & 7))
BIT_GET(pos): covered[pos >> 3] & (1 << (pos & 7))
```

仅用于统计报告，不用于筛选决策。

### 6.4 最终描述长度

（`mdl.c:262-268`）

$$DL_{\text{library}} = L_{\text{int}}(R_{\text{final}}) + \sum_{r \in \text{accepted}} \text{model\_cost}_r$$

$$DL_{\text{total}} = 2N - \text{total\_savings} + DL_{\text{library}}$$

$$\text{compression\_ratio} = \frac{DL_{\text{total}}}{2N}$$

### 6.5 修剪后恢复轮（Post-Prune Recovery Pass）

（`main.c:998-1045`）

修剪减少 $R$ 后，以 $R_{\text{final}} = \max(\text{accepted\_after\_prune}, 2)$ 重新对所有被拒绝家族打分。MDL 分数现在为正的家族被重新接受。

---

## 7. 大基因组支持

### 7.1 基因组采样

**激活条件**：基因组长度 $> \text{sample\_size}$（默认 1 Gb）。（`main.c:765`）

**算法**（`main.c:181-274`）：

1. **瓦片划分**：将每条序列分为大小为 $w$（默认 1 Mb）的不重叠瓦片。
   $$\text{tiles\_per\_seq}[i] = \lfloor \text{seq\_len}_i / w \rfloor$$

2. **目标窗口数**：$n_w = \lfloor \text{sample\_size} / w \rfloor$，截断至 $[1, \text{total\_tiles}]$。

3. **部分 Fisher-Yates 洗牌**（`main.c:221-227`）：
   ```
   srand(seed)
   For i = 0 to n_w - 1:
       j = i + rand() % (total_tiles - i)
       swap(tiles[i], tiles[j])
   ```
   默认 seed = 42。

4. **按基因组坐标排序选中的瓦片**以获得顺序内存访问（`main.c:230`）。

5. **创建采样基因组**——通过 `genome_create_chunk()` 以选中的瓦片为 segment。

**坐标重映射**（`main.c:280-320`）：发现完成后，对每个实例：
1. 计算原始位置：$\text{raw} = \text{position} - \text{PADLENGTH}$
2. 二分搜索找到所属 segment $s$：$\text{raw} \geq \text{seg\_start}[s]$
3. 重映射：$\text{original} = \text{segments}[s].\text{raw\_start} + (\text{raw} - \text{seg\_start}[s])$
4. 更新位置：$\text{position} = \text{original} + \text{PADLENGTH}$
5. 更新 `seq_index` 为原始基因组的序列索引

**内存管理**：采样后释放原始基因组序列缓冲区；保留元数据（boundaries、sequence_ids）用于 BED 输出。（`main.c:783-788`）

### 7.2 分块发现（Chunked Discovery）

**激活条件**：基因组长度 $> \text{chunk\_size}$（默认 200 Mb）。（`main.c:854`）

**序列分段**（`main.c:377-436`）：

1. **分割阈值**：$\tau = 1.8 \times \text{chunk\_size}$

2. 若序列长度 $> \tau$：
   - 分段数：$p = 2^k$，其中 $k = \lceil \log_2(\text{seq\_len} / \tau) \rceil$
   - 分段基本大小：$\text{part\_base} = \text{seq\_len} / p$

3. **重叠**：相邻 segment 在每个边界重叠 $L$ 个碱基（默认 10000）。截断至序列边界。（`main.c:411-416`）

**LPT 装箱**（`main.c:438-487`）：

1. 按长度降序排列 segment（最长处理时间优先，Longest Processing Time first）
2. 箱数：$b = \lceil \text{total\_seg\_size} / \text{chunk\_size} \rceil$，限制在 $[2, n_{\text{segments}}]$
3. 对每个 segment（从最大开始）：分配到当前总量最小的箱
   $$\text{bin}^* = \arg\min_b \text{bin\_size}[b]$$

这最小化了 $\max_b(\text{bin\_size}[b])$，对并行挂钟时间最优。

**并行执行**（`main.c:501-550`）：
- 批大小：$\min(n_{\text{bins}}, n_{\text{threads}})$
- 每个工作线程通过 `genome_create_chunk()` 创建 chunk 基因组（复制序列数据，共享 sequence_id 字符串）
- 每个 chunk 独立调用 `discover_families()`，拥有自己的 `DiscoverContext`
- **每 chunk 的 l-mer 长度**：$l_{\text{chunk}} = \lceil 1 + \log_4(N_{\text{chunk}}) \rceil$ —— 比全基因组的 $l$ 更短，提高种子灵敏度

**结果合并**（`main.c:570-578`）：所有 chunk 结果合并为单个 `CandidateList`，家族 ID 重新编号。

### 7.3 组合方式

1. **采样**（若触发）：将 10+ Gb 基因组缩减至约 1 Gb
2. **分块发现**（若触发）：将约 1 Gb 基因组拆分为约 200 Mb 的并行 chunk
3. **精炼**：在完整（采样后的）基因组上运行，而非按 chunk
4. **输出**：若使用了采样，则坐标重映射回原始基因组

小基因组（$< \text{chunk\_size}$）：不分割，单线程运行。

---

## 8. 数据结构

### 8.1 核心类型（`types.h`）

| 类型 | 定义 | 用途 |
|------|------|------|
| `gpos_t` | `int64_t` | 基因组位置（支持 >2 Gb 基因组） |
| `glen_t` | `int64_t` | 基因组/序列长度 |
| `freq_t` | `int32_t` | K-mer 频率 |
| `uid_t` | `int32_t` | 家族 ID |

DNA 编码：A=0, C=1, G=2, T=3, N=99。

互补：$\overline{c} = 3 - c$（$c \neq 99$ 时）。

### 8.2 Genome（`genome.h`）

```
Genome {
    char    *sequence       // 数值 DNA 编码，起始有填充
    glen_t   length         // raw_length + PADLENGTH
    glen_t   raw_length     // 实际碱基数（无填充）
    gpos_t  *boundaries     // 累积序列边界（填充前坐标）
    int      num_sequences  // FASTA 记录数
    char   **sequence_ids   // 序列标识符
}
```

**边界约定**（$N$ 条序列）：
- `boundaries[i]` = 序列 $i$ 的末尾（= 序列 $i+1$ 的起始），原始坐标
- `boundaries[N-1]` = raw_length + 1（哨兵值）
- `boundaries[N]` = 0（终止符）
- 填充坐标：原始坐标 + PADLENGTH

### 8.3 KmerEntry（`kmer.h`）

```
KmerEntry {
    uint64_t  kmer            // 规范 k-mer，2 bit/碱基压缩
    freq_t    frequency       // 出现计数
    gpos_t    last_plus_occ   // 正向链上次出现位置（TANDEMDIST 用）
    gpos_t    last_minus_occ  // 反向互补链上次出现位置（TANDEMDIST 用）
    gpos_t   *positions       // 位置数组（符号 = 链方向）
    int32_t   num_positions   // 已存储位置数
    int32_t   cap_positions   // 已分配容量（最大 50000）
    KmerEntry *next           // 哈希链指针
}
```

### 8.4 Instance（`candidates.h`）

```
Instance {
    gpos_t position          // 填充基因组坐标中的起始位置
    glen_t aligned_length    // 比对长度
    int    cons_start        // 共识起始（0-based）
    int    cons_end          // 共识末尾（exclusive）
    int    num_edits         // 比对中的编辑数
    float  divergence        // 散度（0.0–1.0）
    int    score             // 比对分数
    int8_t strand            // +1 正向，-1 反向
    int    seq_index         // 所在 FASTA 序列
}
```

### 8.5 CandidateFamily（`candidates.h`）

```
CandidateFamily {
    uid_t    id                // 家族 ID
    char    *consensus         // 数值碱基（0/1/2/3）
    int      consensus_length
    int      component_id
    int      topology          // TOPO_LINEAR=0, TOPO_COMPLEX=1, TOPO_CYCLIC=2
    freq_t   estimated_copies
    Instance *instances
    int      num_instances
    int      cap_instances
    double   mdl_score         // >0 表示被接受
    double   model_cost        // L_int(len) + 2*len
}
```

### 8.6 MDLResult（`mdl.h`）

```
MDLResult {
    double total_model_cost     // 已接受家族模型代价之和
    double total_savings        // 已接受家族节省之和
    double dl_library           // L_int(R) + total_model_cost
    double dl_genome_given_lib  // 2*N - total_savings
    double dl_total             // dl_library + dl_genome_given_lib
    double compression_ratio    // dl_total / (2*N)
    int    num_accepted         // 已接受家族数
    int64_t bases_covered       // 实例覆盖的基因组碱基数
}
```

### 8.7 SeqSegment（`main.c`）

```
SeqSegment {
    int    seq_index      // 原始基因组序列中的索引
    gpos_t raw_start      // 原始基因组原始坐标中的起始
    gpos_t raw_end        // 原始基因组原始坐标中的末尾
    glen_t seg_length     // = raw_end - raw_start
}
```

---

## 9. 完整参数参考

### 9.1 命令行参数

#### 必需参数

| 参数 | 说明 |
|------|------|
| `-sequence <file>` | 输入 FASTA 基因组 |
| `-output <file>` | 输出重复文库（FASTA） |

#### 发现参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-l #` | 自动：$\lceil 1 + \log_4 N \rceil$ | L-mer 长度 |
| `-L #` | 10000 | 每侧最大扩展距离（bp） |
| `-minthresh #` | 2 | 种子最低 l-mer 频率 |
| `-goodlength #` | 30 | 共识最短长度预过滤 |
| `-maxgap #` | 5 | 发现 DP 带偏移量（MAXOFFSET） |
| `-match #` | 1 | DP 匹配分数 |
| `-mismatch #` | -1 | DP 错配分数 |
| `-gap #` | -5 | DP 空位罚分 |
| `-cappenalty #` | -20 | 退出比对的罚分上限 |
| `-minimprovement #` | 3 | 每步最小总分改进 |
| `-stopafter #` | 100 | 连续无进展 N 个位置后停止（WHEN_TO_STOP） |
| `-maxentropy #` | -0.70 | Shannon 熵过滤（自然对数，负值） |
| `-tandemdist #` | 500 | 同链 l-mer 间最小距离 |
| `-maxoccurrences #` | 10000 | 每个种子最大出现数（MAXN） |
| `-maxrepeats #` | 100000 | 最大发现家族数（MAXR） |
| `-freq <file>` | — | 预计算 l-mer 频率表 |
| `-freq-output <file>` | — | 输出 l-mer 频率表 |

#### 精炼参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-threads #` | 1 | 线程数 |
| `-mdl-mode <mode>` | exact | 位置编码模式：`none`、`exact`、`upper` |
| `-max-divergence #` | 0.30 | 实例接受的最大替换率 |
| `-refine-gap #` | -5 | 精炼空位罚分（高 indel 物种推荐 -3） |
| `-refine-maxoffset #` | 12 | 精炼 DP 带半宽（最大 32） |
| `-max-dp-cells #` | 10000000 | 合并比对最大 DP 单元数（约 40 MB） |

#### 大基因组参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-chunk-size #` | 200 | 分块大小（Mb，最小 10） |
| `-sample-size #` | 1000 | 基因组采样阈值（Mb，最小 100） |
| `-window-size #` | 1000 | 采样瓦片大小（kb，范围 [100, 10000]） |
| `-seed #` | 42 | 采样随机种子 |
| `-sample-output <file>` | — | 输出采样基因组 FASTA |

#### 输出参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-instances <file>` | — | 实例 BED6 输出 |
| `-stats <file>` | — | 家族统计 TSV |
| `-v` / `-vv` | 关闭 | 详细程度（1 或 2） |

### 9.2 内部常量

| 常量 | 值 | 来源 | 说明 |
|------|-----|------|------|
| PADLENGTH | 11000 | `types.h:22` | 基因组填充（$\geq L$） |
| DEFAULT_TANDEMDIST | 500 | `types.h:23` | 串联距离 |
| DEFAULT_MAXN | 10000 | `types.h:24` | 每候选家族最大实例数 |
| HASH_SIZE | 16000057 | `discover.c` | 发现哈希表大小（素数） |
| LOG2_C0 | 1.5179605508986484 | `mdl.c:16` | Rissanen 常数 |
| KMER_MAX_POSITIONS | 50000 | `kmer.h:24` | 位置数组上限 |
| NUM_STRIPES | 4096 | `kmer.c:131` | 并行计数锁条带数 |
| POOL_BLOCK_SIZE | 4096 | `kmer.c:10` | KmerEntry 内存池块大小 |
| REFINE_SCREEN_K | 8 | `refine.h:9` | K-mer 筛选用 k 值 |
| PROFILE_BITS | 65536 | `refine.c:18` | $4^8$ profile 位集 |
| REFINE_MIN_JACCARD | 0.15 | `refine.h:10` | 合并预筛选阈值 |
| REFINE_MATCH | +2 | `refine.h:23` | 合并比对匹配分 |
| REFINE_MISMATCH | -3 | `refine.h:24` | 合并比对错配分 |
| REFINE_GAP | -2 | `refine.h:25` | 合并比对空位分 |
| REFINE_MIN_IDENTITY | 0.80 | `refine.h:13` | 80% 一致性 |
| REFINE_MIN_COVERAGE | 0.80 | `refine.h:14` | 80% 覆盖度 |
| REFINE_MIN_ALIGNED | 80 | `refine.h:15` | 80 bp 最低比对长度 |
| REFINE_RELAXED_IDENTITY | 0.70 | `refine.h:19` | 宽松合并一致性 |
| REFINE_RELAXED_COVERAGE | 0.70 | `refine.h:20` | 宽松合并覆盖度 |
| REFINE_OVERLAP_RELAX | 0.50 | `refine.h:18` | 宽松合并实例重叠 |
| REFINE_BIMODALITY_THRESH | 0.20 | `refine.h:35` | 拆分双峰性最低值 |
| REFINE_MIN_SPLIT_INSTANCES | 10 | `refine.h:32` | 拆分最少实例数 |
| REFINE_MIN_CLUSTER_SIZE | 3 | `refine.h:33` | 子家族最少实例数 |
| REFINE_MIN_DIV_GAP | 0.05 | `refine.h:34` | 最小散度差 |
| REFINE_DIV_BINS | 100 | `refine.h:36` | Otsu 直方图 bin 数 |
| ALIGN_MATCH | +1 | `align.h:10` | 精炼匹配分 |
| ALIGN_MISMATCH | -1 | `align.h:11` | 精炼错配分 |
| ALIGN_CAPPENALTY | -20 | `align.h:12` | 分数下限 |
| ALIGN_WHEN_TO_STOP | 100 | `align.h:13` | 停滞阈值 |
| ALIGN_MAX_ITERATIONS | 10 | `align.h:14` | 收敛限制 |
| ALIGN_MAXOFFSET_LIMIT | 32 | `align.h:16` | `-refine-maxoffset` 最大值 |
| ALIGN_MAX_EXTENSION | 10000 | `align.h:17` | 最大扩展 bp |
| EXTENSION_SLACK | 15 | `align.h:15` | 距边缘最小碱基数 |
| MAX_SEED_HITS | 50000 | `align.c:18` | 种子命中上限 |
| MAX_CONS_KMERS | 10000 | `align.c:19` | 共识 k-mer 上限 |
| DISCOVER_SPLIT_THRESHOLD | 200000000 | `main.c:34` | 默认分块大小（200 Mb） |
| 片段组装距离 D | [500, 5000] | `refine.c:1719-1721` | 近邻距离 |
| 片段共现次数 | $\geq 3$ | `refine.c:1834` | 最小共现数 |
| 片段嵌套保护 | 0.50 | `refine.c:1845` | 包含阈值 |
| 片段大小比例 | $\geq 0.10$ | `refine.c:1859` | 最小相对大小 |
| 片段方向一致性 | $\geq 0.80$ | `refine.c:1863` | 同向比例 |
| 片段 MAD | $\leq 100$ | `refine.c:1996` | 间隔一致性 |
| 片段间隔范围 | $[-20, \text{median}+2\text{MAD}]$ | `refine.c:1993-1996` | 允许间隔 |
| 片段链长度上限 | 10 | `refine.c:1925` | 最大并查集链长 |
| 修剪独占过滤 | 25% | `refine.c:1565` | 每实例最低独占碱基 |

---

## 10. 复杂度分析

### 10.1 发现阶段

| 操作 | 复杂度 | 说明 |
|------|--------|------|
| L-mer 计数 | $O(N)$ | 单遍基因组扫描 |
| 单家族扩展 | $O(\text{WHEN\_TO\_STOP} \cdot N_{\text{occ}} \cdot (2\text{MAXOFFSET}+1)^2)$ | 带状 DP |
| 掩蔽 | $O(N_{\text{occ}} \cdot L_c)$ | 每个出现的 1-vs-1 DP |
| 总发现 | $O(R \cdot (N + N_{\text{occ}} \cdot L_c))$ | $R$ 个家族，摊销掩蔽 |

### 10.2 K-mer 表

| 操作 | 复杂度 | 说明 |
|------|--------|------|
| 计数 | 每线程 $O(N)$，$T$ 线程 $O(N/T)$ | 条带锁 |
| 位置索引 | $O(N)$ 计数 + $O(N)$ 填充 | 两阶段 |
| 查找 | 期望 $O(1)$ | 链式法，低装载因子 |
| 修剪 | $O(T_{\text{size}})$ | 单遍扫描表 |

### 10.3 种子比对

| 操作 | 复杂度 | 说明 |
|------|--------|------|
| 种子扫描 | $O(\sum f_j)$ | 每家族，通过位置索引 |
| 聚类 | $O(H \log H)$ | 排序 $H$ 个种子命中 |
| 带状 DP | 每个锚点 $O(L_c \cdot W)$ | $W = 2\text{MAXOFFSET}+1 = 25$ |
| 实例收集 | 每家族 $O(A \cdot L_c \cdot W)$ | $A$ 个锚点 |
| 重建 | $O(n_{\text{inst}} \cdot L_c)$ | 每次迭代 |

### 10.4 精炼

| 操作 | 复杂度 | 说明 |
|------|--------|------|
| K-mer Jaccard | $O(R^2 \cdot P)$ | $P = 1024$ 字 popcount |
| 合并比对 | 最坏 $O(R^2 \cdot L_c^2)$ | DP 单元限制防止最坏情况 |
| Otsu 拆分 | $O(R \cdot n_{\text{inst}} + R \cdot B)$ | $B = 100$ 个 bin |
| 片段组装 | $O(I \log I + I \cdot D_{\text{avg}})$ | $I$ = 总实例数，扫描线 |
| 修剪 | $O(R \cdot n_{\text{inst}} + N)$ | 覆盖数组 |

### 10.5 大基因组

| 操作 | 复杂度 | 说明 |
|------|--------|------|
| 采样 | $O(n_w)$ | Fisher-Yates 部分洗牌 |
| 分块分割 | $O(S)$ | $S$ = 序列数 |
| LPT 装箱 | $O(n_{\text{seg}} \cdot n_{\text{bins}})$ | 简单循环 |
| 并行发现 | $O(\text{largest\_chunk} / T)$ | LPT 平衡 |
| 坐标重映射 | $O(I)$ | 按实例 |

### 10.6 性能分析数据（human3M，2026 年 3 月）

| 函数 | 运行时占比 | 角色 |
|------|-----------|------|
| `align_collect_instances` | 59.4% | 种子比对（主要瓶颈） |
| `semiglobal_align` | 10.8% | 合并阶段比对 |
| `cmp_seed_hits` | 9.6% | 种子命中排序 |
| `candidates_extract` | 1.6% | 图组分提取 |

---

## 输出格式

### FASTA（`-output`）

```
>R=<id> length=<consensus_length> copies=<num_instances> mdl=<mdl_score>
<共识序列 ASCII，每行 80 字符>
```

仅输出 `mdl_score > 0` 的家族。（`output.c:7-35`）

### BED6（`-instances`）

```
<chr>  <local_start>  <local_end>  R=<id>  <score>  <strand>
```

- 坐标：染色体局部坐标（raw_position - chr_offset - PADLENGTH）
- 分数：$\lfloor 1000 \times (1 - \text{divergence}) \rfloor$，截断至 [0, 1000]
- 仅输出已接受家族的实例。（`output.c:37-88`）

### TSV（`-stats`）

```
family_id  consensus_length  num_instances  divergence_mean  mdl_score  model_cost  topology
```

输出所有家族（不过滤）。Topology：`linear`、`cyclic` 或 `complex`。（`output.c:90-123`）

---

## 参考文献

1. Rissanen, J. (1978). Modeling by shortest data description. *Automatica*, 14(5), 465-471.
2. Price, A.L., Jones, N.C., & Pevzner, P.A. (2005). De novo identification of repeat families in large genomes. *Bioinformatics*, 21(suppl_1), i351-i358.
3. Wicker, T. et al. (2007). A unified classification system for eukaryotic transposable elements. *Nature Reviews Genetics*, 8(12), 973-982.
4. Otsu, N. (1979). A threshold selection method from gray-level histograms. *IEEE Trans. SMC*, 9(1), 62-66.
5. Grumbach, S., & Tahi, F. (1994). A new challenge for compression algorithms: genetic sequences. *Information Processing & Management*, 30(6), 875-886.
