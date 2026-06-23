# 通用 de novo 重复序列检测软件 — 方法设计

> 前瞻性方法设计提案（mdl-repeat 后继方向）。所有性能/灵敏度表述均为设计预期，
> 非已实现或已测基准；标注「预期」处需后续实验验证。本文件不依赖、也不修改
> 现有 mdl-repeat 代码，仅作设计参考。

## 0. 设计目标与“漏检”定义

**头条指标**：家族级 recall（cd-hit-est 聚类 + BLAST 命中 truth 库，RepeatModeler2/EDTA 标准）。辅助：家族级 precision（嵌合/假家族率）、consensus 边界完整度（与 truth 全长之比，按 80×80 / 90×80 两档）、大基因组的 RSS/wall-time。

**系统性漏检的五类目标**（决定算法取舍）：

| 类别 | 常规工具失败根因 | 本设计的主攻机制 |
|---|---|---|
| C1 高分歧/古老 TE（>20–30% div） | 拷贝间不再共享精确 l-mer，种子频率塌到阈值下 | 非精确种子（spaced/strobemer/syncmer）+ 蛋白域锚定 |
| C2 低拷贝（2–3 拷贝） | MINTHRESH≥2、贪心高频选种把它们排在长尾外；2-copy 难与偶然区分 | 草图全局配对（不靠频率排序）+ 低拷贝显著性检验 |
| C3 碎片化/被截断 | consensus 重建修剪边界；独立发现的片段不连 | 重复图 + 共现路径连接 + profile 边界界定 |
| C4 嵌套（TE-in-TE） | seed-and-extend 对嵌套不可知，越界或被切碎 | 图式“插入模式”检测 + 宿主/插入分解 |
| C5 非自主 MITE/SINE（短、TIR/无 ORF） | 短于检出长度、无蛋白域、易并入低复杂度过滤 | 结构特征种子（TIR/TSD/末端保守）+ 短元件专用招募 |

**关键假设**：① 无参考库；② recall（尤其 C1–C5）优先于极致 precision，但假家族需可控并可在分类阶段过滤；③ 允许可选外部工具但必须有纯内置回退（保证零依赖可运行）；④ 串联/卫星单独成模块（默认不混入 interspersed 库，但可输出）。

---

## 1. 总体策略

三根支柱，叠加“按规模分层”的执行调度：

- **支柱 A — 多模态灵敏种子**：不再单一依赖“高频精确 l-mer”。并联四条种子通道（精确 minimizer / 非精确 spaced·strobemer / 蛋白域 / 结构特征），各自命中不同漏检类别，统一汇入候选-配对集合。**这是提升 C1/C2/C5 recall 的核心。**
- **支柱 B — 重复图（repeat graph）统一表示**：把“候选重复实例”作为节点、相似/重叠/共现作为边，构建一张图；家族 = 图中稠密子图，嵌套/碎片 = 图上的拓扑模式。用图统一解决 C3/C4，避免 seed-and-extend 的边界与嵌套盲区。
- **支柱 C — 原则化统计模型**：以可正确反映库规模的 MDL/MML 两部分码做家族取舍 + EM 迭代精修成员与边界 + 对低拷贝做显著性检验。修正 mdl-repeat 中“代价 R 无关、单趟、无 recovery”的退化。

**按规模分层**（输入自适应，无人工切换）：

```
小 (≤ ~50 Mb): 全基因组精确自比 + 全索引 + 全图   → 最高灵敏
中 (~50 Mb–2 Gb): minimizer 索引 + 分块发现 + 块间图归并
大 (> ~2 Gb): 草图流式 + 分层覆盖保证(非随机采样) + 增量图
```

与 mdl-repeat 的随机 tile 采样不同：大基因组**不做有损随机采样**（那会丢稀有/低拷贝家族 = 恶化 C2），改用“草图覆盖保证”的分层抽取（§2.6）。

---

## 2. 算法层详解（重点）

### 2.1 多模态种子与索引

四条并联通道，输出统一为 `SeedHit{seq_id, pos, strand, channel, anchor_key}`。

**M1 精确 minimizer（高拷贝快路，~RepeatScout 等价但更省内存）**
- (w,k)-minimizer 索引；canonical（fwd/rc 同桶）。仅保留出现次数 ≥2 的 minimizer 作为 anchor。
- 用途：高拷贝低分歧家族的快速锚定；为支柱 B 提供初始稠密节点。
- 自适应：`k = clamp(ceil(1+log4(N)), 15, 27)`，`w ≈ k`（大基因组增大 w 降密度）。

**M2 非精确种子（攻 C1 高分歧）— 三选一/可并用**
- **Spaced seeds**：多组（如 3–4 个）高灵敏 spaced pattern（weight≈12, length≈18–22），对单点替换鲁棒。比同权连续种子在 25–35% 分歧下命中率显著更高（预期）。
- **Strobemers**（randstrobe，n=2/3）：对 indel 与重排鲁棒，弥补 spaced seed 对 gap 的弱点。
- **Syncmers**（closed syncmer）：相比 minimizer 的“窗口最小”更稳定，分歧下锚点保留率更高。
- 实现要点：M2 仅在 M1 覆盖不足的区域触发（条件式，控成本）——即“M1 已稠密锚定处不重复花 M2 预算”。

```
seed_region(region):
    hits = minimizer_hits(region)            # M1
    if local_anchor_density(hits) < τ_dense: # 仅低密度区升灵敏
        hits += spaced_seed_hits(region)     # M2a
        hits += strobemer_hits(region)       # M2b
    return hits
```

**M3 蛋白结构域锚定（攻 C1 古老编码 TE，常规 DNA 工具完全无能为力）**
- 6-frame 翻译后对**内置精简 TE 蛋白域库**（RT/RH/IN/GAG/AP/ENV、Tnp/转座酶、Helitron Rep/Hel 等，源自 REXdb/Pfam/Gypsy DB 的 HMM 或代表序列）做灵敏搜索。
- 工具：`MMseqs2`（profile/敏感模式，`-s 7`）或 `hmmsearch`；命中给出 anchor + 元件类型先验（直接喂 §2.5 的家族先验与 §5 分类）。
- 价值：DNA 已分歧到 seed 失效的古老 LTR/LINE，其蛋白域仍保守 → 蛋白锚把 DNA 上找不回的家族“拉”出来。这是 RepeatScout/mdl-repeat 结构性缺失的一环。

**M4 结构特征种子（攻 C5 MITE/SINE + 串联分离）**
- TIR/TSD 扫描（反向重复末端 + 侧翼短正向重复）招募 MITE/DNA-TE 末端；polyA/末端保守招募 SINE。
- 串联/低复杂度：内置 TRF-like（或调 `TRF`/`ULTRA`）+ Shannon 熵 + 周期扫描，将串联/卫星**分流到独立模块**而非丢弃（保留 HOR 检测能力）。

**索引层**：
- 小/中基因组：minimizer 索引 +（可选）FM-index/后缀数组支持精确自比定位。
- 大基因组：minimizer 索引分块 + **MinHash 草图 DB**（§2.2）。
- 全程 64-bit 坐标（`int64`），前/后 padding，按记录边界（seq_index）隔离（沿用 mdl-repeat 的正确性教训：跨记录位置比较必须带 seq_index 守卫）。

### 2.2 草图预聚类（把 all-vs-all 从 O(n²) 降到可控）

候选重复区域可能上百万段；直接两两比对不可行。先用 **MinHash/MinHash-bottom-k（`mash` 或内置）** 把候选段粗聚类：

```
for each candidate region r:  sketch[r] = bottom_k_minhash(r)
buckets = LSH_band(sketch)             # 局部敏感哈希分桶
for each bucket B:                     # 只在桶内做精确/半全局比对
    pairs += verify_pairs(B)           # §2.3 的边
```
- 复杂度：从 O(n²) 降到 ≈ O(n·b)（b=桶平均大小）。
- 对 C2 低拷贝友好：草图配对**不依赖频率排序**，2 拷贝只要互相进同桶即被配上，不会因“频率低”被贪心淘汰。

### 2.3 重复图组装与嵌套/碎片消解（攻 C3/C4 — 本设计相对 mdl-repeat 的最大增益）

**建图**：节点 = 候选实例区间；边 = 桶内验证过的相似/重叠关系（带 identity、coverage、对齐坐标、链向）。
- 桶内验证比对：小用 `edlib`/banded WFA；中大用 `wfmash`（WFA + mashmap，对分歧鲁棒）或 `minimap2 -DP`（自比对、`-X` 去自身）。

**家族 = 稠密子图**：在相似图上做社区发现（连通分量 + Louvain/标签传播细分），每个社区 → 一个候选家族的实例集合。比 union-find“传递合并”更抗链式过并（mdl-repeat 的已知风险）。

**嵌套分解（C4）**：在图上识别“插入模式”——实例 X 的中段比对到家族 F_ins、而其两翼比对到家族 F_host：

```
detect_nesting(instance X):
    segs = chain_alignments(X)                  # X 被切成若干对齐段
    if segs = [host_L | ins | host_R] 且
       host_L,host_R ∈ 同一家族 F_host 且 ins ∈ F_ins:
        split X into host(L+R, 记插入断点) 和 ins
        在 F_host 实例上把插入位点记为可变间隔(不计入分歧)
```
- 效果：宿主家族不再因内部插入被切碎或撑长，插入元件独立成家族 → 同时救回 C3 和 C4。

**碎片连接（C3）**：对未被嵌套解释、但**空间共现 + 图上可路径相连**的相邻片段，沿基因组顺序在图中找一致路径并拼接（共现统计 + 方向一致 + 间隔分布 sanity，参照但强于 mdl-repeat 的 sweep-line 共现：这里有图路径证据而非仅邻近计数）。拼接仅在统计模型（§2.5）判定“拼后编码更短”时接受。

### 2.4 一致序列构建（避免边界修剪 = 改善 90×80 完整度）

- 家族实例集 → **POA（partial-order alignment）**一致序列：用 `abPOA`（SIMD，快）或 `spoa`。POA 天然处理 indel 与子结构，优于逐列贪心 majority（mdl-repeat 在 90% identity 下边界被修剪的根因）。
- **边界由 profile 概率界定**：在一致序列两端用 per-column 信息含量/覆盖陡降点定边界，而非“分数不再提升即停”——保留末端 TIR/LTR/polyA 等低信息但生物学关键的边界。
- 输出每家族：consensus + per-column profile（供 §2.5 打分与下游招募）+ 拓扑标记（linear/nested-host/tandem）。

### 2.5 家族统计模型与停止准则（修正 mdl-repeat 的 MDL 退化）

**两部分码 MDL/MML（R-依赖，恢复原则化库规模控制）**：

```
DL(Genome) = DL(Library) + DL(Genome | Library)
DL(Library) = Σ_f [ L_int(len_f) + 2·len_f ]        # 各 consensus
DL(Genome|Library) = Σ_f Σ_i  cost(instance i | consensus_f)
cost(i) = L_int(a_i) + L_int(m_i+1) + log2 C(a_i, m_i)   # 对齐+编辑(exact 模式)
        + [ pointer_term: ceil(log2 R) ]   # ★ 显式 R-依赖：选库时家族 id 开销
        + strand_bit
```
- ★ 关键修正：**保留 `log2 R` 家族指针项**，使每实例代价随库规模 R 变化 → “多收一个家族”的边际成本真实存在 → 停止准则有意义。这正是 mdl-repeat 当前实现丢掉、导致 R-收敛/recovery 全失效的点。
- 选择：family 按边际收益排序，**迭代 EM 直到 R 收敛**（外层 2–3 轮）：
  - E 步：固定当前库，把基因组每段分配到最佳家族（含“literal/不属任何家族”选项）；
  - M 步：重估每家族 consensus（POA）、成员、分歧、边界；重算 DL；接受净收益 > 0 的家族。
- **覆盖记账用 sweep-line 区间**（O(实例数)，非 O(基因组长) bitmap）——沿用 mdl-repeat 已验证的可扩展做法。
- **prune + recovery**：剪掉独占覆盖不足者；R 下降后边际指针成本下降，**重跑 recovery 把先前被拒家族再评一次**（在 R-依赖代价下此步**不再 inert**，与 mdl-repeat 形成对比）。

**低拷贝显著性检验（攻 C2，区分真 2-copy vs 偶然）**：
- 对拷贝数 < 阈值（如 <4）的候选，做 null 模型检验：在与候选同 GC/同长度的随机/打乱序列下，估计“两段达到该 identity×length 比对”的期望次数（类似 BLAST E-value 思路或基于 k-mer 组成的解析近似）。
- 仅当观测显著优于 null（如 E < 1e-3）才接纳为家族。**用统计显著性替代“频率阈值”做低拷贝把关**，既救 recall 又控假阳。

### 2.6 多尺度与大基因组的无损/分层处理

- **小**：全 all-vs-all（草图桶内）+ 全图 + 全 EM。
- **中**：minimizer 分块发现 → 每块本地图 → **块间用草图归并**同家族（跨块 consensus 对齐合并），避免边界家族被切。
- **大（多 Gb）**：
  - **不做随机有损采样**。改“草图覆盖保证”：先对全基因组建轻量草图，按 minimizer/草图把序列空间划成相似性区块，**保证每个相似区块至少有代表性子集进入精算**（稀有家族因此不会被随机丢——直接针对 mdl-repeat 采样恶化 C2 的缺陷）。
  - 流式 + 增量图：分块算草图 → 增量并入全局相似图 → 全局社区发现 → 代表性实例做 POA。内存与家族数成正比，与基因组长度近似解耦。
  - 复杂度目标（预期）：时间 ≈ O(N) 草图 + O(候选²/桶) 验证；峰值内存 ≈ O(家族数 × consensus 长 + 实例索引)。

---

## 3. 系统架构

```
mdl-rep2/
  core (C/C++ 或 Rust)                # 性能关键路径
  ├─ io        FASTA 流式读、64-bit、记录边界、padding
  ├─ index     minimizer / 草图DB / (可选)FM-index
  ├─ seed      M1–M4 多模态种子(插件化通道)
  ├─ sketch    MinHash + LSH 预聚类
  ├─ graph     候选实例图、社区发现、嵌套/碎片拓扑
  ├─ consensus POA 封装(abPOA/spoa) + profile 边界
  ├─ model     MDL/MML + EM + 低拷贝显著性 + sweep-line 覆盖
  ├─ scale     分层调度(小/中/大)、分块、块间归并
  ├─ toolrun   外部工具 fork/exec + 超时 + 回退(参照 mdl-repeat tool_runner)
  └─ output    FASTA 库 / BED 实例 / TSV 统计 / 图(GFA, 供审查)
  orchestrator (主驱动，可 CLI/配置文件)
```
- **并行**：种子/草图/桶验证 = 数据并行（线程池 + 分块）；EM 的 E 步 = 实例并行；图社区发现 = 分块 + 归并。条带锁哈希（参照 mdl-repeat kmer.c）。
- **可复现**：固定随机种子（草图 hash、LSH）、记录工具版本与参数、emit run manifest（命令/环境/校验和）。
- **可扩展**：种子通道与外部工具均插件化（接口统一），便于关掉蛋白通道做纯 DNA 模式或增配新通道。

---

## 4. 外部工具调用（条件式 + 回退）

| 步骤 | 默认调用 | 触发条件 | 选它的理由 | 纯内置回退 |
|---|---|---|---|---|
| 桶内自比对 | `wfmash` 或 `minimap2 -X -DP` | 中/大基因组、候选量大 | WFA 对分歧鲁棒、mashmap 预过滤快 | banded WFA / `edlib`（内置） |
| 草图聚类 | `mash`（或内置 MinHash） | 候选段 > 阈值 | 成熟、快 | 内置 bottom-k MinHash + LSH |
| 蛋白域锚定 (M3) | `MMseqs2`（敏感）/ `hmmsearch` | 启用蛋白通道（默认开，可关） | 敏感、快；HMM 域库现成 | 跳过 M3（仅 DNA 模式，C1 灵敏度下降，显式告警） |
| 一致序列 | `abPOA` / `spoa` | 每家族成员 ≥3 | SIMD POA，质量高于贪心 majority | 内置 banded POA（慢，小家族用） |
| 串联检测 (M4) | `TRF` / `ULTRA` | 高周期/低复杂度区 | 串联识别标准件 | 内置周期扫描 + 熵过滤 |
| 库去冗余 | `cd-hit-est` | 出库前 | 评测口径一致(家族级) | 内置 80-80-80 + 图社区 |
| 终库 QC | `seqkit stats` | `-qc` 开 | 非破坏统计 | 内置统计 |
| (可选)分类 | `RepeatClassifier`/`DeepTE` | `-classify` 开 | TE 家族命名 | 标 Unknown + 蛋白域先验 |

- **铁律**：任何外部工具缺失都有内置回退、绝不静默降级——缺失即显式告警并记录“因 X 不可用，回退为 Y，影响 Z”。`-external-tools require` 模式下缺失为硬错（参照 mdl-repeat external_qc 策略）。
- 工具定位走 PATH 或显式 `--toolname-path`，**不硬编码 conda 路径**（mdl-repeat 的 `find_rmblastn` 硬编码是反例）。

---

## 5. 漏检类别 → 捕获机制 对照表 + 风险/验证

**捕获对照**：

| 漏检类别 | 主捕获机制 | 兜底机制 |
|---|---|---|
| C1 高分歧/古老 | M2 非精确种子 + M3 蛋白域锚定 | wfmash 分歧鲁棒比对 |
| C2 低拷贝 2–3 | §2.2 草图配对(不靠频率) + §2.5 显著性检验 | 图社区最小簇 |
| C3 碎片化/截断 | §2.3 图路径碎片连接 + §2.4 POA 边界 | 共现 sweep-line |
| C4 嵌套 | §2.3 插入模式分解 | 嵌套 veto(防错并) |
| C5 MITE/SINE 短非自主 | M4 TIR/TSD/polyA 结构种子 | 短元件专用招募 + profile 边界 |

**主要风险与缓解**：
- **假家族/嵌合**（多模态种子更激进 → precision 压力）：EM + MDL 净收益门 + 低拷贝显著性 + 图社区（非链式合并）多重把关；可选分类阶段过滤 Unknown 噪声。
- **蛋白通道引入非 TE 基因家族**（如把宿主基因当重复）：限定 TE 专属域库 + 要求 DNA 层也成簇（蛋白锚仅作 anchor，不单独成家族）。
- **大基因组成本**：草图 + LSH 控 all-vs-all；蛋白通道与全图仅在中小或分块内开；提供 `--fast` 档关闭 M2/M3。
- **串联误入 interspersed 库**：M4 显式分流 + 周期/熵过滤。

**验证与基准方案**：
1. 黄金标准基因组：TAIR10、果蝇、人 chr（有人工 curated 库）+ 模拟基因组（已知 TE，可注入 C1–C5 难例，控制分歧/拷贝/嵌套）。
2. 指标：家族级 recall/precision（cd-hit+BLAST，80×80 与 90×80 两档）+ 边界完整度 + 大基因组 RSS/wall。
3. **难例消融**：分别关闭 M2/M3/M4 与图嵌套模块，量化每机制对 C1–C5 的边际贡献（证明各通道确有捕获其目标类别，而非堆叠）。
4. 横向对比 RepeatScout / RepeatModeler2 / EDTA / mdl-repeat，**重点报告“它们漏检而本工具检出”的家族集**（即设计目标的直接证据）。
5. 阴性对照：随机/打乱序列输入应产出近空库（验证 §2.5 显著性把关有效）。

---

## 6. 与 mdl-repeat 的关键差异（一句话定位）

本设计在 mdl-repeat 的 seed-and-extend+MDL 骨架上做四处结构性升级：**① 单一精确高频种子 → 四通道灵敏种子（含蛋白域）** 攻 C1/C2/C5；**② 线性 seed-and-extend → 重复图** 攻 C3/C4 嵌套碎片；**③ 退化的 R-无关单趟 MDL → R-依赖两部分码 + EM + 低拷贝显著性** 恢复原则化停止并控假阳；**④ 大基因组随机有损采样 → 草图覆盖保证的无损分层**，不再以牺牲稀有/低拷贝家族换取扩展性。同时保留其经验正确性资产：64-bit 坐标、seq_index 跨记录守卫、sweep-line 覆盖、tool_runner 式带回退的外部工具纪律。
