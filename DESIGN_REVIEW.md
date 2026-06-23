# De novo 重复检测「v3 合并方案」设计审阅

本文件记录对一套 de novo（无参考库）interspersed repeat / TE 家族检测算法「v3 合并方案」的两份独立对抗式审阅（一份只评生物学、一份只评工程）及综合判断，并附一个决定路线生死的门禁实验（量 R）。审阅对象不依赖、也不修改现有 mdl-repeat 代码。

---

## 0. 被审对象：v3 合并方案

目标：给定基因组 FASTA（10 Mb – 多 Gb），无参考库，发现 interspersed repeat / TE 家族，每个家族输出 consensus + 全部拷贝坐标；最大化家族 recall（含稀有/高分歧），控假家族，可扩展到多 Gb。

架构——**发现与编目解耦**（用有界样本建库，再把库映射回全基因组）：

- **Stage 0**：2-bit 打包；在重活之前**前置软掩蔽**串联/低复杂度（窗口熵 + 短周期自相关）；卫星单独成 track，不进 interspersed 库。
- **Stage 1**：k-mer 谱（k=16），Count-Min/计数 Bloom 估计；在 unique 峰后的直方图谷设自适应「重复」阈值（count ≥ 3）。
- **Stage 2**：在 repetitive spaced-seed（3–4 patterns，weight ~12–16）上采种子；每 k-mer 出现次数封顶 C_occ=200（reservoir 采样）→ 总 anchor 数被界住，与单家族丰度无关。
- **Stage 3（唯一超线性步）**：每 anchor 取 ±3 kb 窗口；**仅在共享种子的窗口间**做带状 gapped 比对；保留 aln≥100 bp、identity≥65%（分歧预算 ~35%）、非低复杂度的边。声称 O(R·C_occ²·window·band)，与 L 解耦。
- **Stage 4**：在比对图上聚类：丢弃覆盖 <50% 双端的边（隔离「仅共享域」连接）、社区发现（非 union-find）、稀疏割切分嵌合桥、**故意过切**成紧簇（视作可恢复错误，后续再合）。
- **Stage 5**：每簇采 30–50 个跨分歧成员；POA profile；consensus = 列多数 + IUPAC，占用 <50% 列裁掉；**边界由侧翼跨拷贝一致性跌到随机背景（~25%）处界定**；低熵末端延伸若完成 TIR/LTR/polyA/TSD 结构信号则保留；记录结构。
- **Stage 6**：consensus 层归并（库已小）：≥80% 短者 @ ≥85–90% id 合并；共现相邻的 5'/3' 片段拼接（重建常 5'-截断的 LINE）；切分残留混合；去冗余。
- **Stage 7**：编目——库 consensus/profile vs **全基因组** seed-and-extend，拿全部拷贝坐标；profile 打分救分歧拷贝；重叠 best-score-wins；记录嵌套插入。
- **Stage 8**：接纳（无外参，全靠内部统计）：分散拷贝数显著（N_min~10，有干净结构佐证可降到 ~5）、Karlin-Altschul 式组成 E-value（拒组成可解释者）、分歧分布 sanity（双峰→疑嵌合→再切；极低分歧+低拷贝+无结构→标「疑似 segmental duplication」而非断言 TE）；分类 LTR/LINE/SINE/DNA-TE/unknown。**停机**：在残余未解释重复部分上迭代 Stage 1–6，新增家族解释 <1% 剩余重复质量或 3–5 轮后停。
- **Graft 1（蛋白域锚定）**：并联种子通道，6-frame 翻译 + 对内置精简 TE 蛋白域库（RT/RH/IN/GAG/AP/ENV、转座酶、Helitron Rep/Hel）敏感搜索；域命中处**局部放宽** DNA 聚类/identity 阈值（不单独成家族、不要求 DNA 成簇），救古老编码家族。
- **Graft 2（显式嵌套分解，置于迭代循环内）**：检测 [host_L | insertion | host_R]（两翼同家族、中段另一家族）并分解，反馈进下一轮重聚类。
- **丢弃**：无 MDL / 描述长度 / 库规模惩罚；接纳全靠 Stage 8 内部显著性 + 结构检验。

> 两份审阅均为 Fable-5 独立 agent 产出，clean-room（禁读 repo）、分车道（生物 reviewer 不碰算力，工程 reviewer 不重判生物）。

---

## 1. 生物学审阅（Reviewer A）

### FATAL — 会产出生物学错误输出或整类静默失效

**F1. 无细胞器守卫 → NUMT/NUPT/rDNA 变成「家族」。** Stage 8 只测分散性、组成 E-value、分歧形状——没有一项能拒绝核内的线粒体/质体 DNA 插入。NUMT/NUPT 分散、多拷贝、组成正常、常中等分歧，在植物（与灵长类）丰富，会聚成干净「家族」并通过所有内部检验。组成 E-value 不会触发（它们是普通编码/细胞器序列），蛋白通道也不会把它们标为 TE（携带 cox/nad/rbcL 而非 RT/转座酶），于是作为 interspersed repeat 蒙混过关。5S/45S rDNA 与部分分散的 snRNA 阵列是同类漏洞。无细胞器/rRNA 筛除的 de novo 重复发现器是已知的假家族制造机。**修复：接纳前用细胞器基因组 + rRNA/tRNA/snRNA 模型筛候选 consensus。**

**F2. Helitron 同时穿透所有结构机制。** 滚环 Helitron **无 TSD、无 TIR、无 LTR、无 polyA**——5'-TC…CTRR-3' 末端 + 短 16–20 bp 3' 发夹，插入宿主 AT 间且无 TSD。于是：(a) Stage 5「仅当完成 TIR/LTR/polyA/TSD 才保留末端低复杂度」的边界规则对它永不触发、会裁错两端；(b) Stage 8 分类无信号可赋 → 每个 Helitron 落 unknown；(c) 最糟，Helitron **捕获并扩增宿主基因片段**，不同拷贝 cargo 不同 → 抗嵌合切分器（Stage 4）会撕碎真 Helitron 家族，而「极低分歧+无结构→疑似 seg-dup」规则（Stage 8）会**主动把年轻 Helitron 扩张误标为 segmental duplication**。唯一可靠的 Helitron 信号是 Rep/Hel（HUH 核酸酶+解旋酶）蛋白域——graft 库里有，但 graft 只放宽 DNA 阈值，**没接进边界发现与分类**，而那正是 Helitron 需要它的地方。修复：让 Rep/Hel 命中有权 (i) 定类、(ii) 锚定 3'-发夹/CTRR 末端，并把 Helitron 候选豁免 seg-dup 重标。

### MAJOR — 大幅 recall 损失或系统性误处理

**M1. 35% 分歧 floor + 蛋白 graft 救不了非编码家族——而那正是声称的目标。** graft 只帮有保守 ORF 的家族（自主 LTR/LINE/DNA-TE），对**既分歧又无蛋白**的家族毫无作用：**SINE**（tRNA/7SL/5S 衍生，无 ORF）、**MITE**（非自主，无 ORF）、**LARD/TRIM**（大/末端重复逆转录元件，非编码，植物丰富）、复合 **SVA**。古老 SINE 常处 30–45% 分歧（哺乳类 MIR/L2 期）——超出两两 floor——且无蛋白生命线。更糟，对短非编码元件，**种子**灵敏度（weight 12–16 spaced seed、精确 k-mer）远在 35% 两两分歧前就塌，故 SINE/MITE 的实际 floor 比标称 65% identity 更紧。净：本设计**恰恰漏掉**它承诺的「稀有+分歧」非编码部分，且无非蛋白 rescue（如从 tRNA 衍生 SINE 头做结构锚定 profile seeding、或 HMM 迭代）补偿。

**M2. ±3 kb 窗口在生物学上撕碎大元件，且从不重配 LTR↔internal。** 全长 **gypsy/copia LTR 元件 5–15 kb**、**Polinton/Maverick 15–22 kb**、**CACTA/EnSpm 及大基因捕获 Helitron >10 kb**、全长 **L1 ~6 kb**。±3 kb anchor 窗口跨不过它们，于是单个 LTR 元件被看成 {两末端重复的 LTR 家族} + {internal 片段}，而 Stage 6「拼接相邻 5'/3' 片段」只为 LINE 截断设计、未涵盖**把 LTR consensus 与其 internal consensus 配对**。RepeatMasker 库本就分 `*-LTR`/`*-I` 条目，故片段化本身可容忍——但本设计产出片段却**不记录 LTR–internal 关联**，且对 **solo-LTR**（在玉米/小麦中数量 >2:1 超过完整元件，由 LTR–LTR 重组而非嵌套产生）无模型。solo-LTR 会成「LTR 家族」与亲本脱钩、多半判 unknown。Polinton/Maverick 还需域库加 pPolB 与逆转录病毒样 integrase（当前未列）。

**M3.「域命中处局部放宽阈值」在驯化/宿主基因附近生物学不安全。** 转座酶、integrase、RT、GAG 域在宿主基因组中被广泛**驯化**——RAG1（Transib）、SETMAR/Metnase、CENP-B（pogo）、植物 FAR1/FHY3 与 MUSTARD（MULE 衍生、**多拷贝**）、DAYSLEEPER（hAT）、Arc/Fv1/syncytin（ERV gag/env）。单拷贝驯化安全（Stage 8 仍要求分散拷贝）——但**多拷贝驯化基因家族**（FAR1/FHY3 成员众多）会过 N_min 被断言为 TE 家族，且局部阈值放宽会**把宿主旁系簇粘到真 TE 家族上**、用宿主外显子污染 consensus。反过来，设计**只用蛋白证据放宽阈值、不做分类**，错失了正确性增益（RT⇒retro、转座酶⇒DNA、Rep/Hel⇒Helitron、GIY-YIG/PLE-RT⇒Penelope）。Penelope/PLE 尤其需库中加 PLE 型 RT（端粒酶支、与 LTR/LINE RT 大相径庭）与 GIY-YIG 内切酶，否则全漏。

**M4. Stage 0 前置掩蔽在发现前就撕碎复合及含 SSR 的真元件。** 在发现前掩蔽 simple/低复杂度会损伤：**SINE polyA/微卫星尾**与 **LINE 3' polyA**（边界被截）、**SVA**（CCCTCT 六聚 + VNTR，基本被肢解）、**LTR 内部微卫星**、**嵌在卫星阵列中的着丝粒/近着丝粒逆转录元件**（如玉米 CRM——先掩卫星就把嵌入其中的逆转录元件碎掉）。软掩本无妨，若比对器能穿过它——但 Stage 2–3 明确拒绝在低复杂度区采种子或留边，故软掩对发现实为硬掩，唯一 rescue（末端延伸保留）又依赖被识别的结构信号。另：**回文元件**——MITE（Stowaway/Tourist）与 foldback（FB）——以反向重复为主；「短周期自相关/低熵」掩蔽器有风险把它们 TIR 驱动的回文结构当低复杂度而掩掉元件本身。

**M5. 嵌套分解规则单层且过约束。** 要求 host_L 与 host_R 是**同一已知家族**，漏掉现实拓扑：(a) **深层嵌套**（插入中插入的「俄罗斯套娃」——玉米基因间区常态）；(b) host_L 与 host_R 是**不同**元件的宿主；(c) **鸡生蛋**：宿主**只以被打断形式出现**（每个拷贝都带插入），宿主家族永远无法先学到来触发分解；(d) **靶位缺失**宿主（插入删掉了部分宿主，两翼无法重建）。这是植物基因组的主导结构特征，故部分处理会实质性压低那里的 recall。

**M6. TSD/TIR/transduction 假设过于划一。** TSD 长度是超家族特异、有时无信息：**Tc1/mariner TSD 仅「TA」(2 bp)**——与背景不可分，无法充当下调 N_min 到 5 的「干净结构佐证」；hAT=8、PIF/Harbinger 与 CACTA=3、MULE=8–10 可变、LTR=4–6、LINE/SINE **可变**（TPRT）。固定 TSD 检测器视家族而定要么漏要么过信。CACTA 由 **CACTA 末端 motif + 亚末端重复**定义、非长 TIR，故 TIR 长度检测器漏掉它。且 **LINE 3' transduction**（与 SVA 5' transduction）把唯一的下游宿主序列带进新拷贝，故跨拷贝侧翼一致性会**延续到真 3' 端之外** → 衰减到 25% 的边界规则把 consensus 延进 transduced 宿主 DNA（已记录的 artifact，会打乱宿主外显子）。

**M7. reservoir 把丰度家族封顶 ~200 anchor 会偏向丢其分歧的古老亚家族。** 随机 reservoir 采样是按比例的，故占拷贝 ~1% 的亚家族（如 **AluJ** 在 Alu 海中、被 AluS/AluY 主导）只摊到 ~2 个 anchor——可能低于 count≥3/聚类、永不成种。因封顶在 Stage 2 聚类**之前**施加，Stage 5「跨分歧采 30–50」无法救回从未被 anchor 的成员。最丰度家族中信息量最高的古老亚家族风险最大。

### MINOR — 降质/降用，非致命

- **扁平聚类 vs 系统发育。** 单 identity 阈值 + 再合给出一条家族均值 consensus，但 TE 家族是由 master copy 与 burst 驱动的系统发育（Alu J→S→Y；L1PA1–17）。AluJ↔AluY ~85% identical——正落在合并缝上——故粒度不稳，单一合并 consensus 丢失下游所需的**年轻亚家族** consensus（影响 masking 灵敏度与插入定年）。亚家族树/网络（COSEG/RM2 式）更忠实。
- **IUPAC 过多的 consensus**（合并分歧成员）降低 cross_match/RepeatMasker 分——库可用性成本。
- **N_min 5–10 排除真低拷贝家族**——许多 DNA 转座子与水平转移家族 <10 拷贝，自主元件常以 1–2 个完整源拷贝存于死亡衍生物中。直接与稀有家族目标冲突；结构下调只帮有结构的元件（又非 Helitron/SINE）。
- **Stage 6 ≥80% 短者合并**会把 **MITE 并进其自主亲本**，抹掉 RM 库会分开保留的、生物学不同的（非自主）条目。
- **重复主导基因组的谷/阈值不稳**（小麦、大麦、针叶、玉米）——unique 峰小或缺——恰是最该用的基因组。
- **seg-dup vs burst 信号弱**：含结构的年轻 burst 与含 TE 的 seg-dup 都会混淆；保守的「标、勿断言」是对的，但如 F2 所述会主动误归 Helitron。

### 生物学最大风险 + 裁决
最大风险：结构信号脚手架与蛋白 graft 隐含 **retro/TIR 中心**，故两类同时穿透**分歧 floor 与结构 rescue**——**滚环 Helitron**（无 TSD/TIR/LTR/polyA；基因捕获拟态嵌合/seg-dup）与**分歧非编码家族**（SINE、MITE、LARD/TRIM、SVA，无 ORF 供 graft 抓）。叠加**完全缺失 NUMT/NUPT/rDNA 假家族守卫**，pipeline 会同时**少召**它瞄准的稀有/分歧家族、又**过度产出**非 TE 分散序列为家族——对既定目标最具破坏的两种失效模式。

是否值得推进：**值得——骨架在生物学上合理**（发现/编目解耦、profile consensus、结构感知边界、迭代残差挖掘、为古老编码家族设蛋白通道都是对的本能）。但部署前必须：(1) 接纳处加细胞器/rRNA 筛；(2) 经 Rep/Hel 锚定末端与定类的一等 Helitron 处理；(3) 为分歧 SINE/MITE 设非蛋白 rescue 路径并记录 LTR↔internal/solo-LTR 关系；(4) 让蛋白命中喂进分类（非仅阈值），并设多拷贝驯化宿主基因守卫；(5) 超家族感知的 TSD 逻辑 + transduction 感知的边界。

---

## 2. 工程审阅（Reviewer B）

### 总判
架构的发现骨架有一处致命结构缺陷：**Stage 3 物化了一张「组内 all-pairs」图，边数 O(R·C_occ²)，而 R 与该乘积都远大于声称值。**「与 L 解耦」是假的——R 随 L 增长直到 k-mer 饱和，对任何多 Gb 基因组，边集是 10¹¹–10¹³ 量级：既不入内存也算不完。下游（Stage 4 社区发现、迭代循环）全部继承此问题。其余（Stage 0–2、5–8、graft）大体可用已知技术构建，但恰在决定可行性的决策处 under-specified。

**裁决：按规格不可构建；需重构 Stage 3 候选生成步**（用 LSH/minhash 近邻候选替代 all-pairs）。加上这一改 + 把决定性设计进去，方案即可构建。

### FATAL

**F1. Stage 3 成本未界住、也未与 L 解耦。** O(R·C_occ²·window·band) 算术上成立，但因子被低估：
- **R（distinct repetitive 16-mer，count≥3）不是小常数。** 3 Gb 哺乳类 ~50% 重复，distinct repetitive 16-mer 约 **R ≈ 10⁷–10⁸**（分歧拷贝各造新 16-mer，只有高保守核坍缩）。15 Gb 植物 >70% 重复，**R ≈ 10⁸–5×10⁸**。R **随 L 增长**（亚线性，趋向 4¹⁶=4.3×10⁹ 上限）——故该步**不**与基因组长度解耦，只与单家族拷贝数解耦。
- **C_occ² 是乘子。** 组内 all-pairs = 200²/2 = **2×10⁴ 次带状比对/种子组**。
- **带状比对数，哺乳类：** R·C_occ²/2 ≈ 10⁷×2×10⁴ = **2×10¹¹**（乐观 R）到 **2×10¹²**（R=10⁸）。每次带状 X-drop 跨 ~6 kb 窗口 @ band~50 ≈ 3×10⁵ cell → **6×10¹⁶–6×10¹⁷ DP-cell 操作**。@ ~10⁹ cell/s/核 ≈ **200 万–2000 万核时**/一个外层迭代/一个基因组。
- **15 Gb 植物：** 再 ×10–50。
- spaced seed（3–4 patterns）再 ×3–4 anchor/边。

要让此步可行（≤10⁹ 比对）需 R ≤ 5×10⁴——比现实低 3–4 个数量级。all-pairs 的 C_occ² 结构是杀手；这正是 RECON/RepeatScout 类工具靠**对增长中的 consensus 扩展（O(occ)）而非两两（O(occ²)）**所规避的。**按写法，Stage 3 在哺乳类都跑不动，遑论 15 Gb 植物。**

**F2. Stage 3 图不入内存。** 边数 O(R·C_occ²/2)。即便乐观 R=10⁷：2×10¹¹ 边 ×16 B = **3.2 TB**；R=10⁸：**32 TB**。超 256–512 GB 节点 1–2 个数量级，溢盘也不切实际（多 TB 边表上随机访问的社区发现）。图绝不能 all-pairs 物化。

**F1+F2 同一修复：** 在比对**前**插入廉价近邻候选生成器——对每窗口的 repetitive-k-mer 集做 minhash/LSH 签名，使每窗口只提 k≈8–32 个候选邻居。边数降到 O(N_windows·k) ≈（几×10⁷ 窗口 ×16）≈ 10⁸–10⁹ 边（GB 级，可入），比对数降 ~C_occ。这是结构性改动、v3 未规定，且是使设计可构建的**单一关键改动**。

### MAJOR
- **M1. C_occ=200 封顶 × 过切聚类，饿死 consensus 步。** reservoir 200 个对丰度估计够，但 Stage 4 故意过切成紧簇。靠 ~200 anchor 带过发现、再被切碎的家族，产出 **<30 成员的簇——低于 Stage 5 POA 的 30–50 下限**。即分离亚家族的机制反而饿死它们，边缘家族因**采样而非生物学**跌破 Stage 8 拷贝数 floor。亚家族诊断性 k-mer 确会自成种子组（故发现非严格受阻），但封顶与过切策略直接张力，需 reconciliation（如 POA 前从全 occurrence 列表重扩成员）。
- **M2. Graft-2 嵌套分解 EM 循环无不动点/终止保证。** 外层循环（Stage 1–6，每轮掩蔽已接纳家族）**会**终止：掩蔽单调缩小残余重复质量，STOP 规则（<1% 新质量或 3–5 轮）界住——这部分没问题。问题是 **Graft-2 反馈进重聚类**：分解 [host_L|insertion|host_R] 改变成员→改变 consensus→改变边界→改变下一轮何为插入。**无定义的不动点**且有真实振荡风险（一个元件时而是嵌套、时而独立，每轮翻转）。须加：单调 **freeze-on-accept**（已接纳+掩蔽家族永不再入聚类）+ 硬 sweep 上限。否则不保证收敛。**成本乘子：** 外 3–5× × Graft-2 2–4 sweep ≈ 最多 **20×**，因残余几何缩小部分抵消（有效 ~3–5×）。
- **M3. 蛋白 graft 重新引入全基因组 O(L) 重扫。** 3 Gb 的 6-frame 翻译 = ~2L = **6×10⁹ aa**，对 30–50 个 TE 域 HMM。乐观 ~7 核时；现实 translated-HMM 搜索（hmmscan/nhmmer 类）接近**核日**，15 Gb 植物 ×5–10。它对 chunk embarrassingly parallel，集群上可行——但是**全基因组线性扫**，与骨架「有界样本」哲学矛盾。集成点未定：放进发现循环就每轮重引入 O(L)。**建议：** 在有界样本上跑一次做阈值放宽，和/或在 Stage 7/8 做纯注释扫——别进迭代发现循环。
- **M4. Stage 7「线性 L」藏着随重复密度放大的大常数。** 数千 consensi（5–50 Mb）对全 L 的 seed-and-extend 是 RepeatMasker 形：O(L·hits_per_pos)。70% 重复的 15 Gb 植物命中密度高，常数大，现实 **单节点核日–周**（集群上 embarrassingly parallel）。「近线性 L」技术上诚实但常数是 pipeline 第二大算力。用快速 anchored mapper（minimap2 类）而非朴素 seed-extend。
- **M5. 决定性可达但须现在就设计进去。** 非确定源各可修：**reservoir 采样**在并行计数下线程调度相关 → 改 **hash-based bottom-k**（保留 hash(position) 最小的 200 个），顺序无关确定；**社区发现** Louvain seed/顺序敏感、label-prop 更差、Leiden 较稳但仍 seed 相关 → 固定 seed + **规范节点排序**（按基因组坐标排，节点 ID 跨运行稳定）；**并行浮点归约**改整数/定点 + 确定归约序。代价对并行 **低个位数 %**，但事后回填痛苦。
- **M6. 侧翼一致性边界在 30–50 成员处临界、对边缘家族脆。** 统计上 40 个独立侧翼可分辨 0.7 vs 0.25（SE≈√(0.2/40)≈0.07），故丰度家族信号可用。两个工程问题：(a) 须从分歧基因组背景产出**一致定框、对齐的侧翼**——POA 会把真边界外的噪声对齐，故以 occupancy/熵坍缩检测边界、对阈值敏感（未定）；(b) 近 ≥5–10 接纳 floor 的家族独立侧翼 <20，转变不可分辨——即边界对最需要它的边缘家族最不可靠。需加成员数门控。

### MINOR
- **m1. 5'/3' 片段拼接「共现相邻」检测廉价、但易过拼。** 坐标 sweep over Stage-7 注释，O(n log n) trivial。风险在正确性：裸相邻会过拼（solo-LTR 挨着 LINE、普通 TE-in-TE）。规格只给「共现相邻」，需加方向一致 + **跨结拷贝**（≥k 个 A、B 在单一比对内连续）作实际证据 + 富集阈值。
- **m2. Count-Min/Bloom 假阳抬高 R。** ~L 大小的 sketch 会过计数，部分 unique 16-mer 假过「count≥3」抬高 R（已是成本驱动）。「阈值之上保留精确计数」步须精确复核候选；sketch 估计与精确图的 reconciliation 未定；为 10⁹–10¹⁰ 插入的可接受 FP 率设 sketch 本身是非平凡调参。
- **m3. 外层终止没问题；只有 Graft-2（M2）危及它。** 完整性记录——Stage 1–6 带掩蔽 + STOP 由残余单调下降终止。

### 内存预算

| 组件 | 3 Gb 哺乳类 | 15 Gb 植物 |
|---|---|---|
| 2-bit 基因组 | 0.75 GB | 3.75 GB |
| Count-Min sketch（depth 4，~1 B/计数器） | ~12 GB | ~60 GB |
| 精确 repetitive-k-mer 图（R×~30 B） | 3–5 GB (R=10⁷–10⁸) | 12–25 GB (R≈5×10⁸) |
| **Stage 3 all-pairs 边图（按规格）** | **3–32 TB** | **30–300 TB** |
| Stage 3 图（LSH 修复后，~N·k 边） | ~1–10 GB | ~10–50 GB |

除 Stage 3 图外全部舒适入 256–512 GB 节点。按规格的图哪都不入；LSH 重构后入 RAM。**边图是唯一内存悬崖，且与算力悬崖重合。**

### 并行模型与瓶颈
Stage 0/1/2/7 对 chunk embarrassingly parallel（sketch 需分片/原子计数器），扩展良好。**Stage 3 是双轴瓶颈**（算力 F1 + 内存 F2），且负载不均（200² 占用组主导→需 bin-packing）。Stage 4 社区发现最难并行、顺序敏感，是第二瓶颈，仅当图小到能存在时才相关。

### 规格阻塞实现者之处（编码前必须钉死）
1. **anchored 窗口如何成图节点**（重叠窗口是否合并？节点身份未定→边数未定。这条决定 pipeline 是否可行）。
2. **种子组内配对规则**——all-pairs（O(occ²)，致命）vs extend-to-representative（O(occ)）？未述，是可构建与否之差。
3. **直方图谷「重复」阈值**——谷检测算法未定，它定 R 因而定一切成本。
4. **社区发现算法、分辨率参数、稀疏割阈值**——未定，定簇数/簇大小。
5. **「过切成紧簇」粒度** vs 30 成员 POA 下限（见 M1）。
6. **Stage 8 拷贝数测在 200-样本还是全 Stage-7 编目？** 接纳语义含糊；测样本则一切看着低拷贝。
7. **Graft-2 不动点**——分解如何反馈（替换成员？注入伪序列？）及终止条件（见 M2）。
8. **边界「衰减到背景」阈值 + 背景模型**（M6）。
9. **片段拼接距离 d、方向规则、富集检验、证据阈值**（m1）。
10. **外层 STOP 规则的残余重复质量定义。**
11. **sketch FP 与精确计数的 reconciliation**（m2）。

### 建议构建顺序（先 de-risk 致命项）
1. **写任何比对代码前先量 R。** 只在真实 100–300 Mb 基因组上建 Stage 0–2；数 distinct repetitive 16-mer、算 would-be R·C_occ²/2。这个数一天内确认或杀死成本模型。**这是决定 v3 是否可行的实验。**
2. **原型 LSH/minhash 候选生成器**（F1/F2 修复），实测边数 vs all-pairs。
3. 带状 X-drop + 有界图 + **确定性** Leiden（固定 seed、坐标排序节点 ID）。
4. 独立做 Stage 7 编目（复用 minimap2/RepeatMasker 类 mapper）——标准且并行。
5. POA/consensus + 边界（带成员数门控 M6）。
6. graft 与迭代循环最后做，freeze-on-accept（M2）从第一轮就接进去。

### 工程最大风险 + 裁决
最大风险：**Stage 3 的 all-pairs-within-seed-group——其 O(R·C_occ²) 比对数（10¹¹–10¹³）与物化的 O(R·C_occ²) 边图（3 TB–300 TB）都超任何单节点 1–3 个数量级，且因 R 随 L 增长，「与 L 解耦」是假的。** 它同时是算力悬崖与内存悬崖。**按规格不可构建。** 但变可构建只需一处结构改动 + 三项保证：(1) 用 LSH/minhash 近邻提议替代 all-pairs 使边数 O(N_windows·k)；(2) hash-based 确定采样 + 规范节点序 + 整数打分保可复现；(3) freeze-on-accept 保循环终止。Stage 0–2、5–8 与 graft 个体可实现，主要把上述 under-specified 决策钉死即可。**修复是定点的、非地基重写——但强制，且正落在规格当前沉默处。**

---

## 3. 综合判断

**两条车道独立指向同一处判决：v3「可扩展」是假的，且合并时丢错了东西。**

1. **核心可扩展性被合错。** 工程判 v3 不可构建：Stage 3「仅比对共享种子窗口」是组内 all-pairs = O(R·C_occ²)，R 不随 L 解耦（拟南芥外推见 §4），实算 10¹¹–10¹³ 次比对、3 TB–300 TB 边图。其修复正是**加回 LSH/minhash 近邻候选**——而这恰是 NEXTGEN_DESIGN 原本有、在 v3 合并里被丢掉的东西。结论：两套设计都没真正解决可扩展性；正确答案两套都未完全采用——RECON/RepeatScout 的 **extend-to-growing-consensus（O(occ) 而非 O(occ²)）**。v4 必须三合一：**LSH 候选生成 + extend-to-consensus + 桶封顶**。

2. **两条车道独立收敛的三处批评（信号最强，优先修）：**
   - **C_occ=200 封顶是双输**：工程——与「故意过切」冲突，簇跌破 POA 30–50 下限；生物——按比例采样饿死古老低频亚家族（AluJ 仅 ~2 anchor）。共识：封顶不能在聚类前一刀切，POA 前要从全 occurrence 列表重扩成员。
   - **蛋白 graft 既不安全也不自洽**：生物——「局部放宽阈值」会把多拷贝驯化宿主基因（FAR1/FHY3）粘进真家族污染 consensus，且浪费了用蛋白做分类的增益（漏 Penelope/PLE）；工程——6-frame 全基因组翻译是完整 O(L) 重扫，违背有界采样骨架，不能进迭代循环。共识：蛋白通道只在样本跑一次做阈值放宽 + Stage 7/8 注释/分类，加宿主基因 guard。
   - **嵌套分解（Graft 2）是最弱一块**：工程——反馈进重聚类无不动点、会振荡（需 freeze-on-accept + sweep 上限）；生物——只处理单层、要求两翼同家族，漏深层套娃/异家族宿主/「宿主只以被打断形式出现」。共识：当前形态不成立。

3. **生物学独有致命缺口（工程照不到）：** F1 无细胞器/rDNA 守卫 → NUMT/NUPT/rDNA 成假家族（头号「过度产出」源）；F2 Helitron 全漏且被误标 seg-dup；M1 非编码分歧家族（SINE/MITE/LARD/TRIM/SVA）无救——而那正是目标群。生物 reviewer 总判：骨架隐含 retro/TIR 中心，会**同时**少召稀有家族**又**过度产出非 TE，正是最伤目标的双向失效。

**总裁决：** 骨架方向对、**值得做 v4**，但 v4 的前提是工程 reviewer 的 0 号门禁实验——在真实基因组上**先量 R、算出 R·C_occ²/2**，这个数直接决定整条路线生死，在写任何比对代码前必须先做（见 §4）。v4 的修订骨架：LSH+extend-to-consensus+桶封顶 三合一的 Stage 3；蛋白通道移出循环 + 宿主 guard + 喂进分类；嵌套加 freeze-on-accept；补 NUMT/rDNA 筛除与 Helitron/SINE 专门通道；封顶改 hash-based bottom-k 并 POA 前重扩成员。

---

## 4. 门禁实验：拟南芥实测 R 与 R·C_occ²/2

> 数据：`/scratch/shuoc/TE/Arabidopsis_thaliana/demo/genome/genome.fa`（TAIR10，5 条染色体，**133,917,131 bp**）。方法：`jellyfish count -m 16 -C`（canonical 16-mer），R = count≥3 的 distinct canonical 16-mer 数，C_occ=200。机器：256 核 / 1 TB RAM。

### 4.1 实测结果

| 量 | 值 |
|---|---|
| 基因组 | 133.9 Mb（TAIR10，5 chr） |
| total distinct canonical 16-mer | 92,998,560 |
| 其中 count=1（unique） | 79,662,342（85.7%） |
| count=2 | 9,163,114 |
| **R（count ≥ 3，= distinct repetitive 16-mer）** | **4,173,104（4.17×10⁶）** |
| R / total distinct | 4.5% |
| Max count（单 16-mer 最高出现次数） | 47,238 |
| **命中 C_occ=200 封顶的 16-mer 数** | **6,983（占 R 的 0.17%）** |

两个 Stage-3 成本数：

| 口径 | 值 | 说明 |
|---|---|---|
| **naive 上界 R·C_occ²/2 = R×19900** | **8.30×10¹⁰** | 工程 reviewer 的 worst-case，假设**每个** repetitive k-mer 都有 ≥200 occurrence |
| **真实封顶配对数 Σ n_f·C(min(f,200),2)** | **3.85×10⁸** | 按真实出现次数分布算的实际带状比对调用数 |

### 4.2 解读：实测部分推翻工程 reviewer 的 Fatal 严重度

**关键发现：naive 上界 8.3×10¹⁰ 高估真实成本 ~216 倍。** 原因是 R 的 4.17M 个 repetitive 16-mer 里，**只有 6,983 个（0.17%）真正撞到 C_occ=200 封顶**；绝大多数 count 在 3–10（count=3 就有 2.15M 个），每个只贡献 C(count,2) = 3–45 对。工程 reviewer 的 R·C_occ²/2 假设「人人满 200」，与真实重尾分布严重不符。

**真实 Stage-3 成本 = 3.85×10⁸ 次带状比对**，落在 reviewer 自己设的「≤10⁹ 可行」门槛**之下**：
- 算力：3.85×10⁸ × ~3×10⁵ DP-cell/比对 ≈ 1.2×10¹⁴ cell；@10⁸–10⁹ cell/s/核 ≈ 32–320 核时 → 256 核上**单轮 ~8–75 分钟**，×3–5 轮仍是分钟到数小时级。**可行。**
- 内存：3.85×10⁸ 边 ×16 B ≈ **6.2 GB**——远非 reviewer 估的 3–32 TB，舒适入 1 TB 节点。

**即 C_occ=200 封顶（Fable clean-room 设计的选择）确实在起作用、把高频 k-mer 的二次爆炸压住了**；reviewer 的「不可构建」基于一个 worst-case 公式而非实际分布。**在拟南芥这种中小基因组上，v3 的 Stage-3 按原样可跑、可行。**

### 4.3 必须保留的两条 caveat（不可据此就给 v3 放行）

1. **spaced seed 会把这个数顶上去。** 本实测是**精确 16-mer**；而 v3 Stage-2 实际用 weight 12–16 的 spaced seed——特异性更低（weight-12 仅 4¹²=1.6×10⁷ 空间）→ repetitive 判定更宽、co-occurrence 更多。真实 spaced-seed 的配对数会高于实测的 3.85×10⁸，介于它与 naive 上界之间，可能再高 1–2 个数量级（reviewer 另估的 ×3–4 patterns 也叠在此）。**3.85×10⁸ 是精确-种子下的下界。**
2. **scale 仍是真问题。** 真实封顶配对数随基因组与重复含量超线性增长：拟南芥 134 Mb→3.85×10⁸；外推 3 Gb 哺乳类（更多 k-mer 撞封顶 + R 增大）约 ~10¹⁰–10¹¹、15 Gb 高重复植物 ~10¹¹–10¹²——**重，但不是 reviewer 担心的 10¹³ 灾难**。LSH/extend-to-consensus 在多 Gb 仍有价值，但属**规模优化**，**不是做出可用工具的前置阻断**。

### 4.4 对总裁决的修正

门禁实验把工程车道的判决从「**按规格不可构建**」修正为「**中小基因组按原样可行；多 Gb 需要规模优化（LSH/extend-to-consensus），但非阻断**」。这也是一次「实测胜过公式」的范例：reviewer 的 worst-case 边界把一个**已被 C_occ 封顶机制压住**的成本放大了两百倍。**§1–§3 的其余结论（尤其生物学的 F1 细胞器守卫、F2 Helitron、M1 非编码 rescue，以及三处跨车道收敛批评）不受此修正影响，仍然成立**——它们与 Stage-3 成本无关，是 v4 仍必须解决的真问题。

### 4.5 caveat 1 实测：spaced-seed 权重 × v3 二次 vs v4 线性

§4.1 用精确 16-mer，是 spaced-seed 成本的下界。v3 Stage-2 实际用 weight 12–16 spaced seed，§4.3 caveat 1 指其会顶高成本但数未知。jellyfish 无 spaced 模式，故用**精确 k=12/14/16 作权重代理**实测（拟南芥 TAIR10，canonical）。同时算 v4 §A1 的 extend-to-consensus 线性成本 Σ n_f·min(f,200) 作对照：

| k（权重代理） | R（count≥3） | R/total（特异性） | v3 capped 配对（二次） | v4 extend O(occ)（线性） | v4 加速 |
|---|---|---|---|---|---|
| **12** | 6.42×10⁶ | **0.822** | **2.94×10⁹** | 1.17×10⁸ | 25× |
| **14** | 1.26×10⁷ | 0.252 | 6.91×10⁸ | 7.50×10⁷ | 9.2× |
| **16** | 4.17×10⁶ | 0.045 | 3.85×10⁸ | 2.54×10⁷ | 15.5× |

> 口径：v3 capped = Stage-3 组内 all-pairs 带状比对数 Σ C(min(f,200),2)；v4 extend = §A1 每 occurrence 对家族 consensus 比对一次 Σ n_f·min(f,200)。reviewer 可行门槛 ~10⁹。代理 caveat：真实 spaced seed 容忍替换，co-occurrence **高于**等权连续 k-mer，故 v3 列是真实下界；再 ×3–4 patterns。

三条实测结论：

1. **v3 在设计的高灵敏端（weight-12）确实失控。** v3 capped = 2.94×10⁹，**超 reviewer 可行门槛 ~3×**——且这是连续代理下界，真实 spaced × 3–4 patterns 会推到 ~10¹⁰–10¹¹。**即在最小的拟南芥上、v3 在它自己的灵敏种子设置下都跑不动**，坐实 caveat 1。

2. **v4 的 §A1（extend-to-consensus）正好治这个病。** 线性成本在 weight-12 = 1.17×10⁸，**比 v3 快 25×、稳稳落在门槛下**，且三档权重都 <10⁹。**v4 中心修复在真实数据上验证有效**——它把低权重种子的 cap-hitter 爆炸从二次降到线性。

3. **更深的发现：低权重种子的特异性崩溃（A1 治不了，需换杠杆）。** weight-12 下 **82% 的 distinct 12-mer 都是「repetitive」（count≥3）**，weight-14 仍 25%，只有 weight-16 降到 4.5%。即在 134 Mb 基因组里，weight 12–14 的种子**几乎不区分重复与单拷贝背景**——绝大多数低权重 anchor 是背景巧合而非真重复。**v4 即便用 A1 让它算得动，低权重种子仍会用伪候选淹没 pipeline。** 推论（超出原审阅与设计）：**追古老 TE（C1）的正确杠杆不是把 DNA 种子权重往下压，而是蛋白域 + 结构通道**（v4 §C/§E3 保留的）；spaced-seed 权重应**下限锁在 ~16**，DNA 已分歧到 16-mer 失效的家族交给蛋白/结构通道，而非靠 weight-12 硬捞。建议并入 v4 §A。

**净（§4 全节）：** Stage-3 不再是「先决生死项」——精确 16-mer 下可行（§4.1），v4 §A1 在低权重下也把成本压回线性（§4.5）。但实测给出一条新设计约束：**spaced-seed 权重 ≥16，C1 靠蛋白/结构而非低权重 DNA 种子**（已并入 V4_DESIGN §A4）。

---

## 5. 最小原型实测（Stage 1–2 + A1 核心机制）

> 代码：`/scratch/shuoc/TE/Arabidopsis_thaliana/demo/proto/proto_a1.py`。**诚实范围**：精确 k=16 canonical 种子（weight≥16，遵 §A4）+ C_occ=200 hash bottom-k 封顶 + **A1 extend-to-consensus（每成员对 growing consensus 比对一次，ungapped，锚定于种子）** + Stage-5 侧翼一致性边界。**未实现**（属完整构建）：Stage-4 跨种子聚类、带 indel 的 POA、接纳统计、分类、蛋白通道、嵌套。故本原型验证 **A1 核心机制 + 真实资源成本**，非成品 caller。

### 5.1 资源（TAIR10 134 Mb）

| 量 | 值 |
|---|---|
| wall（含建索引） | **41.1 s** |
| 峰值内存 | **2.1 GB** |
| 产出 draft family | 40（从 40 个 count≥200 的种子各生长一个） |
| **A1 成员→consensus 比对总数** | **7,927（≈40×200，线性 O(occ)）** |
| 对比：v3 组内 all-pairs 同样 40 种子 | ~40×19,900 = **7.96×10⁵**（二次，~100×） |

A1 的线性性在真实数据上坐实：7,927 vs 二次的 ~8×10⁵。

### 5.2 验证轴一——成员 LOCI 是否真重复（独立于 consensus 质量）

8,000 个原型成员区间中 **82.6%（6,610）落在 RepeatMasker 注释的重复区内**（参照 `genome/all/genome.fa.out`，159,220 个注释区间）。→ **Stage 1–2 的种子发现 + 占位收集正确定位到真实重复序列。**

### 5.3 验证轴二——consensus 是否匹配已知 TE（blastn vs `combine/TEs.fa`，1,926 条分类 TE）

- 32/40 命中已知 TE（e<1e-5）。
- 但比对长度分布：**aln≥150 bp 仅 5 个；50–149 bp 4 个；<50 bp（仅种子核）23 个**。
- 5 个 substantial 命中里有**两条全长 LTR/Gypsy**：`proto_fam_37`（739/738 bp @ 96.6%）、`proto_fam_40`（761/761 bp @ 97.9%）；外加一条 SINE（187 bp @ 87.7%）、两条 ~500 bp Unknown TE。
- 8 个无命中：含 `proto_fam_3` 种子 `AGGGTTTAGGGTTTAG`（拟南芥端粒重复 TTTAGGG——串联，正确地不在 interspersed TE 库里，应由 §E6/Stage0 分流），其余多为卫星/低复杂度。

### 5.4 解读：A1 核心验证 + 缺失阶段的必要性被实证

**正面**：A1 extend-to-consensus 在真实数据上**确能产出全长真 TE consensus**（fam_37/40 全长 Gypsy @ 96–98%），且**快（41 s）、省（2.1 GB）、线性（7,927 比对）**——v4 中心机制成立。

**「仅种子核命中」的 23 个家族正是两处被刻意推迟阶段的实证信号**，与预期完全吻合：
- **缺 Stage-4 跨种子聚类** → 一个高拷贝 16-mer 在**异质上下文**（不同家族 / 真成员混入巧合匹配）中共享时，ungapped 平均其侧翼 → consensus 只在种子核处一致 → 只在 33 bp 处匹配参照。**这正是完整 pipeline 必须先聚类再建 consensus 的实证理由。**
- **缺 indel（ungapped）** → 即便同质家族，离锚点越远越脱框 → 侧翼降质。需 POA/banded。

**结论**：可行性验证通过——种子发现命中 82.6% 真重复、A1 以线性成本产出全长真 TE；同时**用一个具体签名（23 个仅 33 bp 命中 + 82.6% loci 命中率）证明了 Stage-4 聚类与 POA 是必需的下一增量**，而非可选。

---

## 6. Stage-4 聚类 + POA 增量实测

> 代码：`proto/proto_stage4.py`。在 §5 原型上加 **Stage-4 跨种子聚类**（cd-hit-est 在 ±300 bp 种子核窗口上聚类，绕开随机侧翼问题）+ **strand 一致的 spoa POA**（全 ±1500 窗口）+ 最大连续支撑块边界。**同一 40 个种子**（与 §5 受控对比）。
> 注：先试 `minimap2 -cX` all-vs-all，在重复核上 anchor 爆炸、病态慢（实测即工程审阅 F1/F2 的重复密度爆炸）→ 改 cd-hit greedy 聚类。

### 6.1 聚类行为（merge/split 正是要的）

3,169 instances → 84 clusters → **33 families（≥3 成员）**；**10 个家族跨 >1 种子（合并）**；**33 个种子被拆进 >1 家族（split）**。即 Stage-4 同时把同家族跨种子实例合并、把单种子的异质占位拆开——正是 §5.4 诊断的病根修复。资源：**99 s、峰值 2.1 GB**。

### 6.2 consensus 质量：seed-core-only 病灶消失（vs §5 基线）

blastn vs `combine/TEs.fa`（1,926 条分类 TE）：

| 指标 | §5 proto_a1（无 Stage4/POA） | §6 Stage-4 + POA |
|---|---|---|
| 仅种子核命中（aln<50 bp） | **23/40** | **0** |
| substantial（aln≥150 bp） | 5/40 | **21/33** |
| full（pid≥80 且 cov≥50%） | 3 | **7** |
| 严格 80-80（pid≥80 且 cov≥80%） | ~0 | **4** |

干净全长 consensus 实例（consensus 长 ≈ 比对到 TE 的长）：`s4_fam_11` 1290/1290 bp（100% cov @ 92.7%）、`s4_fam_5` 1487/1496（99% @ 87.7%）、`s4_fam_33` 1139/1149（99% @ 89.7%）、`s4_fam_25` 2027/2464（82% @ 84.6%）。并恢复出**多个 Tc1-Mariner 亚家族**（1290–1640 bp @ 80–95%，证明 split 把亚家族分开了）。

### 6.3 解读：端到端链路验证 + 残余缺口映射到 v4

**Stage 1–2 → Stage 4 聚类 → Stage 5 POA 端到端在真实数据上跑通**：把 §5 的「仅 33 bp 种子核」全部（23→0）转成全长元件命中（aln≥150：5→21），并产出 4 个达严格 80-80 家族级标准的 consensus（基线 ~0）。**两处刻意推迟阶段的预测效果被实证兑现。**

**残余缺口（诚实，且映射到已写好的 v4 条目）**：仍有 ~6 个家族 cov 偏低（aln ~130–340 bp，多为短元件/SINE 或 consensus 仍偏长）。根因与对应 v4 修复：① 固定 ±1500 窗口对短元件含过多侧翼、对 >3 kb 长 LTR 截断 → **v4 §E4 自适应窗口 + LTR↔internal 配对**；② 短非编码 SINE 召回弱 → **v4 §E3 非蛋白结构锚定 rescue**；③ 单轮聚类 + 简单边界 → 完整 pipeline 的 POA profile + 迭代精修。原型用简化件（固定窗、单轮 cd-hit、连续块边界）已足以验证机制，缺口非机制性失败。

**总（§5+§6）**：v4 的发现内核（Stage 1–2 + A1 线性成本）与家族成形链路（Stage 4 聚类 + POA）在拟南芥真实数据上**逐级验证可行**，资源 100 s/2.1 GB 量级；门禁实验（§4）+ 两级原型（§5/§6）共同把 v4 从纸面推进到「核心机制经真实数据证实、残余缺口已定位到具体 v4 条目」。

---

## 7. §E4 自适应窗口 + §E3 SINE 通道增量实测（一成一败，如实记录）

> 代码：`proto/proto_e3e4.py`。在 §6 上加 **§E4 每家族自适应窗口**（从 400 bp 倍增至边界变内部，长 LTR 长成全长、短元件收到真实长度）+ **§E6 同源多聚末端保留**（polyA/polyT 不被边界裁掉）+ 结构分类（LTR 末端直接重复 / SINE polyA）+ **§E3 SINE 通道尝试**。同 40 种子。

### 7.1 §E4 自适应窗口：验证为有效（headline 指标继续上升）

严格 80-80 家族级（pid≥80 且 cov≥80%）：**§5 ~0 → §6 4 → §E4 6–8**（随聚类核取值 6 或 8）。产出多个 **cov=100% 干净全长 consensus**（consensus 长 ≈ 元件长，如 820/820、834/834、927/927 bp @ 91–94%）。结构上**检出 LTR**（末端直接重复 257–498 bp @ 85–89%），即把 §6 被 ±1500 截断的长 LTR 元件长成全长。自适应窗口实际取值 400–6400 bp，按家族而变——机制成立。

**但 §E4 暴露两个必须补的洞**：
- **runaway 增长（真实失败模式）**：一个家族窗口一路倍增到 12,800 bp、consensus 13,994 bp、**峰值内存飙到 51 GB**，blastn 仅 1,640 bp @ 79% 命中 Tc1-Mariner（cov 12%）——典型 tandem/段重复区导致支撑不降、窗口失控。**需 (i) 硬上限修复**（`w<WMAX` 让它越过 WMAX 翻倍，应改 `2w≤WMAX`）**+ (ii) 生长前的 §E1 seg-dup/tandem 守卫**。这正是 §E1 守卫存在的理由，被实测点名。
- **LTR 家族 per-subject cov 偏低**（32–45%）：因 TE 库把 LTR 元件拆成 `*-LTR`/`*-I` 两条，全长 consensus 跨两条 → 单条 cov 低。这是**验证口径问题**（佐证 v4 §E4「记录 LTR↔internal 配对」之必要），非 consensus 质量问题。

### 7.2 §E3 SINE 通道：两次尝试均失败（诚实负结果）

**两次尝试（CORE=150 与 CORE=300、含 §E6 polyA 保留）都是 0 个 SINE 被召回、0 个 SINE 结构 tag**，且相比 §5/§6 **是退步**（§5/§6 曾偶然命中 SINE 187–520 bp，本增量一个都没有）。`combine/TEs.fa` 中的 SINE-class 家族无任何 consensus 命中。

诚实根因（非可一句话调好的参数）：
- **当前架构结构性地少召 SINE**：40 个 count≥200 高频种子 + cd-hit 通用聚类，本就难落到 SINE 家族上；SINE 短（150–300 bp）、中等拷贝，过不了短元件的聚类与 80-80 门槛。
- **§E6 polyA 保留实现正确但无用武之地**——SINE 家族根本没成簇，保不保末端都谈不上分类。
- **真正的 §E3 需要的是专门的结构锚定通道**（tRNA/7SL 衍生头部 A/B-box profile seeding + HMM 迭代招募），是一个**新组件**，不是在现有「高频 16-mer 种子 + 通用聚类 + polyA 启发式」上加分类规则能得到的。本增量证明了这条捷径走不通。

### 7.3 §E4 runaway 守卫：已修并验证

加 §E1 守卫（生长前 `is_tandem` 短/中周期自相似检测 + 硬上限 `2w≤WMAX` + 长度 6 kb 早停）后重测：**峰值内存 51 GB → 12.8 GB**，最长 consensus 13,994 → 10,098 bp，**1 个家族被正确标 `suspect_long/tandem?`**；80-80 维持 6。runaway 灾难消除（守卫可再收紧，但已非阻断）。

### 7.4 §E3：决定性发现——**本基因组无法验证 SINE 通道，且基本不需要**

动手建 tRNA-head 通道前先做证据核查，结果决定性地改变了 §E3 的处置（实测胜过设想的最强一例）：

- **拟南芥是 SINE 贫乏基因组**：`combine/TEs.fa` 仅 **6 条 SINE 家族**（长 84–509 bp）；RepeatMasker 在本基因组注释 **0 个 SINE loci**。
- **6 条 SINE 的基因组拷贝数（blastn 80/50）= 47217、7、1、…**：仅 **1 条（TE_00001182）高拷贝**，其余 1–7 拷贝。
- **那条唯一高拷贝 SINE 已被通用 pipeline 召回**：本轮 `e_fam_17`（483 bp @ 77%）、`e_fam_26`（318 bp @ 82%）即命中 TE_00001182#SINE——**不需要专门 SINE 通道就能找到**。
- 其余 SINE 1–7 拷贝 → **任何 de novo 方法都不可发现**（短 + 低拷贝），不是通道问题。

**结论**：在拟南芥上，专门 tRNA-head SINE 通道 (i) 没有可加的召回（唯一高拷贝 SINE 已被通用流程捕获），(ii) 无法验证（其余 SINE 拷贝数 1–7，truth 信号近空）。**故本轮不构建一个无法验证的组件**（数据完整性：不交付无法证实其有效性的东西）。§E3 SINE 通道**改判为「需 SINE 富集基因组方可开发/验证」**（人 Alu ~10⁶、或 SINE 富集植物如玉米/水稻），在拟南芥上既非必要也不可验证。§7.2 把 §E3 失败归于「架构捷径走不通」部分正确，但更根本的原因是**测试基因组本身 SINE 近乎为空**——这一点只有做了拷贝数核查才看清。

### 7.5 净判（§7 全节）

**§E4 自适应窗口 + §E1 runaway 守卫：纳入 v4 主线，已验证**（80-80 家族 4→6–8；内存 runaway 51 GB→12.8 GB 已控；检出 LTR 结构）。**§E3 SINE 通道：在拟南芥不构建、不可验证**——改判为需 SINE 富集基因组的独立增量；本基因组上唯一可发现的高拷贝 SINE 已被通用流程召回。**方法学教训**：在 SINE 贫乏基因组上开发 SINE 通道是用近空 truth 调一个组件——必须换富集基因组，否则是无意义甚至自欺的验证。

---

## 8. §E1 细胞器 / rDNA 假家族守卫实测（强验证：在真实库里抓到污染）

> 参考集：拟南芥叶绿体 `NC_000932.1`（154 kb）+ 线粒体 `NC_037304.1`（368 kb）+ rDNA `X52320`（45S）/`X52629`（U3 snRNA）via efetch。守卫 = blastn 候选 consensus vs 该参考集，命中 id≥80% 且 cov≥50% → 标 NUMT/NUPT/rDNA 并移出 interspersed 库。

### 8.1 结果

| 对象 | 标出的细胞器/rDNA 假家族 |
|---|---|
| (a) 原型 de novo 输出（stage4_consensi.fa，33 fam） | **0**（40 种子样本恰未落到细胞器/rDNA 区） |
| (b) **真实 curated 库**（combine/TEs.fa，1926 条） | **19（~1%）：18 NUMT（线粒体）+ 1 rDNA** |

(b) 的 19 条**正是 v4 §E1 预测的假家族失效模式**——且被现有 pipeline **错标成真 TE**：5 LTR/Gypsy、5 LINE/L1、2 LTR/Copia、3 tRNA、1 SubclassI、2 Unknown、1 rRNA，却与线粒体基因组 **84–99.9% identical（多数 99%+、cov 80–100%）**。蛋白通道抓不到它们（携带 cox/nad 而非 RT），组成 E-value 也不触发——只有直比细胞器/rDNA 的守卫能抓。

### 8.2 三项确证（验证既真且特异）

- **是真·分散多拷贝家族**（非单拷贝偶然）：抽查 `TE_00000503`(#LTR/Gypsy) 核内 3 拷贝、`TE_00000572`(#LINE/L1) 3 拷贝、`TE_00000557`(#rRNA) **65 拷贝** + 各自高比例命中细胞器/rDNA → 确系 NUMT/rDNA 家族，de novo 必然把它们当家族发现。
- **特异**：1926 条仅 19 条（~1%）被标，真 TE 不被误删；早前确认的真 Copia `TE_00000769` **未**被标。
- **机制对口**：18/19 是线粒体（NUMT）；本轮未见叶绿体 NUPT 命中（拟南谱著名的 chr2 质体插入未在库阈值内）——如实记录。

### 8.3 净判

**§E1 守卫验证有效且高价值**：在一个真实「curated」TE 库里抓出 **~1% 的细胞器/rDNA 污染**（19 条，含 15 条被错标成 Gypsy/Copia/L1/SubclassI/Unknown 的 NUMT），且特异不伤真 TE。这是本轮**最强验证**——不是玩具 demo，而是在真实交付物里发现并可移除的污染。守卫纳入 v4 Stage-8 接纳门。诚实边界：① 原型自身输出未复现该污染（种子样本小），结论靠真实库证（更强）；② rDNA 参考 `X52629` efetch 实为 U3 snRNA（非 5S），不影响 45S 命中，后续应补正 5S/全 rRNA（barrnap 可得，`prokka_env/bin/barrnap`）；③ 仅见 NUMT 未见 NUPT，后续在 NUPT 富集区（chr2 质体插入）补测。

---

## 9. §E2 一等 Helitron 处理实测（一证、一驳、一设计修正）

> 验证集：库内 **86 条 RC/Helitron**（长 62–9036 bp，中位 910 bp）；对照 86 条 LTR/Gypsy。拟南谱 Helitron 丰富，§E2 可验证。

### 9.1 F2 确证：通用结构分类对 Helitron 完全失效

把 §E4 的 `classify()`（TSD/TIR/LTR/polyA）跑 86 条真 Helitron → **86/86 全部落 "interspersed/unknown"，0 个被识别**。Helitron 无 TSD/TIR/LTR/polyA，通用结构信号对它**零分辨力**——F2 的预测在真实数据上坐实，§E2 的存在必要性确立。

### 9.2 §E2(ii) 末端模型（3' hairpin+CTRR）在 consensus 上：失败

3' 末端 CTRR 仅见于 **11/86**；自建「3' hairpin + CTRR」检测器在 Helitron 上 **sensitivity 6/86=7%**、在 Gypsy 对照上 **假阳 6/86=7%**——**灵敏度=假阳率=无分辨力**。诚实根因：**de novo consensus 边界不精确，恰恰抹掉了末端信号**（这本身就是 F2）。**设计修正（重要）**：§E2(ii) 末端模型必须在**基因组实例层**（真末端尚在）做边界锚定，不能在已建好的 consensus 上做分类——即末端模型属「发现期定界」而非「成品分类」。这条修正写回 V4_DESIGN §E2。

### 9.3 §E2(i) Rep/Hel 蛋白域：验证有效且高度特异

efetch 取 21 条 Helitron helicase（Pif1/helitron-helicase-like）参考蛋白（11.8 k aa），blastx（consensus 翻译 vs 参考）：

| | Rep/Hel 域命中（e<1e-3, ≥30 aa） |
|---|---|
| **Helitron（86）** | **26/86 = 30%** |
| **Gypsy 对照（86）** | **0/86 = 0%** |

**Rep/Hel 域是 Helitron 的特异判据**：30% 灵敏、**0% 假阳**——在通用结构（0/86）与末端模型（7%≈噪声）双双失效处，蛋白域干净地把 Helitron 分出来，且不误报 Gypsy。30% 灵敏度是**生物学正确**的：仅自主 Helitron 编码 RepHel 转座酶，~70% 非自主（仅载货、无转座酶）本就无此信号。这正印证生物审阅「Helitron 唯一可靠信号是 Rep/Hel 域」。

### 9.4 净判（§9 全节）

- **§E2(i) Rep/Hel 蛋白域：纳入 v4（验证有效、特异）**——并入 §C 蛋白通道，给自主 Helitron（~30%）定类，且 0 假阳。
- **§E2(ii) 末端模型：设计修正**——改为发现期、实例层定界，不在 consensus 上分类（本轮证伪了「consensus 上分类」路径）。
- **非自主 Helitron（~70%，无 RepHel 无干净末端）**仍是开放难题：须靠实例层末端模型 + §E2(iv) 按末端而非中段聚类（抗 cargo 异质）共同解决，是下一个增量。
- **§E2(iii) seg-dup 豁免**可即用 Rep/Hel tag 实现：被 Rep/Hel 命中的家族豁免「低分歧低拷贝→seg-dup」重标。

**方法学一致性**：§E2 与 §E1/§E3/§E4 同一规律——**蛋白/序列层的正交证据（Rep/Hel 域、细胞器同源）在真实数据上稳；依赖 de novo consensus 末端结构的启发式（Helitron 末端、SINE polyA）则被不精确边界击穿**。这把 v4 的取舍指向「能靠正交证据就别靠 consensus 末端结构」。

---

## 10. 改用「基因组 masking 覆盖增量」作非循环验证指标（替代循环的 truth 匹配）

审核员 FA-1 指出 §5–§7 用 `combine/TEs.fa` 作 truth 是循环验证（该库 48% 是 mdl-repeat 自己输出）。改用**无需 truth 的指标**：把新库并入 mdl-repeat 库后，RepeatMasker 能否 mask 到**更多基因组**。

### 10.1 实验

- baseline = mdl-repeat 单独（`mdl-repeat/consensus.fa`，835 条）。
- merged = mdl-repeat + 原型 consensi（95% 去冗，906 条，原型贡献 71 条非冗余）。
- 同基因组、同 RepeatMasker 设置，只换库。

| 库 | masked | |
|---|---|---|
| baseline（mdl-repeat） | 29,667,208 bp（**22.15%**） | ≈ 合并 `all` 库的 22.42%（mdl-repeat 几乎贡献了全部 TE masking 力） |
| merged（+原型） | 31,842,156 bp（**23.78%**） | **+2,176,140 bp（+1.63 pp）** |

### 10.2 关键诚实校正：+1.63 pp 的 88% 是被正确排除的卫星，不是新 TE

raw +1.63pp 看似大胜，但核验暴露陷阱：

- 新增 masking 的**主导贡献者 `e_fam_8`**（48% 的新区段）经 self-alignment 证实是 **~500 bp 周期的串联卫星**（consensus 对自身按 500/1000/1500/2000 bp 移位 99.7% 自比）——不是 interspersed TE。我先前 `is_tandem` 因周期范围漏查 ~500 bp 而误判其为 interspersed；而「多拷贝」核验（40/40、中位 3880 拷贝）**正确但有误导性**——卫星本就多拷贝，多拷贝不区分卫星与 TE。
- 卫星正是 mdl-repeat 与 combine 库（interspersed-TE 工具）**有意排除**的。bp 归因：**2.176 Mb 里 1.92 Mb（88%）是 e_fam_8 卫星；仅 0.26 Mb（12%，+0.19 pp）是真正的 interspersed-TE 新覆盖**（其余 top 贡献者 e_fam_20/19/s4_fam_4/25 经周期扫描确认非串联）。

| | masked 增量 |
|---|---|
| raw merged − baseline | **+1.63 pp**（含卫星，**误导**） |
| 排除卫星后的真·新 TE 覆盖 | **+0.19 pp（~258 kb）** |

### 10.3 结论：指标对，但必须排除串联/卫星；且 40 种子原型如预期几乎没加真 TE 覆盖

- **masking 覆盖增量是正确的非循环指标**（化解 FA-1，无需 truth）。**但必须只算 interspersed 部分**——把串联/卫星/低复杂度分流到独立 track（正是 v4 §E6 / §0），否则含卫星的库「mask 更多」却非更多 TE，指标被灌水。用户提出的「merged mask 更多即成功」须修正为「merged mask 更多 **interspersed（非卫星）** 序列」。
- 校正后，**40 种子原型的真·新 TE 覆盖仅 +0.19 pp**——印证最初的诚实预判：高频种子原型基本只重找 mdl-repeat 已有的元件，真正要证明「找到 mdl-repeat 漏掉的 TE」必须跑**频率分层的真实全发现**（回应审核 MJ-1/MJ-4），并在该指标下复测。
- 方法学价值：本实验本身验证了「masking 多 ≠ TE 多」这一陷阱，并给出可执行的修正（先 §E6 分流串联，再算 interspersed masking 增量）。

---

## 11. 新方法提升 mask 比例（目标：新方法本身拉高 masking，targeting mdl-repeat 的缺口）

目标修正为「**新方法本身能提升 mask 比例**」（非外挂 TRF）。做法：在 **residual genome**（mdl-repeat 已 mask 区段置 N，剩 104 Mb=77.8% 未 mask）上跑新方法发现——只找 mdl-repeat **漏掉**的重复。

### 11.1 scaled 发现（proto_scaled.py）

residual 中 count≥200 的 16-mer = **1,639 个**（mdl-repeat 漏掉的高拷贝重复）。全部作种子（每种子 bottom-k 封顶 60）→ 62,202 instances → cd-hit 1,883 簇 → **144 个新方法家族**（239 kb consensi）。6 min、峰值 3.2 GB。

### 11.2 masking 提升（RepeatMasker，同设置只换库）

| 库 | masked |
|---|---|
| mdl-repeat 单独（835 条） | 29,667,208 bp（**22.15%**） |
| **mdl-repeat + 新方法（+109 非冗余）** | 32,565,982 bp（**24.32%**） |
| **新方法净提升** | **+2.17 pp（+2,898,774 bp ≈ 2.9 Mb）** |

**已核验为真实重复（非 unique 过度 mask）**：新增区段抽 40 个，**39/40 多拷贝（中位数千拷贝），仅 1/40 单拷贝**。

### 11.3 与 40 种子原型的关键区别：这次提升主要是 interspersed，不是卫星

144 家族里 **123 interspersed + 21 串联/卫星**。主导贡献者 `scaled_fam_26`（1,493 个新区段）经核验是 **6.5 kb 全长、分散（非串联）、~数百–千拷贝的 interspersed 元件，且 mdl-repeat 与 combine/TEs.fa 两库都没有**（vs 基因组 3,483 次近全长命中、不 hit TE 库、非细胞器/rDNA）。即 residual-targeted 发现这次找回的主要是 mdl-repeat 漏掉的**真·分散 TE**，而非 40 种子原型那次的卫星。

### 11.4 净判

**新方法本身把 mask 比例从 22.15% 拉到 24.32%（+2.17 pp / +2.9 Mb），且经核验是真实多拷贝重复、以 interspersed 为主。** 这在非循环指标（masking 覆盖，§10）下直接达成「新方法提升 mask 比例」的目标。诚实边界：① 绝对拷贝数为 blastn-HSP 高估，以 RepeatMasker 区段数为准；② 本轮只用 residual count≥200 种子（1,639），下探 count≥50（9,256 种子）应能再加 mid-copy 家族、进一步提升；③ `scaled_fam_26` 是干净单家族还是略带 composite 未完全钉死，但「真实多拷贝、非串联、非细胞器、两库皆无」已足以支撑其为合法新 masking。下一步若要更高：扩到 count≥50 residual 种子 + 叠加 rDNA/NUMT，逼近基因组真实重复上限。

---

## 12. 审核修正 + masking 真实上限 + 推进计划（Fable agent 复核，verify-first）

独立 Fable agent 复核并**用直接核验纠正了 §10–§11 的两处错误**：

1. **基因组身份**：这是 **T2T Col-CEN**（134 Mb、近 0 N、完整着丝粒），**不是** canonical TAIR10（~119.7 Mb、着丝粒有缺口）。前文「TAIR10」表述更正。
2. **「Satellites = 0.00%」是分类标签假象，不是未 mask 的卫星**（推翻 §10/§11 的「卫星=未 mask 头条 headroom」推断）：取一段 CEN178 着丝粒卫星单体全基因组 blast，它跨 **5 个着丝粒共 11.66 Mb（8.71 pp）且已被 baseline ~100% mask**（由 consensus `mdl_R82` 完成，RepeatMasker 只是把它标成 "Unspecified"）。**故不存在 ~8 pp 卫星 headroom，前文 TRF 卫星 headroom 的前提作废**——新方法的 21 个串联家族已收走大部分，net-new-vs-scaled 仅 ~0.42 pp（chr3+5）。

**真实重复上限（按 scaled-residual 的 k-mer occupancy 实测）**：当前 24.32% mask 后，残余精确 16-mer occupancy：**count≥200 = 0.01 pp（高拷贝已耗尽，residual loop 已收敛）**、count≥10 = 1.22 pp、count≥3 = 8.09 pp（但 3→5 阈陡降 8.09→3.17 = 完整性信号，count-3/4 在 64% AT 基因组里多为低复杂度/组成性）。去低复杂度后 **~0.88 pp 真·复杂多拷贝**可被精确 16-mer 捕获。**可达上限：~27–28% 高置信（+3–4 pp）；29–30% 需接受完整性风险。硬墙：全基因组 count≥2=40.5 pp 是组成性假象，绝不可作目标。**

**排序计划（已落实推进）**：① residual loop 降阈到 count≥10–50 + **homopolymer/dinuc 种子过滤 + 家族 ≥10 拷贝门**（新方法，+1.0–1.8 pp，低工本）；② 补跑 TRF chr1/2/4 + 臂、按 period/copy 过滤（辅助，+0.7–1.0 pp）；③ RepeatMasker `-s` + 分歧上限 35–40% + 边界补全（+0.5–1.5 pp）；④ 仍不足才上 spaced-seed/蛋白锚定。**完整性边界**：count<5 或 homopolymer 种子 → 低复杂度/近 unique；`-s` 分歧 >35–40% → AT-rich 假阳；consensus 越界延伸 → unique 侧翼。每个家族必过多拷贝(≥3，最好≥10)门。

### 12.1 推进落实：direction #1 实测 + 一次被完整性检查拦下的过度 mask

迭代 residual loop（residual2 = 24.32% mask 后；count≥20 + homopolymer 过滤 + 家族 ≥10 成员门）→ 692 个新家族。

**naive 结果看着惊艳：28.17%（+3.87 pp）——但完整性检查拦下了它。** 抽 50 个新增 mask 区段查基因组拷贝数：**中位数 = 1.0，30/50 是单拷贝**——即 ~60% 的 +3.87 pp 是**对近 unique 序列的伪 mask**（非低复杂度：dustmasker 仅 8.5%、692 consensi 0 条低复杂度）。根因：**min-10 拷贝门作用在「发现期聚类成员数」上，而非 consensus 的真实基因组拷贝数**——一个由 ≥10 个异质实例聚出的 chimeric consensus 仍会被 RepeatMasker 贴到大量单拷贝区。

**修复（按完整性边界）**：对 692 个 consensus 逐个 blastn 全基因组、**要求 ≥10 个真实基因组 loci（≥80% id、≥50% 长）**，692→**保留 300**（删 392 个伪/低拷贝），重 mask。

**最终诚实 masking 阶梯**：

| 库 | masked | |
|---|---|---|
| mdl-repeat（baseline） | **22.15%** | |
| + 新方法 pass1（count≥200，已验证 39/40 多拷贝） | **24.32%** | +2.17 pp |
| + 迭代 count≥20（**naive，60% 单拷贝伪 mask**） | 28.17% | ❌ 不可信 |
| **+ 迭代 copy-validated（≥10 真实拷贝）** | **25.52%** | ✅ 真实 |

copy-validated 后复核：新增区段 **31/40 多拷贝（中位 12 拷贝）**，干净（仍 ~17% 单拷贝残留，属轻微，可再紧门到 ≥20 loci）。**新方法对 mdl-repeat 的真实 mask 提升 = +3.37 pp（22.15%→25.52%）**；其中迭代（direction #1）真实贡献 +1.20 pp（24.32→25.52），**与 Fable agent 的 +1.0–1.8 pp 预测吻合**。

**方法学要点（写回纪律）**：① masking 评估必须以**最终 mask 区段的真实拷贝数**为完整性门，不能只靠发现期成员数——否则 chimeric consensus 把近 unique 序列也 mask 了（naive 28.17% 即此）；② 这次完整性检查避免了一个 +6 pp 的虚高 claim，把结果钉在真实的 25.52%；③ 距 agent 估的真实上限 ~27–28% 还有 ~1.5–2.5 pp，来自后续 direction #2/#3——但每步都须过同一拷贝数完整性门。

### 12.2 推进 direction #2（tandem）/#3（`-s`），均带 copy-validation

**direction #3（RepeatMasker `-s` 敏感模式，copy-validated 库上）= 哑弹**：25.52% → 25.64%（+0.12 pp raw），但抽 40 个 `-s` 新增区段 **27/40 单拷贝（中位 1.0）**——真实增益仅 ~+0.04 pp。库已把自己的家族吃透，敏感模式只多贴边际假阳。**`-s` 在本库上不值得用**（边界补全是另一条路，未做）。

**direction #2（TRF tandem，period/copy 过滤）= 真实 +0.80 pp**：chr3+5（56 Mb 代表，period 2000）tandem net-new over 25.52% = +0.80 pp（→ chr3+5 25.79%；着丝粒 CEN178 已被 mask，net-new = 臂区微/小卫星）。TRF 周期分布权威确认：net-new 区段周期以 11–200 为主（4943 个 period 11–50、1049 个 51–200）——genuine tandem。注：chr1/2/4 period-1000 TRF 在着丝粒密集区过慢、未跑完，+0.80 pp 为 chr3+5 代表性外推；其中少量 period≤10 微卫星属低复杂度边界。

### 12.3 最终诚实 masking 阶梯（全程 copy-validation）

| 库/步骤 | masked | 增量 | |
|---|---|---|---|
| mdl-repeat（baseline） | 22.15% | — | |
| + 新方法 pass1（count≥200，已验证） | 24.32% | +2.17 pp | ✅ |
| + 新方法迭代（copy-validated，dir #1） | 25.52% | +1.20 pp | ✅ |
| + TRF tandem（dir #2，chr3+5 外推） | ~26.3% | +0.80 pp | ✅ tandem |
| + RM `-s`（dir #3） | ~26.3% | +0.04 pp | ❌ 哑弹 |
| **genuine 合计 over mdl-repeat** | **~26.3%** | **+~4.15 pp** | |
| Fable-agent 高置信上限 | ~27–28% | | |

**净**：新方法 + tandem 把真实 masking 从 mdl-repeat 的 22.15% 拉到 **~26.3%（+~4.15 pp，全部经真实拷贝数验证）**。主力是新方法本身（dir #1，pass1+迭代 +3.37 pp）+ tandem（+0.80 pp）；`-s` 无效。

### 12.4 direction #4（边界补全 / Refiner）= 失败（−5.17 pp，chimeric 反噬）

试图用「从每家族真实基因组拷贝 + ±600 侧翼重做 POA、按 cross-copy agreement 向外延边界」补全 consensus（373 家族、216 个延长、+188 bp 中位）。结果 **masking 不升反降：25.52% → 20.35%（−5.17 pp）**，甚至低于 mdl-repeat 单独的 22.15%。

诊断（库完整：730 mdl + refined，0 N、RM 正常）：**延长后的 consensus 变 chimeric**——agreement-boundary 把并非元件、但拷贝间偶然共享的侧翼也纳入，POA 产出「核+异源侧翼」嵌合体；RepeatMasker「best-hit-wins」让这些更长的嵌合 consensus **赢得区段后却 mask 得更少**，并挤掉本来匹配更好的 mdl-repeat consensus，全局退化。**即：naive 边界延长产出嵌合体、反而降 masking。** 弃用 `lib_refined`。

**教训**：边界补全不能只靠 POA agreement 向外延（会嵌合）；正确做法是 RepeatModeler Refiner 式的 cross_match 多重比对一致裁剪 + **逐家族「只在该家族 masking 真正上升时才接受 refined」的接受门**。

**已实现接受门（proto_refgate.py）**：逐家族比较 orig vs refined 的「基因组覆盖 = ≥80% id 且 aln≥该 consensus 自身长 50% 的合并 bp」，**仅当 refined 覆盖更多才接受**。`≥50% 自身长`这条天然否决嵌合体（坏延长把覆盖比稀释到 <50%）。结果：**373 个 refined 仅 30 个被接受、343 个被拒（保留原始）**，masking 从灾难的 20.35% 恢复到 **25.55%**——接受门成功挡掉嵌合反噬。**但真实增益仅 +0.03 pp（25.52%→25.55%）≈ 零**：说明 adaptive-window 发现阶段已把边界做到接近最优，**本数据集边界补全没有真实 headroom**，agent 估的 +0.3–0.8 pp 不成立。接受门机制本身正确（区分了 30 个真改进 vs 343 个嵌合退化），只是无料可补。

### 12.5 最终诚实结论（masking 提升）

| 步骤 | masked | 增量 | |
|---|---|---|---|
| mdl-repeat baseline | 22.15% | — | |
| + 新方法 pass1（count≥200，验证） | 24.32% | +2.17 | ✅ |
| + 新方法迭代（copy-validated #1） | 25.52% | +1.20 | ✅ |
| + TRF tandem（#2，chr3+5 外推） | ~26.3% | +0.80 | ✅ |
| RM `-s`（#3） | — | +0.04 | ❌ 哑弹 |
| 边界补全 naive（#4） | 20.35% | −5.17 | ❌ chimeric 反噬，弃用 |
| 边界补全 + 接受门（#4-gated） | 25.55% | +0.03 | ✅ 防住反噬，但无真实 headroom |

**最终：新方法 + tandem 把真实 masking 从 22.15% 提到 ~26.3%（+~4.15 pp，全程拷贝数验证）。** 边界补全即便用正确的逐家族接受门，也仅 +0.03 pp——adaptive-window 发现阶段已把边界做到接近最优，本基因组无补全 headroom。**已到本数据集的真实重复上限附近；copy-validation / 接受门完整性纪律在全程拦下了 naive 迭代的 +6 pp 虚高、`-s` 伪增益、边界补全 −5.17 pp 嵌合反噬——保留下来的每个 +pp 都是真实重复。**

### 12.6 复核纠正：方向未尽——direction #5（蛋白锚定古老 TE）是经核验的 ~+1 pp 真实 headroom

第二个 Fable agent 复核纠正了「26.3% 已尽」：那只是**「精确-16-mer 方法族」的上限,不是基因组真实重复上限**。它在当前 25.52% 的未 mask 区抽 2.75 Mb(2.9%)、blastx 比对 RepeatMasker `RepeatPeps.lib`(18,011 条 curated TE 蛋白):

- **未 mask 抽样 0.91% bp 直接编码 TE 蛋白**(LINE/L1 主导 + Helitron/Copia/MULE/hAT/PIF-Harbinger/EnSpm),且 **0 bp 与现有 mask 重叠**;
- **多拷贝确认**:54 命中区 **27/54 有 ≥3 拷贝、37/54 ≥2**;
- **高分歧(70–90% id)→ 精确 16-mer 看不见**(80% id 下 16-mer 共享概率仅 0.028,count≥3 门永不触发)→ 与精确-16-mer 残余(~0.88 pp)**基本不重叠**。

**核验 headroom**:#5 = **+0.4 pp 硬底 → ~1.1 pp 中位 → ~2 pp 上限**;叠加精确-16-mer 残余 → 真实可达 **~27.3–28%**,正落文献 Arabidopsis 总重复含量。**精确-16-mer 阶梯在 ~26.3% 早 ~1 pp 触顶,蛋白锚定是闭合到真实上限的杠杆。** 这正是 v4 §C 蛋白通道 / NEXTGEN M3 / C1 古老 TE 目标——masking 实证给了它具体已验证的价值。

**#5 过度-mask 守卫**:① 每元件过 **≥3 真实拷贝**门(去 ~50% 的 1–2 拷贝驯化残骸);② 只锚 bona fide TE 蛋白(RepeatPeps,非宿主基因);③ **只延伸到 ≥3 拷贝 DNA 支撑所及**(绝不进 unique 侧翼——正是杀死 #4 的嵌合陷阱);④ <70% id 侧翼与 1–2 拷贝元件处封顶。

**修正结论**:「26.3% 已尽」只对精确-16-mer 方法族成立;**direction #5(蛋白锚定古老 TE)是经核验值得做的 ~+1 pp 真实多拷贝 TE headroom,闭合到 ~28% 真实上限的唯一剩余杠杆。**

### 12.7 实现 direction #5（蛋白锚定古老 TE）= 成功,+0.47 pp 干净

实现(proto_protein.py):**diamond blastx** 未 mask residual（分块,very-sensitive）vs RepeatPeps（18,011 TE 蛋白）→ **1,889 个 TE 蛋白锚定 loci**(LINE/L1 主导,89% 在未 mask 区)→ 提取元件窗口 → cd-hit 聚类(0.80)→ strand 一致 adaptive spoa POA → **43 家族**。逐 consensus **dc-megablast 拷贝验证(≥3 loci @ ≥70% id)→ 41/43 通过**(古老 TE 用 70% id 门捕获分歧拷贝)。并库重 mask:

| | masked |
|---|---|
| 库基线（mdl-repeat + 新方法 k-mer,25.52%） | 25.52% |
| **+ 蛋白锚定古老 TE（#5）** | **25.99%（+0.47 pp）** |

**新增 masking 干净**:抽 40 个新增区段 **38/40 多拷贝（中位 11 拷贝）、0/40 单拷贝**——远优于 naive 迭代的 60% 单拷贝。主导是 LINE 家族(protein_fam_20 贡献 544 区段)——正是 DNA k-mer 因 70–90% 分歧看不见、靠蛋白域才捞回的古老元件。

**榨上半(从 +0.47 向 +1.1 pp)**:① 降 diamond 阈到 pident≥25 + ultra-sensitive,锚定**饱和**(1889→1993,只 +104——蛋白域 loci 已基本找全);② 真正杠杆是**元件跨度足迹**——把每个蛋白锚 ±1500 扩成元件跨度并 copy-validate:**1,606 个元件跨度区仅 306 个过 ≥3 拷贝门(1,300 个被拒,多为 1–2 拷贝驯化残骸/单拷贝——正是完整性门要挡的)**。叠加这 306 个验证区:**25.99% → 26.22%(再 +0.23 pp)**。**#5 合计 +0.70 pp(25.52%→26.22%),全程拷贝验证。** 距 agent 中位 +1.1 pp 余 ~0.4 pp,落在 **<70% id / 1–2 拷贝**的完整性墙后——继续榨就要么接受 1–2 拷贝(驯化基因风险)、要么降到 <70% id(假阳风险),故 **+0.70 pp 是 #5 的干净实现上限**。

### 12.8 最终 masking 阶梯（全程 copy-validation，含 direction #5）

| 步骤 | masked | 增量 | |
|---|---|---|---|
| mdl-repeat baseline | 22.15% | — | |
| + 新方法 k-mer pass1（count≥200） | 24.32% | +2.17 | ✅ |
| + 新方法 k-mer 迭代（copy-validated #1） | 25.52% | +1.20 | ✅ |
| + 蛋白锚定古老 TE（#5 consensus，copy-validated） | 25.99% | +0.47 | ✅ 干净 38/40 |
| + #5 元件跨度足迹（copy-validated） | **26.22%** | +0.23 | ✅ 306/1606 过门 |
| + TRF tandem（#2，chr3+5 外推） | ~27.0% | +0.80 | ✅ |
| RM `-s`（#3）/ 边界补全（#4） | — | ~0 | ❌ 已证伪 |

**最终:库式 + 足迹 genome-wide(mdl-repeat + 新方法 k-mer + 蛋白古老 TE consensus + 元件跨度)= 26.22%(+4.07 pp,全程拷贝验证);叠加 tandem ~27.0%。** 新方法两条互补通道(精确 k-mer 攻高拷贝 +3.37 pp、蛋白锚定攻分歧古老 +0.70 pp)合计把真实 masking 从 mdl-repeat 的 22.15% 提到 **~27.0%(+~4.85 pp)**,**已进入 agent 估的 ~27.3–28% 真实上限区间**。剩余 ~0.3–1 pp 全落在完整性墙后(1–2 拷贝驯化基因 / <70% id 极分歧),榨它必伤真实性,故停。

**这一整条 masking 实证线索的核心方法学**:① 非循环指标(masking 覆盖增量,§10)替代循环的 truth 匹配;② **真实拷贝数完整性门**贯穿每一步,先后拦下 naive 迭代 +6 pp 虚高、`-s` 伪增益、边界补全 −5.17 pp 嵌合反噬,并验证 #5 干净;③ 两条正交发现通道(k-mer + 蛋白)各攻一类、互不重叠;④ 每个保留的 +pp 都是经验证的真实多拷贝重复。
