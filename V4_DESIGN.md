# de novo 重复检测 v4 修订提案

本提案是对 [v3 合并方案](DESIGN_REVIEW.md#0-被审对象v3-合并方案) 的修订，逐项关闭两份独立审阅（生物 + 工程）与门禁实验暴露的问题。**只写「相对 v3 的改动」**，未提及处沿用 v3。每条改动标注它关闭的审阅项（F/M 编号见 DESIGN_REVIEW.md §1/§2）。

设计仍是纸面设计，所有性能表述为预期，需实测验证。

---

## A. 可扩展性内核重构（关闭 工程 F1/F2，吸收门禁实验 §4 修正）

门禁实测把工程「不可构建」判决修正为「中小基因组可行、多 Gb 需规模优化」：拟南芥真实封顶配对数 3.85×10⁸（naive 上界高估 216×）。据此 v4 不做「LSH 必须替换 all-pairs」的一刀切，而是**分两层、且引入两套设计都缺的正确内核**：

**A1. 主机制改为 extend-to-growing-consensus（O(occ)，取代组内 all-pairs O(occ²)）。** 这是 RECON/RepeatScout 的核心、两套设计都没完全采用的那一招：不在种子组内两两比对，而是对每个新生家族维护一条**增长中的 consensus**，新成员只与 consensus 比对一次（O(occ)），而非与所有已有成员比对（O(occ²)）。即便 C_occ=200 全满，单家族也从 19,900 次降到 ~200 次比对。这一条直接把 §4 的 caveat 2（多 Gb 外推 10¹¹–10¹²）压回线性级。

**A2. 候选生成用 LSH/minhash（仅多 Gb 启用，作为 A1 的邻居提议）。** 对每窗口的 repetitive-k-mer 集做 bottom-k minhash 签名，每窗口只提 k≈8–32 个候选邻居 → 边数 O(N_windows·k)。**桶大小封顶 + 高频 k-mer 频率掩蔽**防 repeat 富集基因组的桶爆。中小基因组（实测可行）默认走 A1 + C_occ 封顶即可，不必 LSH。

**A3. 分层触发（自适应、无人工切换）**：
- ≤ ~500 Mb：C_occ 封顶 + extend-to-consensus，全量，无 LSH（拟南芥实测此路可行）。
- 0.5–3 Gb：加 LSH 候选生成 + 频率掩蔽。
- \> 3 Gb：LSH + 分片 + 块间 consensus 归并。

> 净：v4 同时拿到 NEXTGEN 的对（LSH）、Fable 的对（C_occ 封顶）、和**两套都缺的对**（extend-to-consensus）。

**A4. spaced-seed 权重下限锁 ≥16；C1 古老 TE 改由蛋白/结构通道捕获，不靠压低 DNA 种子权重（门禁实验 DESIGN_REVIEW.md §4.5 直接驱动）。**

实测（TAIR10）暴露一条 v3、两份审阅都没点破的硬约束：**低权重 spaced seed 在中小基因组里特异性崩溃**——weight-12 下 82% 的 distinct 12-mer 被判「repetitive」（count≥3），weight-14 仍 25%，只有 weight-16 降到 4.5%。即 weight 12–14 的种子几乎不区分重复与单拷贝背景，会用伪候选淹没 Stage 3/4。这与 v3「spaced seed weight ~12–16 攻 C1 高分歧」的初衷直接冲突：

- **种子权重下限锁 ≥16**（spaced seed 的 weight，不是 span）。weight 12–14 不作为发现种子（§A1 的 extend-to-consensus 虽把成本从二次降到线性、算得动，但治不了特异性——伪候选照样进），仅在已知家族内部招募时可临时降权。
- **C1（DNA 已分歧到 16-mer 失效的古老家族）的捕获杠杆改为蛋白域（§C）+ 结构通道（§E3）**，而非压低 DNA 种子权重。理由：分歧到精确 16-mer 都失效的家族，其 DNA 信号已淹在背景里，低权重种子捞回的多是噪声；而蛋白域（编码家族）与结构信号（TIR/LTR/末端，非编码家族）是**正交于 DNA 分歧**的证据，特异性不随分歧塌。
- 若确需 DNA 层非精确种子补 16–25% 分歧带，用**多 spaced pattern @ weight ≥16 + 显著性过滤**，而非降 weight；且只在高拷贝核 mask 后的残余区触发（控伪候选）。

> 实测佐证（DESIGN_REVIEW.md §4.5）：weight-12 的 v3 二次成本 2.94×10⁹（超可行门槛 3×，真实 spaced ×patterns 更高）；§A1 把它降到线性 1.17×10⁸（可行），**但特异性 0.822 仍不可用**——故算力可解不等于该用，权重下限是独立的正确性约束。

---

## B. 封顶与聚类的双输修复（关闭 工程 M1 + 生物 M7）

**B1. C_occ 采样改 hash-based bottom-k**（保留 hash(position) 最小的 200 个，顺序无关、决定性）——同时关闭 工程 M5 的可复现性问题。

**B2. POA 前从全 occurrence 列表重扩成员。** 封顶只用于发现阶段控成本；一旦草簇成形，**从该家族 k-mer 的完整出现列表召回全部成员**再建 consensus，使 Stage 5 不被 200-封顶饿死（关闭 工程 M1 簇跌破 POA 30–50 下限）。

**B3. 古老低频亚家族保护**：reservoir 按比例采样会让 ~1% 占比的古老亚家族（如 AluJ）只摊到 ~2 anchor。v4 在 Stage 1 对**亚家族诊断性 k-mer**（中频、与高频核共线但分歧）单独保底纳入种子集，不被高频核淹没（关闭 生物 M7）。

---

## C. 蛋白通道重构（关闭 生物 M3 + 工程 M3）

**C1. 移出发现循环。** 6-frame 翻译 + HMM 搜索是全基因组 O(L) 重扫，违背有界采样骨架。v4：在**有界样本上跑一次**做局部阈值放宽；在 **Stage 7/8 做一次纯注释/分类扫**。绝不放进迭代循环（关闭 工程 M3）。

**C2. 用蛋白做分类，不只放宽阈值**：RT⇒retro、转座酶⇒DNA-TE、Rep/Hel⇒Helitron、PLE型RT+GIY-YIG⇒Penelope、pPolB+逆转录病毒IN⇒Polinton/Maverick（关闭 生物 M3 的「错失分类增益」）。

**C3. 扩充内置域库**：在 v3 的 RT/RH/IN/GAG/AP/ENV/转座酶/Rep-Hel 上**增补** PLE型RT、GIY-YIG（Penelope）、pPolB、逆转录病毒IN（Polinton/Maverick）。

**C4. 宿主基因守卫。** 「局部放宽阈值」会把多拷贝驯化宿主基因家族（FAR1/FHY3 等 MULE 衍生、CENP-B、DAYSLEEPER）粘进真家族污染 consensus。v4：蛋白命中处放宽阈值的前提是**该位点另有 TE 结构佐证（TSD/TIR/LTR/末端）或邻接已知 TE**；纯蛋白域命中、无结构、且与已注释宿主基因共线者**不放宽、不并簇**（关闭 生物 M3 的污染风险）。

---

## D. 嵌套分解重构（关闭 生物 M5 + 工程 M2）

**D1. freeze-on-accept + 硬 sweep 上限。** 已接纳并掩蔽的家族永不再入聚类；嵌套分解反馈进重聚类设硬 sweep 上限（如 ≤3），消除工程 M2 的振荡/无不动点（关闭 工程 M2）。

**D2. 推广嵌套拓扑**（关闭 生物 M5）：
- 不再要求 host_L、host_R 同家族——允许**异家族宿主**；
- 递归处理**深层嵌套**（插入中插入，玉米基因间区常态）；
- 处理**「宿主只以被打断形式出现」**：当某家族所有拷贝都带同位插入时，先以「带间隔的 profile」建 host consensus（把插入位点记为可变间隔，不计入分歧），再分解；
- 处理**靶位缺失宿主**（插入删掉部分宿主、两翼无法精确重建）：用 profile 比对而非精确拼接判定。

---

## E. 生物学致命缺口补全（关闭 生物 F1/F2/M1/M2/M4/M6）

**E1. 细胞器 / rRNA 假家族守卫（关闭 生物 F1，头号「过度产出」源）。** Stage 8 接纳前，把候选 consensus 比对**细胞器基因组（线粒体+质体）+ rRNA/tRNA/snRNA 模型**；命中者标为 NUMT/NUPT/rDNA/ncRNA-array 并**移出 interspersed 库**。这些序列分散、多拷贝、组成正常、蛋白通道不识别（携带 cox/nad/rbcL），不加此守卫必成假家族。

**E2. 一等 Helitron 处理（关闭 生物 F2）。** 滚环 Helitron 无 TSD/TIR/LTR/polyA，且基因捕获使拷贝 cargo 各异。v4：(i) **Rep/Hel 蛋白命中有权定类为 Helitron**——拟南谱实测验证有效且特异（26/86 Helitron 命中、0/86 Gypsy 假阳；30% 灵敏 = 自主 Helitron 占比，生物学正确），并入 §C 蛋白通道；(ii) 用 5'-TC…CTRR-3' 末端 + 3' 短发夹锚定边界，**必须在基因组实例层（发现期定界）做，不能在已建好的 consensus 上做分类**——实测「consensus 上 hairpin+CTRR 检测」灵敏度=假阳率=7%（无分辨力），因 de novo 边界不精确恰好抹掉末端信号（DESIGN_REVIEW §9.2）；(iii) Helitron 候选**豁免「低分歧低拷贝→疑似 seg-dup」重标**（可即用 Rep/Hel tag 实现）；(iv) 对 cargo 异质拷贝，按末端而非中段聚类，防抗嵌合器撕碎。**开放难题**：~70% 非自主 Helitron（无 RepHel、无干净末端）须靠 (ii) 实例层末端 + (iv) 末端聚类共同解决，是独立增量。
> 方法学（实测归纳）：正交证据（Rep/Hel 域、细胞器同源）在真实数据上稳；依赖 de novo consensus 末端结构的启发式（Helitron 末端、SINE polyA）被不精确边界击穿——能靠正交证据就别靠 consensus 末端结构。

**E3. 非编码分歧家族 rescue（关闭 生物 M1）。** SINE/MITE/LARD/TRIM/SVA 无 ORF，蛋白 graft 救不了。v4 增**非蛋白结构锚定通道**：
- SINE：从 tRNA/7SL/5S 衍生**头部保守模块**做 profile seeding，HMM 迭代招募分歧拷贝；
- MITE：TIR + 短长度 + TSD 招募，专用短元件路径（不被低复杂度过滤吞掉）；
- 末端重复非编码逆转录元件（LARD/TRIM）：按末端直接重复 + 内部低保守招募。

**E4. 大元件与 LTR–internal 关联（关闭 生物 M2）。** ±3 kb 窗口跨不过 gypsy/copia(5–15kb)、Polinton(15–22kb)、全长 L1(6kb)。v4：(i) 窗口**自适应可扩展**（沿种子链延伸，不固定 ±3kb）；(ii) 显式**记录 LTR consensus 与其 internal consensus 的配对关系**（输出 `*-LTR`/`*-I` 并保留 linkage）；(iii) **solo-LTR 模型**：识别由 LTR–LTR 重组产生、与亲本脱钩的孤立 LTR，归回亲本家族而非新建 unknown。

**E5. 超家族感知 TSD + transduction 感知边界（关闭 生物 M6）。**
- TSD 长度按超家族判定：**mariner 的「TA」(2bp) 无信息 → 不作为下调 N_min 的结构佐证**；hAT=8、CACTA=3（且 CACTA 由末端 motif+亚末端重复定义、非长 TIR）等分别处理；
- **transduction 感知**：LINE 3'/SVA 5' transduction 把唯一下游宿主序列带进新拷贝，使侧翼一致性延续到真末端之外 → 当延伸跑进**单拷贝/低重复下游**时截断边界，防 consensus 吞入 transduced 宿主外显子。

**E6. 前置掩蔽改为发现感知（关闭 生物 M4）。** 不在发现前硬掩 SSR/低复杂度（会肢解 SVA、SINE polyA 尾、LTR 内部微卫星、嵌在卫星中的着丝粒逆转录元件、MITE/FB 回文）。v4：低复杂度仅**降权不剔除**，比对器允许穿过；卫星阵列内仍尝试发现嵌入的 interspersed 元件；回文结构（TIR 驱动）不被周期掩蔽器误杀。

---

## F. 其余工程修复

- **F1. Stage 7 用快速 anchored mapper（minimap2 类），成本如实记为 O(L·library_seeds·hit_density)**（关闭 工程 M4）。
- **F2. Stage 8 拷贝数测在全 Stage-7 编目、非 200-样本**（关闭 工程含糊项 6）。
- **F3. 边界调用加成员数门控**：<20 独立侧翼的边缘家族不强行调用侧翼一致性边界，转用结构信号或标记 boundary-uncertain（关闭 工程 M6）。
- **F4. 片段拼接需 junction-spanning 证据**：≥k 个拷贝中 A、B 在单一比对内连续 + 方向一致 + 富集检验，不靠裸相邻（关闭 工程 m1）。
- **F5. 全程决定性**：hash-based 采样 + 规范节点序（按基因组坐标）+ 整数/定点打分 + 确定归约序（关闭 工程 M5）。
- **F6. sketch FP 复核**：Count-Min 估计的 repetitive 候选须精确复核计数再入 R（关闭 工程 m2）。

---

## G. 维持 v3 的正确决策（不改）

发现/编目解耦骨架、profile consensus、迭代残差挖掘 + 递减收益停机、**丢弃 MDL**（两套独立设计 + 审阅一致认为 MDL 非必需且 log2(R) 自相矛盾）、内部显著性 + 结构接纳。

---

## H. 实现前的门禁与构建顺序

**H0（先决，部分已做）：量真实种子成本。** 精确 16-mer 已测（DESIGN_REVIEW.md §4：R=4.17M，真实封顶配对 3.85×10⁸，可行）。**仍缺：真实 spaced-seed（weight 12–16）的配对数**——精确 16-mer 是下界，spaced seed 特异性更低会顶上去（§4.3 caveat 1）。此数定 A1/A2 是否够、是否必须 LSH。→ 见 DESIGN_REVIEW.md §4.5 实测。

构建顺序（先 de-risk）：
1. H0 实测 spaced-seed 成本（本轮做）。
2. Stage 0–2 + A1 extend-to-consensus 最小原型，拟南芥实测真实运行时/内存。
3. Stage 4 决定性聚类 + Stage 5 POA（B2 重扩成员）。
4. Stage 7 编目（复用 minimap2 类）。
5. E1 细胞器守卫、E2 Helitron、E3 非编码 rescue（生物致命三项，优先于 graft）。
6. C 蛋白通道、D 嵌套，最后接入。
7. 全程对照 TAIR10 curated TE 库做家族级 recall/precision。

**注**：「v3/v4 端到端跑通」= 实现上述 1–7，是一个从零的多周系统工程，非现成可运行物。当前可立即做的是 H0 实测与逐步原型，不可伪造完整 pipeline 输出。
