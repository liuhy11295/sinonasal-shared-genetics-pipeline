# PROJECT MEMORY

## Current Project

CRSwNP-centered shared genetic architecture, pQTL/mQTL MR, bulk transcriptomics, single-cell analysis, spatial transcriptomics, network pharmacology and virtual perturbation project.

中文定位：

以慢性鼻-鼻窦炎伴鼻息肉（CRSwNP / nasal polyps）为中心，先挖掘其与哮喘、过敏性鼻炎、IgE、嗜酸粒细胞和其他免疫相关疾病的共享遗传架构，再用 pQTL/mQTL MR 筛选蛋白质和代谢物，最后通过 bulk、单细胞、空间转录组、网络药理学和虚拟扰动形成纯生信闭环。

## Main Files

1. `sinonasal_mr_spatial_virtual_cell_stepwise_protocol.md`
   - 最新、最完整的逐步执行方案。
   - 后续分析应以此文件为主。

2. `mr_multiomics_spatial_virtual_cell_plan.md`
   - MR 多组学-单细胞-空间-虚拟细胞方案初版。

3. `ent_external_database_feasibility.md`
   - 外部数据库可行性，包含 GEO、TCGA、空间转录组等。

4. `nhanes_ent_topic_screening.md`
   - NHANES 耳鼻喉方向筛选记录。
   - 当前 NHANES 已降级为可选补充，不是主证据链。

5. `comprehensive_bioinformatics_analysis_plan.md`
   - 早期综合生信分析大方案。

6. `zotero_bioinformatics_method_innovation_report.md`
   - 最早 Zotero 文献方法创新总结。

## Key Decisions

1. 研究重点是纯生信，不做湿实验。
2. 主疾病锚点是 CRSwNP / nasal polyps。
3. 甲状腺疾病排除，不属于本项目耳鼻喉范围。
4. NHANES 不作为主线；只有当候选分子有清晰代理变量时才作为补充。
5. 空间转录组疾病场景固定为 CRS/CRSwNP，尤其是鼻息肉组织。
6. 共病不是装饰，而是上游机制分层工具。
7. MR 结局必须分层：
   - 主结局：nasal polyps / CRSwNP、CRS、allergic rhinitis。
   - 共病结局：只纳入 Phase 0 筛出来共享强的表型，例如 asthma。
   - 内表型：IgE、eosinophil、neutrophil、CRP、cytokines。
8. 最终药物筛选不是宣称预防共病，而是寻找能打到 CRSwNP 与共病共享炎症通路的候选靶点或药物。

## Final Workflow

```text
Phase 0  CRSwNP-centered shared genetic architecture
        ↓
Phase 1  pQTL/mQTL MR screening
        ↓
Phase 2  colocalization / reverse MR / replication
        ↓
Phase 3  bulk transcriptome validation
        ↓
Phase 4  single-cell localization, pseudotime, communication, regulons
        ↓
Phase 5  spatial transcriptomics niche validation
        ↓
Phase 6  network pharmacology + signature reversal + virtual perturbation
        ↓
Phase 7  integrated prioritization
```

## Expected Final Results

1. CRSwNP-centered immune comorbidity genetic map。
2. 有 MR、colocalization 和 reverse MR 支持的候选蛋白/代谢物。
3. CRSwNP bulk 组织层面的通路和模块验证。
4. 单细胞层面的来源细胞、作用细胞、拟时序、通讯轴和 regulon。
5. 空间层面的 immune-epithelial-stromal remodeling niche。
6. 网络药理学、signature reversal 和虚拟扰动支持的候选治疗轴或再利用药物。
7. 最终综合优先级表。

## Current Next Step

从 Phase 0 开始，建立 `gwas_inventory.tsv`：

1. nasal polyps / CRSwNP。
2. chronic rhinosinusitis。
3. allergic rhinitis / hay fever。
4. asthma。
5. total IgE。
6. eosinophil count。
7. neutrophil count / CRP / cytokines。
8. 免疫共病扩展表型，例如 eczema、food allergy、urticaria、IBD、psoriasis、RA、SLE 等。

优先任务：

1. 确认每个 GWAS summary statistics 的 accession、release、ancestry、sample size、case/control、genome build、download URL。
2. 锁定主分析版本。
3. 再进入 LDSC/LAVA shared genetic architecture。

## How To Resume In Codex

下次继续时，对 Codex 说：

```text
工作目录切到 D:\生信。先读取 PROJECT_MEMORY.md 和 sinonasal_mr_spatial_virtual_cell_stepwise_protocol.md，然后从 Phase 0 的 gwas_inventory.tsv 开始继续。
```

