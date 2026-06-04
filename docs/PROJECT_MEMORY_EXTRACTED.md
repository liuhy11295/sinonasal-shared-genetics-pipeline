# Project Memory Extracted

## Active project

CRSwNP / nasal polyps centered pure bioinformatics project.

Core evidence chain:

```text
shared genetic architecture
-> pQTL/mQTL MR
-> colocalization / reverse MR / replication
-> bulk transcriptome validation
-> single-cell localization, pseudotime, communication, regulons
-> spatial transcriptomics niche validation
-> network pharmacology + signature reversal + virtual perturbation
-> integrated prioritization
```

## Main source files

1. `sinonasal_mr_spatial_virtual_cell_stepwise_protocol.md`
   - Current main protocol. Future analysis should follow this file first.
2. `mr_multiomics_spatial_virtual_cell_plan.md`
   - Earlier MR multi-omics / single-cell / spatial plan.
3. `ent_external_database_feasibility.md`
   - External database feasibility notes.
4. `nhanes_ent_topic_screening.md`
   - NHANES screening notes. NHANES is optional support, not the main evidence chain.
5. `comprehensive_bioinformatics_analysis_plan.md`
   - Earlier broad bioinformatics plan.
6. `zotero_bioinformatics_method_innovation_report.md`
   - Zotero-derived method innovation summary.

## Key decisions

1. The project is pure bioinformatics; no wet-lab experiment is required.
2. The disease anchor is CRSwNP / nasal polyps.
3. Thyroid disease is excluded from the ENT scope.
4. NHANES is not a main line; it is only supplementary when a candidate molecule has a clear measurable proxy.
5. Spatial transcriptomics disease context should be CRS/CRSwNP, especially nasal polyp tissue.
6. Comorbidity analysis is a mechanistic stratification tool, not decoration.
7. MR outcomes must be stratified:
   - Main outcomes: nasal polyps / CRSwNP, CRS, allergic rhinitis.
   - Comorbidity outcomes: only Phase 0-supported shared phenotypes, such as asthma.
   - Endophenotypes: IgE, eosinophils, neutrophils, CRP, cytokines.
8. Drug prioritization should target shared inflammatory pathways between CRSwNP and supported comorbidities.

## Current next step

Start Phase 0 by building `gwas_inventory.tsv`.

Priority metadata for each GWAS summary statistics dataset:

```text
accession, release, ancestry, sample size, case/control count,
genome build, download URL, analysis role, notes
```

Initial target phenotypes:

1. nasal polyps / CRSwNP
2. chronic rhinosinusitis
3. allergic rhinitis / hay fever
4. asthma
5. total IgE
6. eosinophil count
7. neutrophil count / CRP / cytokines
8. immune comorbidity expansion phenotypes, e.g. eczema, food allergy, urticaria, IBD, psoriasis, RA, SLE

