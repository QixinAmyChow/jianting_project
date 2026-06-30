# gsea_families_shared.R
# Unified semantic-dedup family list used by both
#   gsea_alldb_stacking.R (slide 12, 05e)
#   gsea_cl1_combined.R   (slide 14, 05f)
# First-match-wins: specific patterns first, generic ones last.
# Sourced from jianting_project/: source("gsea_families_shared.R")

shared_families <- list(
  # ── HALLMARK-specific concepts (Hallmark-only biological processes) ──────
  H_MYC          = "MYC_TARGET|HALLMARK_MYC",
  H_MTORC1       = "MTORC1",
  H_OXPHOS       = "OXIDATIVE_PHOSPH",
  H_GLYC         = "GLYCOLYSIS",
  H_HYPOXIA      = "HYPOXIA",
  H_APOP         = "APOPTOSIS",
  H_G2M          = "G2M_CHECKPOINT",
  H_E2F          = "E2F_TARGET",
  H_ANDROGEN     = "ANDROGEN",
  H_ESTROGEN     = "ESTROGEN",
  H_P53          = "P53|TP53",
  H_PI3K         = "PI3K|AKT_MTOR",
  H_COAGUL       = "COAGULATION",
  H_COMPLEM      = "COMPLEMENT",
  H_FATTY        = "FATTY_ACID",
  H_CHOLEST      = "CHOLESTEROL|BILE",
  H_ANGIO        = "ANGIOGENESIS",
  H_UNFOLDED     = "UNFOLDED_PROTEIN",
  H_DNA_REP      = "DNA_REPAIR",
  H_MITOTIC      = "MITOTIC_SPINDLE",
  H_UV           = "UV_RESPONSE",
  H_XENO         = "XENOBIOTIC",
  H_INFLAM       = "INFLAMMATORY_RESPONSE",

  # ── C6 oncogenic ──────────────────────────────────────────────────────────
  C6_KRAS        = "KRAS",
  C6_PRC1        = "BMI1|MEL18|PRC1",
  C6_PRC2        = "PRC2|SUZ12|EZH2",
  C6_ESC         = "ESC_V6|ESC_J",
  C6_MORF        = "MORF_",
  C6_STK33       = "STK33",
  C6_AKT         = "^AKT",
  C6_CYCLIN      = "CYCLIN",
  C6_ERBB        = "ERBB2",
  C6_EGFR        = "EGFR",

  # ── ECM / Structural ──────────────────────────────────────────────────────
  COLLAGEN       = "COLLAGEN",
  ECM_ORG        = "ECM_ORGAN|EXTRACELLULAR_MATRIX_ORGAN|EXTRACELLULAR_MATRIX_ASSEMB|ECM_ASSEMB|ECM_REMODEL|EXTRACELLULAR_STRUCT",
  MMP            = "MMP|MATRIX_METALLOPROTEIN|METALLOENDOPEPTI",
  BASEMENT_MEM   = "BASEMENT_MEMBRANE",
  INTEGRIN       = "INTEGRIN",
  FIBRONECTIN    = "FIBRONECT",
  LAMININ        = "LAMININ",
  ECM_GEN        = "EXTRACELLULAR_MATRIX",

  # ── TGF-b / EMT ───────────────────────────────────────────────────────────
  TGFB           = "TGFB|TGF.?B|TRANSFORMING_GROWTH|SMAD",
  EMT            = "EPITHELIAL_MESENCH|EMT_TRANS|MESENCHYMAL_TRANS|\\bEMT\\b",

  # ── Wnt / osteogenic niche ────────────────────────────────────────────────
  WNT            = "WNT|BETA_CATENIN",
  OSTEOBLAST     = "OSTEOBLAST|OSTEOGENIC_DIFF|BONE_FORM",
  OSTEOCLAST     = "OSTEOCLAST|BONE_RESORPT",
  OSSIFICATION   = "OSSIF|BONE_MINERAL|BONE_REMODEL|BONE_DEVELOP",
  BMP            = "\\bBMP\\b|BONE_MORPHOGEN",
  NOTCH          = "NOTCH",
  HEDGEHOG       = "HEDGEHOG",

  # ── Immunosuppression / T cell ────────────────────────────────────────────
  CD8_EXHAUST    = "CD8.*EXHAUST|EXHAUST.*CD8|T_CELL_EXHAUST",
  TREG           = "TREG|REGULATORY_T_CELL|T_REGULATORY",
  CHECKPOINT     = "CHECKPOINT|\\bPD.?1\\b|\\bPDL1\\b|CTLA4|LAG3|TIM3",
  IL10           = "\\bIL10\\b|\\bIL_10\\b",
  T_CELL_ACTIV   = "T_CELL_ACTIV|T_CELL_STIMUL|T_CELL_PROLIF|T_CELL_RECEPTOR_SIGNAL",
  CD8_GEN        = "\\bCD8\\b",
  CD4_GEN        = "\\bCD4\\b",
  T_CELL_GEN     = "T_CELL|TCELL",

  # ── Macrophage polarization ───────────────────────────────────────────────
  PPARG          = "PPARG|PPAR.?GAMMA",
  WBP7           = "WBP7",
  MDSC           = "MDSC|MYELOID_DERIVED_SUPPRESS",
  M2_POLAR       = "M2_POLARIZ|ALTERNATIVE_ACTIV|ALTERNATIVELY_ACTIV|ANTI.INFLAM.*MAC",
  M1_POLAR       = "M1_POLARIZ|CLASSIC_ACTIV|CLASSICALLY_ACTIV|INFLAM.*MAC",
  TAM            = "\\bTAM\\b|TUMOR.*MACRO|MACRO.*TUMOR|TUMOR_ASSOC.*MACRO",
  MYELOID_DIFF   = "MYELOID_DIFF|MYELOID_CELL_DIFF|MYELOID_LEUK_DIFF",
  MONOCYTE       = "MONOCYTE",
  MYELOID_GEN    = "MYELOID",
  MACRO_GEN      = "MACROPHAGE|\\bMACRO\\b",

  # ── Cytokine / NFkB / JAK-STAT / interferons ──────────────────────────────
  # TNFA before NFKB so HALLMARK_TNFA_SIGNALING_VIA_NFKB stays in TNFA.
  TNFA           = "\\bTNF\\b|TNF_ALPHA|TNFA",
  NFKB           = "NFKB|NF.KB|NFKAPPAB",
  JAK_STAT       = "JAK.?STAT|JAK_STAT",
  STAT3          = "STAT3",
  STAT_GEN       = "\\bSTAT\\b",
  IFN_GAMMA      = "INTERFERON_GAMMA|IFN_GAMMA|IFNG|TYPE_II_IFN",
  IFN_ALPHA      = "INTERFERON_ALPHA|IFN_ALPHA|IFNA|TYPE_I_IFN",
  IFN_GEN        = "INTERFERON|\\bIFN\\b",
  IL6            = "\\bIL6\\b|\\bIL_6\\b|INTERLEUKIN_6",
  IL4_IL13       = "\\bIL4\\b|\\bIL13\\b|IL_4|IL_13",
  CHEMOKINE      = "CHEMOKINE",
  CYTOKINE_PROD  = "CYTOKINE_PRODUCT|CYTOKINE_SECRET|CYTOKINE_BIOSYN",
  CYTOKINE_SIG   = "CYTOKINE_SIGNAL|CYTOKINE_MEDIAT|CYTOKINE_NETWORK",
  CYTOKINE_GEN   = "CYTOKINE"
)
