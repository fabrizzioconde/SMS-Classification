# ==============================================================================
# Computa gVI e pVI para as 3 bases do TCC e gera tabelas LaTeX
# Base 1, Base 2 (dados.csv) e Simulação 6 (dados_sim6_base_unica.csv)
# ==============================================================================
suppressPackageStartupMessages(library(randomForest))

raiz <- "D:/Dropbox/SMS - Classificação/Códigos em R/Códigos usados no TCC/TCC - Izabela"

# ------------------------------------------------------------------------------
# Função: treina RF e retorna tabela de importância (gVI + pVI)
# ------------------------------------------------------------------------------
calc_importancia <- function(dados_path, causais, seed = 42, ntree = 4000) {
  d <- read.csv(dados_path, stringsAsFactors = FALSE)

  # Selecionar colunas SNP (prefixo "data.SNP" ou "SNP") + resposta (data.fenotipo ou fenotipo)
  snp_cols  <- grep("^data\\.SNP|^SNP", colnames(d), value = TRUE)
  resp_col  <- if ("data.fenotipo" %in% colnames(d)) "data.fenotipo" else "fenotipo"
  df <- d[, c(snp_cols, resp_col)]
  # Renomear: remover prefixo "data."
  colnames(df) <- gsub("^data\\.", "", colnames(df))
  df$fenotipo <- as.factor(df$fenotipo)
  col_resposta <- "fenotipo"

  p <- length(snp_cols)
  n <- nrow(df)
  cat(sprintf("  n=%d, p=%d SNPs | Distribuição: %s\n",
              n, p, paste(table(df[[col_resposta]]), collapse=" / ")))

  set.seed(seed)
  rf <- randomForest(as.formula(paste(col_resposta, "~ .")),
                     data = df, ntree = ntree,
                     importance = TRUE, keep.forest = FALSE)

  imp     <- importance(rf)
  gvi_col <- "MeanDecreaseGini"
  pvi_col <- "MeanDecreaseAccuracy"
  snp_names <- gsub("^data\\.", "", rownames(imp))

  df_imp <- data.frame(
    SNP      = snp_names,
    gVI      = round(imp[, gvi_col], 3),
    pVI      = round(imp[, pvi_col], 3),
    rank_gVI = as.integer(rank(-imp[, gvi_col])),
    rank_pVI = as.integer(rank(-imp[, pvi_col])),
    causal   = snp_names %in% causais,
    stringsAsFactors = FALSE
  )
  df_imp <- df_imp[order(df_imp$rank_gVI), ]

  sp <- cor(df_imp$rank_gVI, df_imp$rank_pVI, method = "spearman")
  cat(sprintf("  Correlação Spearman gVI vs pVI (p=%d SNPs): %.3f\n\n", p, sp))

  df_imp
}

# ------------------------------------------------------------------------------
# Função: imprime tabela apenas com SNPs causais
# ------------------------------------------------------------------------------
tabela_causais <- function(df_imp, causais) {
  dc <- df_imp[df_imp$causal, ]
  dc <- dc[order(dc$rank_gVI), ]
  print(dc[, c("SNP","gVI","rank_gVI","pVI","rank_pVI")], row.names=FALSE)
  cat("\n")
  invisible(dc)
}

# ------------------------------------------------------------------------------
# Função: gera código LaTeX da tabela
# ------------------------------------------------------------------------------
latex_causais <- function(dc, caption, label) {
  cat(sprintf("\n%%---- LaTeX: %s ----\n", label))
  cat("\\begin{table}[H]\n\\centering\n")
  cat(sprintf("\\caption{\\novo{%s}}\n", caption))
  cat(sprintf("\\label{%s}\n", label))
  cat("\\small\n\\renewcommand{\\arraystretch}{1.3}\n")
  cat("\\begin{tabular}{L{2.2cm}|C{1.5cm}|C{1.5cm}|C{1.5cm}|C{1.5cm}}\n")
  cat("\\hline\n")
  cat("\\textbf{SNP} & \\textbf{gVI} & \\textbf{Rank gVI} & \\textbf{pVI} & \\textbf{Rank pVI} \\\\\n")
  cat("\\hline\n")
  for (i in seq_len(nrow(dc))) {
    r <- dc[i, ]
    cat(sprintf("%s & %.3f & %d & %.3f & %d \\\\\n",
                r$SNP, r$gVI, r$rank_gVI,
                r$pVI, r$rank_pVI))
  }
  cat("\\hline\n")
  cat(paste0("\\multicolumn{5}{L{10cm}}{\\footnotesize\\textit{",
             "gVI = \\textit{Mean Decrease in Gini}; pVI = \\textit{Mean Decrease in Accuracy}. ",
             "Rank calculado sobre todos os 100 SNPs. ",
             "RF com \\texttt{ntree=4000}, \\texttt{set.seed(42)}, base bruta sem SMOTE-NC.}} \\\\\n"))
  cat("\\hline\n")
  cat("\\end{tabular}\n\\end{table}\n")
}

# ==============================================================================
# BASE 1 — Balanceada (500/500)
# ==============================================================================
cat("=== BASE 1 BALANCEADA ===\n")
causais_b1 <- paste0("SNP", 1:8)
df1bal <- calc_importancia(
  dados_path = file.path(raiz, "SMS_Completo_Acuracia_1_Base_Balanceada_Validacao_Balanceada/dados.csv"),
  causais    = causais_b1
)
dc1bal <- tabela_causais(df1bal, causais_b1)
latex_causais(dc1bal,
  caption = "gVI e pVI dos 8 SNPs causais da Base de Dados~1 (efeitos aditivos independentes, base balanceada 500/500, $p=100$, \\texttt{ntree=4000}).",
  label   = "tab:gvi_pvi_base1_bal")

# ==============================================================================
# BASE 1 — Desbalanceada (800/200)
# ==============================================================================
cat("=== BASE 1 DESBALANCEADA ===\n")
df1des <- calc_importancia(
  dados_path = file.path(raiz, "SMS_Completo_Acuracia_1_Desbalanceada_CV_estrat_SMOTE/dados.csv"),
  causais    = causais_b1
)
dc1des <- tabela_causais(df1des, causais_b1)
latex_causais(dc1des,
  caption = "gVI e pVI dos 8 SNPs causais da Base de Dados~1 (efeitos aditivos independentes, base desbalanceada 800/200, $p=100$, \\texttt{ntree=4000}).",
  label   = "tab:gvi_pvi_base1")

# ==============================================================================
# BASE 2 — Balanceada (500/500)
# ==============================================================================
cat("=== BASE 2 BALANCEADA ===\n")
causais_b2 <- paste0("SNP", 1:8)
df2bal <- calc_importancia(
  dados_path = file.path(raiz, "SMS_Completo_Acuracia_2_Balanceada_Validacao_Balanceada/dados.csv"),
  causais    = causais_b2
)
dc2bal <- tabela_causais(df2bal, causais_b2)
latex_causais(dc2bal,
  caption = "gVI e pVI dos 8 SNPs causais da Base de Dados~2 (efeitos agrupados e interações epistáticas, base balanceada 500/500, $p=100$, \\texttt{ntree=4000}).",
  label   = "tab:gvi_pvi_base2_bal")

# ==============================================================================
# BASE 2 — Desbalanceada (800/200)
# ==============================================================================
cat("=== BASE 2 DESBALANCEADA ===\n")
df2des <- calc_importancia(
  dados_path = file.path(raiz, "SMS_Completo_Acuracia_2_Desbalanceada_CV_estrat_SMOTE/dados.csv"),
  causais    = causais_b2
)
dc2des <- tabela_causais(df2des, causais_b2)
latex_causais(dc2des,
  caption = "gVI e pVI dos 8 SNPs causais da Base de Dados~2 (efeitos agrupados e interações epistáticas, base desbalanceada 800/200, $p=100$, \\texttt{ntree=4000}).",
  label   = "tab:gvi_pvi_base2")

# ==============================================================================
# SIMULAÇÃO 6 — SNP1-SNP8 (usa comparacao_gVI_pVI.csv já gerado)
# ==============================================================================
cat("=== SIMULAÇÃO 6 ===\n")
csv6 <- file.path(raiz, "SMS_Oliveira2015_Sim6_Desbal_SMOTE_F1_1base10exec/comparacao_gVI_pVI.csv")
if (file.exists(csv6)) {
  df6 <- read.csv(csv6, stringsAsFactors = FALSE)
  cat("  [OK] lido comparacao_gVI_pVI.csv\n")
  cat(sprintf("  Correlação Spearman: %.3f\n\n",
              cor(df6$rank_gVI, df6$rank_pVI, method="spearman")))
  dc6 <- df6[df6$causal, ]
  dc6 <- dc6[order(dc6$rank_gVI), ]
  print(dc6[, c("SNP","gVI","rank_gVI","pVI","rank_pVI")], row.names=FALSE)
  cat("\n")
  latex_causais(dc6,
    caption = "gVI e pVI dos 8 SNPs causais da Simulação~6 de \\citeonline{Oliveira2015} (efeitos aditivos e epistáticos, base desbalanceada 862/138, $p=100$, \\texttt{ntree=4000}).",
    label   = "tab:gvi_pvi_sim6")
} else {
  cat("  [!] comparacao_gVI_pVI.csv não encontrado. Execute o script 1base10exec.\n")
}
