# ============================================================
# SMS - CLASSIFICAÇÃO: Tabela Padrão de Saída
# Base: desbalanceada  |  CV estratificada + SMOTE
# Gera snps_output_2_DESBALANCEADA_SMOTE_ACC_F1.tex
# Dois ambientes: acurácia (SMS_acuracia_2.RData)
#                 F1 Score  (SMS_f1_2.RData)
# ============================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ------------------------------------------------------------
# 1. Carrega os ambientes em objetos separados
# ------------------------------------------------------------
env_acc <- new.env(parent = emptyenv())
env_f1  <- new.env(parent = emptyenv())

load("SMS_acuracia_2.RData", envir = env_acc)
load("SMS_f1_2.RData",       envir = env_f1)

# ------------------------------------------------------------
# 2. Extrai os fitness do GA diretamente dos ambientes salvos
#    GA[[i]]@fitnessValue = aptidão da melhor solução encontrada
#    (calculada com CV estratificada + SMOTE, igual ao pipeline principal)
#    Kernels: [[1]] linear  [[2]] radial 0.001
#             [[3]] radial 0.01  [[4]] radial 0.1  [[5]] radial 1.0
# ------------------------------------------------------------
extract_fitness <- function(env) {
  GA <- env$GA
  vapply(seq_along(GA), function(i) GA[[i]]@fitnessValue, numeric(1))
}

cat("Extraindo fitness do experimento Acurácia...\n")
fitness_acc <- extract_fitness(env_acc)

cat("Extraindo fitness do experimento F1...\n")
fitness_f1 <- extract_fitness(env_f1)

# ------------------------------------------------------------
# 4. Funções auxiliares de formatação
# ------------------------------------------------------------
snps_causais <- paste0("SNP", 1:8)

process_snps <- function(snps_causais, snps_selecionados) {
  nums_caus <- as.integer(sub("SNP", "", snps_causais))
  nums_sel  <- sort(as.integer(sub("SNP", "", snps_selecionados)))
  causais_presentes <- nums_sel[nums_sel %in% nums_caus]

  snp_labels <- sapply(nums_sel, function(n) {
    if (n %in% nums_caus) paste0("\\textbf{", n, "}") else as.character(n)
  })

  # Quebra de linha a cada 5 SNPs
  chunks <- split(snp_labels, ceiling(seq_along(snp_labels) / 5))
  snp_list_str <- paste(sapply(chunks, paste, collapse = ", "), collapse = ",\\\\\n    ")

  list(snp_list_str     = snp_list_str,
       contagem_presentes = length(causais_presentes))
}

fmt_pct <- function(x) sprintf("%.2f", 100 * x)

# ------------------------------------------------------------
# 5. Prepara os resultados de SNPs
# ------------------------------------------------------------
prepare_results <- function(env) {
  sr  <- env$snps_selec_ref
  uni <- env$uniao_final
  int <- env$intersecao_final
  pb  <- env$pvals_bruto_selec
  pa  <- env$pvals_ajustado_selec

  list(
    causais = process_snps(snps_causais, snps_causais),
    r1  = process_snps(snps_causais, sr[[1]]),
    r2  = process_snps(snps_causais, sr[[2]]),
    r3  = process_snps(snps_causais, sr[[3]]),
    r4  = process_snps(snps_causais, sr[[4]]),
    r5  = process_snps(snps_causais, sr[[5]]),
    uni = process_snps(snps_causais, uni),
    int = process_snps(snps_causais, int),
    pb  = process_snps(snps_causais, pb),
    pa  = process_snps(snps_causais, pa),
    sr  = sr, uni_raw = uni, int_raw = int, pb_raw = pb, pa_raw = pa
  )
}

res_acc <- prepare_results(env_acc)
res_f1  <- prepare_results(env_f1)

# Formata os valores de aptidão em percentual
acc_vals <- fmt_pct(fitness_acc)
f1_vals  <- fmt_pct(fitness_f1)

# ------------------------------------------------------------
# 6. Bloco auxiliar para uma linha de tabela
# ------------------------------------------------------------
table_row <- function(method, kernel_str, gamma_str, snp_str, count_str, metric_val) {
  paste0(
    method, " & ", kernel_str, " & ", gamma_str, " &\n",
    "  \\begin{tabular}[c]{@{}c@{}}\n",
    "    ", snp_str, "\n",
    "  \\end{tabular} & ", count_str, " & ", metric_val, " \\\\\n",
    "\\hline\n"
  )
}

count_str <- function(snps_raw, res_obj) {
  paste0(length(snps_raw), " (", res_obj$contagem_presentes, ")")
}

# ------------------------------------------------------------
# 7. Monta o conteúdo LaTeX
# ------------------------------------------------------------
header <- "\\documentclass{article}
\\usepackage[utf8]{inputenc}
\\usepackage[T1]{fontenc}
\\usepackage[brazil]{babel}
\\usepackage{amsmath}
\\usepackage{graphicx}
\\usepackage{array}
\\usepackage{geometry}
\\geometry{a4paper, margin=2cm}

\\begin{document}

"

build_table <- function(res, metric_vals, caption, label, metric_col_header) {
  tab_header <- paste0(
    "% ", paste(rep("-", 70), collapse = ""), "\n",
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\renewcommand{\\arraystretch}{1.2}\n",
    "\\begin{tabular}{c|c|c|c|c|c}\n",
    "\\hline\n",
    "\\textbf{Método} &\n",
    "\\textbf{Kernel} &\n",
    "\\textbf{$\\gamma$} &\n",
    "\\textbf{\\begin{tabular}[c]{@{}c@{}}\\emph{SNPs} selecionados\\end{tabular}} &\n",
    "\\textbf{\\begin{tabular}[c]{@{}c@{}}\\# \\emph{SNPs}\\\\(V)\\end{tabular}} &\n",
    "\\textbf{\\begin{tabular}[c]{@{}c@{}}", metric_col_header, "\\end{tabular}} \\\\\n",
    "\\hline\n"
  )

  rows <- paste0(
    table_row("SNPs causais", "--", "--",
              res$causais$snp_list_str, "--", "--"),
    table_row("SMS", "Linear", "--",
              res$r1$snp_list_str,
              count_str(res$sr[[1]], res$r1), metric_vals[1]),
    table_row("SMS", "Radial", "$0{,}001$",
              res$r2$snp_list_str,
              count_str(res$sr[[2]], res$r2), metric_vals[2]),
    table_row("SMS", "Radial", "$0{,}01$",
              res$r3$snp_list_str,
              count_str(res$sr[[3]], res$r3), metric_vals[3]),
    table_row("SMS", "Radial", "$0{,}1$",
              res$r4$snp_list_str,
              count_str(res$sr[[4]], res$r4), metric_vals[4]),
    table_row("SMS", "Radial", "$1{,}0$",
              res$r5$snp_list_str,
              count_str(res$sr[[5]], res$r5), metric_vals[5]),
    table_row("União", "--", "--",
              res$uni$snp_list_str,
              count_str(res$uni_raw, res$uni), "--"),
    table_row("Interseção", "--", "--",
              res$int$snp_list_str,
              count_str(res$int_raw, res$int), "--"),
    table_row("Valor-$p$ bruto", "--", "--",
              res$pb$snp_list_str,
              count_str(res$pb_raw, res$pb), "--"),
    table_row("Valor-$p$ corrigido", "--", "--",
              res$pa$snp_list_str,
              count_str(res$pa_raw, res$pa), "--")
  )

  tab_footer <- "\\end{tabular}\n\\end{table}\n\n"

  paste0(tab_header, rows, tab_footer)
}

table1 <- build_table(
  res_acc, acc_vals,
  caption = paste0("SVM com kernels Linear e Radial para base de dados 2 desbalanceada\n",
                   "         (SMOTE-NC, acurácia como critério de corte)."),
  label   = "tab:sms_acc_desbal",
  metric_col_header = "Acurácia\\\\média (\\%)"
)

table2 <- build_table(
  res_f1, f1_vals,
  caption = paste0("SVM com kernels Linear e Radial para base de dados 2 desbalanceada\n",
                   "         (SMOTE-NC, F1 Score como critério de corte)."),
  label   = "tab:sms_f1_desbal",
  metric_col_header = "F1\\\\médio (\\%)"
)

legend <- paste0(
  "\\noindent\n",
  "\\textbf{Legenda:} Os números em negrito correspondem aos \\emph{SNPs} causais conhecidos\n",
  "(\\emph{SNPs} 1 a 8). A coluna \\# \\emph{SNPs} (V) indica o total de \\emph{SNPs}\n",
  "selecionados e, entre parênteses, a quantidade de \\emph{SNPs} causais corretamente\n",
  "identificados. A coluna de métrica (acurácia ou F1) refere-se ao valor de aptidão\n",
  "(\\emph{fitness}) obtido pelo Algoritmo Genético com validação cruzada estratificada\n",
  "de 10 \\emph{folds}. As linhas de Valor-$p$ utilizam limiar $\\alpha = 0{,}05$\n",
  "com correção de Bonferroni para os valores ajustados.\n\n",
  "\\end{document}\n"
)

latex_content <- paste0(header, table1, table2, legend)

# ------------------------------------------------------------
# 8. Escreve o arquivo
# ------------------------------------------------------------
out_file <- "snps_output_2_DESBALANCEADA_SMOTE_ACC_F1.tex"
writeLines(latex_content, out_file)
cat("Arquivo '", out_file, "' gerado com sucesso.\n", sep = "")
