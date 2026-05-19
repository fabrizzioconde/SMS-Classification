# ==============================================================================
# Gera tabela LaTeX equivalente à Tabela 8.22 de Oliveira (2015)
# a partir dos resultados da execução 1 base × 10 execuções algorítmicas.
# Substitui AUC (Oliveira) por F1-Score da classe minoritária (controles).
# ==============================================================================

.script_dir <- tryCatch({
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path)) {
    dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    .args     <- commandArgs(trailingOnly = FALSE)
    .file_arg <- sub("^--file=", "", .args[grep("^--file=", .args)])
    if (length(.file_arg) > 0) dirname(normalizePath(.file_arg)) else getwd()
  }
}, error = function(e) getwd())
setwd(.script_dir)

# Ler CSV de resultados detalhados por execução
df <- read.csv("Resultados_Detalhados_1base10exec.csv", stringsAsFactors = FALSE)

# Ler CSV de resultados por kernel (união)
df_uniao <- read.csv("Resultados_Uniao_1base10exec.csv", stringsAsFactors = FALSE)

kernels <- c("Linear", "Radial_0.001", "Radial_0.01", "Radial_0.1", "Radial_1")
kernel_labels <- c(
  "Linear"       = "Linear",
  "Radial_0.001" = "$\\gamma=0{,}001$",
  "Radial_0.01"  = "$\\gamma=0{,}01$",
  "Radial_0.1"   = "$\\gamma=0{,}1$",
  "Radial_1"     = "$\\gamma=1$"
)

# Ler F1_min no ponto de corte a partir do log
# (extraído do log: F1_min no grupo do ponto de corte ótimo por exec/kernel)
log_lines <- readLines("execucao_Sim6_1base10exec.log")

# Extrair F1_min no ponto de corte ótimo para cada exec e kernel
f1_corte_mat <- matrix(NA_real_, nrow = 10, ncol = length(kernels),
                       dimnames = list(1:10, kernels))

exec_atual <- NA
kernel_atual <- NA
corte_grupos <- c()
kernel_map <- c(
  "Linear"       = "Linear",
  "Radial_0.001" = "Radial_0.001",
  "Radial_0.01"  = "Radial_0.01",
  "Radial_0.1"   = "Radial_0.1",
  "Radial_1"     = "Radial_1"
)

for (line in log_lines) {
  # Detectar execução
  m <- regmatches(line, regexpr("EXECUÇÃO (\\d+) / 10", line))
  if (length(m) > 0) {
    exec_atual <- as.integer(sub("EXECUÇÃO (\\d+) / 10", "\\1", m))
    corte_grupos <- c()
  }
  # Detectar kernel
  for (k in kernels) {
    if (grepl(paste0("Kernel ", k), line)) {
      kernel_atual <- k
      corte_grupos <- c()
      break
    }
  }
  # Coletar F1_min dos grupos de corte
  m2 <- regmatches(line, regexpr("F1_min=([0-9.]+)", line))
  if (length(m2) > 0) {
    corte_grupos <- c(corte_grupos, as.numeric(sub("F1_min=([0-9.]+)", "\\1", m2)))
  }
  # Detectar ponto de corte e registrar F1_min do grupo ótimo
  m3 <- regmatches(line, regexpr("Ponto de corte: (\\d+) SNPs", line))
  if (length(m3) > 0 && !is.na(exec_atual) && !is.na(kernel_atual)) {
    pt_corte   <- as.integer(sub("Ponto de corte: (\\d+) SNPs", "\\1", m3))
    idx_grupo  <- pt_corte / 10 - 1   # grupo anterior ao ponto de corte
    if (idx_grupo >= 1 && idx_grupo <= length(corte_grupos)) {
      f1_corte_mat[exec_atual, kernel_atual] <- corte_grupos[idx_grupo]
    } else if (length(corte_grupos) > 0) {
      f1_corte_mat[exec_atual, kernel_atual] <- max(corte_grupos, na.rm = TRUE)
    }
  }
}

# Montar data.frame por execução e kernel com SNPs(V) e F1_min
cat("\n=== Tabela equivalente à Tabela 8.22 de Oliveira (2015) ===\n")
cat("Protocolo: 1 base × 10 execuções algorítmicas\n")
cat("Métrica: F1-Score da classe minoritária (controles, classe 0) no ponto de corte\n\n")

# Cabeçalho da tabela
cat(sprintf("%-6s", "Iter"))
for (k in kernels)
  cat(sprintf("  %-14s %-8s", "SNPs (V)", "F1 min"))
cat("\n")

means_snps <- setNames(numeric(length(kernels)), kernels)
means_f1   <- setNames(numeric(length(kernels)), kernels)
sds_snps   <- setNames(numeric(length(kernels)), kernels)
sds_f1     <- setNames(numeric(length(kernels)), kernels)
snps_vals  <- setNames(vector("list", length(kernels)), kernels)
f1_vals    <- setNames(vector("list", length(kernels)), kernels)
for (k in kernels) { snps_vals[[k]] <- c(); f1_vals[[k]] <- c() }

for (exec in 1:10) {
  cat(sprintf("%-6d", exec))
  for (k in kernels) {
    row_k <- df[df$Exec == exec & df$Kernel == k, ]
    if (nrow(row_k) == 0) {
      cat(sprintf("  %-14s %-8s", "---", "---"))
      next
    }
    n_ga <- row_k$N_GA
    v_ga <- row_k$V_GA
    f1   <- if (!is.na(f1_corte_mat[exec, k])) round(f1_corte_mat[exec, k], 3) else NA
    cat(sprintf("  %4d (%-2d)      %-8s", n_ga, v_ga, ifelse(is.na(f1), "---", sprintf("%.3f", f1))))
    snps_vals[[k]] <- c(snps_vals[[k]], n_ga)
    f1_vals[[k]]   <- c(f1_vals[[k]],   ifelse(is.na(f1), NA_real_, f1))
  }
  cat("\n")
}

cat(sprintf("%-6s", "Média"))
for (k in kernels) {
  mn_snps <- round(mean(snps_vals[[k]], na.rm=TRUE), 1)
  mn_v    <- round(mean(sapply(1:10, function(e) {
    r <- df[df$Exec==e & df$Kernel==k,]; if(nrow(r)>0) r$V_GA else NA }), na.rm=TRUE), 1)
  mn_f1   <- round(mean(f1_vals[[k]], na.rm=TRUE), 3)
  cat(sprintf("  %4.1f (%-3.1f)   %-8s", mn_snps, mn_v, sprintf("%.3f", mn_f1)))
}
cat("\n")

cat(sprintf("%-6s", "σ"))
for (k in kernels) {
  sd_snps <- round(sd(snps_vals[[k]], na.rm=TRUE), 1)
  sd_v    <- round(sd(sapply(1:10, function(e) {
    r <- df[df$Exec==e & df$Kernel==k,]; if(nrow(r)>0) r$V_GA else NA }), na.rm=TRUE), 1)
  sd_f1   <- round(sd(f1_vals[[k]], na.rm=TRUE), 3)
  cat(sprintf("  %4.1f (%-3.1f)   %-8s", sd_snps, sd_v, sprintf("%.3f", sd_f1)))
}
cat("\n\n")

# --- Gerar código LaTeX -------------------------------------------------------
cat("\n% ---- CÓDIGO LATEX ----\n")

linha_hdr <- paste0(
  "\\begin{tabular}{c|",
  paste(rep("rr", length(kernels)), collapse = "|"),
  "}\n\\hline\n"
)

# Linha de nomes dos kernels (multicolumn 2)
nomes_mc <- paste(sapply(kernels, function(k) {
  sprintf("\\multicolumn{2}{c|}{%s}", kernel_labels[k])
}), collapse = " & ")
linha_hdr <- paste0(linha_hdr, "\\textbf{Iter} & ", nomes_mc, " \\\\\n")

# Sub-cabeçalho SNPs(V) | F1 min
sub_hdr <- paste(rep("SNPs (V) & F1 min", length(kernels)), collapse = " & ")
linha_hdr <- paste0(linha_hdr, " & ", sub_hdr, " \\\\\n\\hline\n")

cat(linha_hdr)

for (exec in 1:10) {
  linha <- sprintf("%d", exec)
  for (k in kernels) {
    row_k <- df[df$Exec == exec & df$Kernel == k, ]
    if (nrow(row_k) == 0) { linha <- paste0(linha, " & --- & ---"); next }
    n_ga <- row_k$N_GA; v_ga <- row_k$V_GA
    f1   <- if (!is.na(f1_corte_mat[exec, k])) round(f1_corte_mat[exec, k], 3) else NA
    # Negrito na melhor exec por kernel (maior F1)
    cel_snps <- sprintf("%d (%d)", n_ga, v_ga)
    cel_f1   <- if (is.na(f1)) "---" else sprintf("%.3f", f1)
    linha <- paste0(linha, " & ", cel_snps, " & ", cel_f1)
  }
  cat(paste0(linha, " \\\\\n"))
}

cat("\\hline\n")

# Linha de média
linha_media <- "Média"
linha_sd    <- "$\\sigma$"
for (k in kernels) {
  mn_snps <- mean(snps_vals[[k]], na.rm=TRUE)
  mn_v    <- mean(sapply(1:10, function(e) {
    r <- df[df$Exec==e & df$Kernel==k,]; if(nrow(r)>0) r$V_GA else NA }), na.rm=TRUE)
  mn_f1   <- mean(f1_vals[[k]], na.rm=TRUE)
  sd_snps <- sd(snps_vals[[k]], na.rm=TRUE)
  sd_v    <- sd(sapply(1:10, function(e) {
    r <- df[df$Exec==e & df$Kernel==k,]; if(nrow(r)>0) r$V_GA else NA }), na.rm=TRUE)
  sd_f1   <- sd(f1_vals[[k]], na.rm=TRUE)
  linha_media <- paste0(linha_media,
    sprintf(" & %.1f (%.1f) & %.3f", mn_snps, mn_v, mn_f1))
  linha_sd <- paste0(linha_sd,
    sprintf(" & %.1f (%.1f) & %.3f", sd_snps, sd_v, sd_f1))
}
cat(paste0(linha_media, " \\\\\n"))
cat(paste0(linha_sd,    " \\\\\n"))
cat("\\hline\n\\end{tabular}\n")

cat("\nTabela gerada com sucesso.\n")
