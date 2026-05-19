# ============================
# SMS - CLASSIFICAÇÃO (SVM e1071)
# Kernels: Linear e Radial (γ = 0.001, 0.01, 0.1, 1)
# ============================

# Definindo a trilha de dados
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ------------------------------------------------------------------
# Função de validação cruzada p/ classificação (mantém assinatura)
# Retorno: vetor onde o 5º elemento é a ACURÁCIA MÉDIA (para [5])
# ------------------------------------------------------------------
validacao_cruzada_cls <- function(df, folds = 10, gamma = 0.01, cost = 1.0, epsilon = 0.1, kernel = c("linear","radial")) {
  # fenotipo binário como fator
  df$fenotipo <- as.factor(df$fenotipo)
  kernel <- match.arg(kernel)

  n <- nrow(df)
  set.seed(123)
  idx <- sample(rep(1:folds, length.out = n))

  accs <- numeric(folds)
  sens <- numeric(folds)
  spes <- numeric(folds)

  # define a classe positiva (se existir "1", usa "1")
  pos <- if ("1" %in% levels(df$fenotipo)) "1" else levels(df$fenotipo)[1]

  for (k in 1:folds) {
    test_idx  <- which(idx == k)
    train_idx <- which(idx != k)

    train <- df[train_idx, , drop = FALSE]
    test  <- df[test_idx,  , drop = FALSE]

    if (kernel == "linear") {
      fit <- e1071::svm(fenotipo ~ ., data = train,
                        type = "C-classification",
                        kernel = "linear",
                        cost = cost, scale = TRUE)
    } else {
      fit <- e1071::svm(fenotipo ~ ., data = train,
                        type = "C-classification",
                        kernel = "radial",
                        cost = cost, gamma = gamma, scale = TRUE)
    }

    pr <- predict(fit, newdata = test)
    # métricas simples
    tab <- table(pred = pr, obs = test$fenotipo)
    # acurácia
    accs[k] <- sum(diag(tab)) / sum(tab)

    # sensibilidade e especificidade (se ambas as classes aparecerem no fold)
    if (all(c(pos, setdiff(levels(df$fenotipo), pos)[1]) %in% colnames(tab)) &&
        all(c(pos, setdiff(levels(df$fenotipo), pos)[1]) %in% rownames(tab))) {
      TP <- as.numeric(tab[pos, pos])
      FN <- as.numeric(sum(tab[, pos])) - TP
      TN <- as.numeric(tab[setdiff(levels(df$fenotipo), pos)[1],
                           setdiff(levels(df$fenotipo), pos)[1]])
      FP <- as.numeric(sum(tab[setdiff(levels(df$fenotipo), pos)[1], ])) - TN

      sens[k] <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
      spes[k] <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
    } else {
      sens[k] <- NA_real_
      spes[k] <- NA_real_
    }
  }

  acc_mean <- mean(accs, na.rm = TRUE)
  # Retorna 5 posições, com a 5ª sendo a métrica principal (compatível com [5])
  c(acc_mean, sd(accs, na.rm = TRUE), mean(sens, na.rm = TRUE), mean(spes, na.rm = TRUE), acc_mean)
}

# -------------------------------------------
# Parâmetros gerais (mantidos no mesmo estilo)
# -------------------------------------------
cost    <- 1.0
epsilon <- 0.1
folds   <- 10

# -------------------------
# Métricas por cada kernel
# -------------------------
# Linear
kernel <- "linear"
acc_linear <- validacao_cruzada_cls(dados[[1]][c(snps_selec_ref[[1]], "fenotipo")],
                                    folds, gamma = 0.01, cost, epsilon, kernel)[5]

# Radial γ = 0.001
kernel <- "radial"
acc_rbf_g1 <- validacao_cruzada_cls(dados[[1]][c(snps_selec_ref[[2]], "fenotipo")],
                                    folds, gamma = 0.001, cost, epsilon, kernel)[5]

# Radial γ = 0.01
acc_rbf_g2 <- validacao_cruzada_cls(dados[[1]][c(snps_selec_ref[[3]], "fenotipo")],
                                    folds, gamma = 0.01, cost, epsilon, kernel)[5]

# Radial γ = 0.1
acc_rbf_g3 <- validacao_cruzada_cls(dados[[1]][c(snps_selec_ref[[4]], "fenotipo")],
                                    folds, gamma = 0.1, cost, epsilon, kernel)[5]

# Radial γ = 1
acc_rbf_g4 <- validacao_cruzada_cls(dados[[1]][c(snps_selec_ref[[5]], "fenotipo")],
                                    folds, gamma = 1.0, cost, epsilon, kernel)[5]

# -----------------
# Vetores e helpers
# -----------------
snps_causais <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

process_snps <- function(snps_causais, snps_selecionados) {
  numeros_causais <- as.numeric(sub("SNP", "", snps_causais))
  numeros_selecionados <- as.numeric(sub("SNP", "", snps_selecionados))
  numeros_selecionados_ordenados <- sort(numeros_selecionados)
  causais_presentes <- numeros_selecionados_ordenados[numeros_selecionados_ordenados %in% numeros_causais]
  snp_list <- sapply(numeros_selecionados_ordenados, function(num) {
    if (num %in% causais_presentes) paste0("\\textbf{", num, "}") else as.character(num)
  })
  snp_list_str <- paste(snp_list, collapse = ", ")
  snp_list_str <- gsub("(([^,]+, ){4}[^,]+), ", "\\1\\\\\\\\", snp_list_str)
  contagem_presentes <- length(causais_presentes)
  list(snp_list_str = snp_list_str, contagem_presentes = contagem_presentes)
}

# Processando cada lista de SNPs selecionados
result_1  <- process_snps(snps_causais, snps_causais)
result_2  <- process_snps(snps_causais, snps_selec_ref[[1]]) # Linear
result_3  <- process_snps(snps_causais, snps_selec_ref[[2]]) # RBF 0.001
result_4  <- process_snps(snps_causais, snps_selec_ref[[3]]) # RBF 0.01
result_5  <- process_snps(snps_causais, snps_selec_ref[[4]]) # RBF 0.1
result_6  <- process_snps(snps_causais, snps_selec_ref[[5]]) # RBF 1.0
result_8  <- process_snps(snps_causais, uniao_final)
result_9  <- process_snps(snps_causais, intersecao_final)
result_10 <- process_snps(snps_causais, pvals_bruto_selec)
result_11 <- process_snps(snps_causais, pvals_ajustado_selec)

# Formatando a acurácia (em %)
fmt_pct <- function(x) sprintf("%.2f", 100 * x)
acc_linear   <- fmt_pct(acc_linear)
acc_rbf_g1   <- fmt_pct(acc_rbf_g1)
acc_rbf_g2   <- fmt_pct(acc_rbf_g2)
acc_rbf_g3   <- fmt_pct(acc_rbf_g3)
acc_rbf_g4   <- fmt_pct(acc_rbf_g4)

# ---------------------------
# LaTeX (tabela no mesmo estilo)
# ---------------------------
latex_content <- "\\documentclass{article}\n\\usepackage{amsmath}\n\\usepackage{graphicx}\n\\usepackage{array}\n\\begin{document}\n"
latex_content <- paste0(latex_content, "\\begin{table}[h]\n\\caption{SVM com kernels Linear e Radial para base de dados 1 balanceada.}\n\\label{tab:sms_cls}\n")
latex_content <- paste0(latex_content, "\\begin{tabular}{c|c|c|c|c|c}\n\\hline\n")
latex_content <- paste0(latex_content, "\\textbf{Método} & \\textbf{Kernel} & \\textbf{$\\gamma$} & \\textbf{\\begin{tabular}[c]{@{}c@{}}\\emph{SNPs} selecionados\\end{tabular}} & \\textbf{\\# \\emph{SNPs} (V)} & \\textbf{\\begin{tabular}[c]{@{}c@{}}Acurácia\\\\média (\\%)\\end{tabular}} \\\\\n\\hline\n")

add_to_table <- function(method, kernel, gamma, snps, snps_count, acc = "-") {
  latex_content <<- paste0(latex_content, method, " & ", kernel, " & ", gamma,
                           " & \\begin{tabular}[c]{@{}c@{}}", snps, "\\end{tabular} & ",
                           snps_count, " & ", acc, " \\\\\n\\hline\n")
}

# Adicionando as linhas
add_to_table("SNPs causais", "-", "-", result_1$snp_list_str, "-")
add_to_table("SMS", "Linear", "-", result_2$snp_list_str,
             paste(length(snps_selec_ref[[1]]), "(", result_2$contagem_presentes, ")", sep=""), acc_linear)

add_to_table("SMS", "Radial", "0.001", result_3$snp_list_str,
             paste(length(snps_selec_ref[[2]]), "(", result_3$contagem_presentes, ")", sep=""), acc_rbf_g1)

add_to_table("SMS", "Radial", "0.01", result_4$snp_list_str,
             paste(length(snps_selec_ref[[3]]), "(", result_4$contagem_presentes, ")", sep=""), acc_rbf_g2)

add_to_table("SMS", "Radial", "0.1", result_5$snp_list_str,
             paste(length(snps_selec_ref[[4]]), "(", result_5$contagem_presentes, ")", sep=""), acc_rbf_g3)

add_to_table("SMS", "Radial", "1.0", result_6$snp_list_str,
             paste(length(snps_selec_ref[[5]]), "(", result_6$contagem_presentes, ")", sep=""), acc_rbf_g4)

add_to_table("União", "-", "-", result_8$snp_list_str,
             paste(length(uniao_final), "(", result_8$contagem_presentes, ")", sep=""))

add_to_table("Interseção", "-", "-", result_9$snp_list_str,
             paste(length(intersecao_final), "(", result_9$contagem_presentes, ")", sep=""))

add_to_table("Valor-p bruto", "-", "-", result_10$snp_list_str,
             paste(length(pvals_bruto_selec), "(", result_10$contagem_presentes, ")", sep=""))

add_to_table("Valor-p corrigido", "-", "-", result_11$snp_list_str,
             paste(length(pvals_ajustado_selec), "(", result_11$contagem_presentes, ")", sep=""))

latex_content <- paste0(latex_content, "\\end{tabular}\n\\end{table}\n\\end{document}")

# Saída LaTeX
writeLines(latex_content, "snps_output_1_GRUPO_1_CLASSIFICACAO_SVM_Balanceada.tex")
cat("Arquivo 'snps_output_1_GRUPO_1_CLASSIFICACAO_SVM.tex' foi gerado com sucesso.\n")
