# ==============================================================================
# SMS aplicado à Simulação 6 de Oliveira (2015)
# PROTOCOLO: 1 BASE SIMULADA × 10 EXECUÇÕES DO PIPELINE
# ------------------------------------------------------------------------------
# Referência: OLIVEIRA, Fabrízzio Condé de. Um método para seleção de atributos
#   em dados genômicos. Tese (Doutorado em Modelagem Computacional) – UFJF, 2015.
#   Seção 7.1.6 (descrição da Simulação 6) e Seção 8.8 (resultados).
#
# MODELO ORIGINAL (Equação 7.6 de Oliveira 2015):
#   Y = β₀ + β₁·I[SNP1=1] + β₂·I[SNP2=2] + β₃·I[SNP3=3]
#         + β₄·I[SNP4≠1]·I[SNP5=3]
#         + β₅·I[SNP6=1]·I[SNP7=2]·I[SNP8=3] + ε
#
#   β₀=0, β₁=2.0, β₂=1.3, β₃=0.9, β₄=2.0, β₅=3.0
#   n=1000, 100 SNPs, MAF~U[0.1, 0.4]
#   Distribuição das classes (Oliveira 2015): 862 casos (1) × 138 controles (0)
#
# PROTOCOLO DESTA VERSÃO (comparação direta com Oliveira 2015):
#   - A base de dados é simulada UMA ÚNICA VEZ (semente 42), idêntico a Oliveira.
#   - O pipeline SMS (RF → Corte → GA) é executado 10 vezes sobre essa mesma base,
#     com sementes algorítmicas distintas em cada rodada (RF, CV, GA).
#   - Ao final, calcula-se a UNIÃO dos SNPs selecionados pelo GA por kernel,
#     exatamente como Oliveira reportou na Seção 8.8.
#   - A variação entre as 10 execuções é exclusivamente algorítmica (igual a
#     Oliveira), permitindo comparação direta e sem distorção.
#
# INOVAÇÕES em relação a Oliveira (2015):
#   1. SMOTE-NC aplicado exclusivamente no treino de cada fold (sem data leakage)
#   2. Validação cruzada k-fold ESTRATIFICADA (k=10)
#   3. F1-Score da classe minoritária (controles) como métrica de corte e fitness
#   4. SVM binário explícito (C-classification)
#   5. Avaliação dual: nível do indivíduo e da variável
# ==============================================================================

# --- Definir pasta de trabalho ------------------------------------------------
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
cat("[setwd] pasta de trabalho:", getwd(), "\n")

# --- Pacotes ------------------------------------------------------------------
packages <- c("scrime", "e1071", "kernlab", "randomForest",
              "doParallel", "GA", "ggplot2", "recipes", "themis")
installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) install.packages(packages[!installed_packages])
library(scrime)

# ==============================================================================
# SEÇÃO 1 – FUNÇÕES AUXILIARES
# ==============================================================================

.prop_ones_given_beta0 <- function(beta0, n.obs, n.snp, list.snp, list.ia,
                                   beta, maf, err.fun = rnorm, rand = 123,
                                   reps = 1) {
  acc <- 0
  for (i in 1:reps) {
    sim <- simulateSNPglm(n.obs = n.obs, n.snp = n.snp,
                          list.snp = list.snp, list.ia = list.ia,
                          beta0 = beta0, beta = beta, maf = maf,
                          err.fun = err.fun, rand = rand + i)
    acc <- acc + mean(sim$y > 0)
  }
  acc / reps
}

find_beta0_for_target <- function(target_prop, n.obs, n.snp, list.snp, list.ia,
                                  beta, maf, err.fun = rnorm, rand = 123,
                                  reps = 3, lower = -20, upper = 20,
                                  tol = 1e-3, max_expand = 5) {
  if (target_prop <= 0 || target_prop >= 1) stop("target_prop deve estar em (0,1).")
  f <- function(b0) .prop_ones_given_beta0(b0, n.obs, n.snp, list.snp, list.ia,
                                            beta, maf, err.fun, rand, reps) - target_prop
  L <- lower; U <- upper
  fL <- f(L); fU <- f(U)
  k <- 0
  while (sign(fL) == sign(fU) && k < max_expand) {
    L <- L - 10; U <- U + 10; fL <- f(L); fU <- f(U); k <- k + 1
  }
  if (sign(fL) == sign(fU)) stop("Não foi possível fazer bracketing.")
  uniroot(f, interval = c(L, U), tol = tol)$root
}

simulate_scrime_binary_with_ratio <- function(n1, n0, n.snp, list.snp, list.ia,
                                              beta = rep(2, length(list.snp)),
                                              maf = c(0.1, 0.4), err.fun = rnorm,
                                              rand = 123, reps_beta0 = 3,
                                              exact_by_rank = FALSE) {
  n.obs  <- n1 + n0
  target <- n1 / n.obs
  beta0_star <- find_beta0_for_target(target, n.obs, n.snp, list.snp, list.ia,
                                       beta, maf, err.fun, rand, reps_beta0,
                                       lower = -20, upper = 20, tol = 1e-3)
  sim    <- simulateSNPglm(n.obs = n.obs, n.snp = n.snp,
                            list.snp = list.snp, list.ia = list.ia,
                            beta0 = beta0_star, beta = beta,
                            maf = maf, err.fun = err.fun, rand = rand)
  X      <- as.data.frame(sim$x)
  linpred <- sim$y
  prob    <- plogis(linpred)
  if (!exact_by_rank) {
    Y <- as.numeric(linpred > 0)
  } else {
    ord <- order(linpred, decreasing = TRUE)
    Y   <- integer(length(linpred))
    Y[ord[seq_len(n1)]] <- 1
  }
  base <- data.frame(X, linpred = linpred, prob = prob, fenotipo = Y)
  list(data = base, beta0 = beta0_star,
       achieved_prop_1 = mean(Y == 1), counts = table(Y))
}

sanitize_predictors <- function(df, target = "fenotipo") {
  pcols <- setdiff(names(df), target)
  for (nm in pcols) {
    if (is.list(df[[nm]]))    df[[nm]] <- vapply(df[[nm]], function(z) if (length(z) == 1) as.numeric(z) else NA_real_, numeric(1))
    if (is.factor(df[[nm]]))  df[[nm]] <- as.numeric(as.character(df[[nm]]))
    else if (is.character(df[[nm]])) suppressWarnings(df[[nm]] <- as.numeric(df[[nm]]))
  }
  df
}

prepara_cls <- function(df) {
  if ("prob" %in% names(df)) df$prob <- NULL
  df$fenotipo <- factor(df$fenotipo, levels = c(0, 1))
  df
}

aplica_smote_treino <- function(trainset, over_ratio = 1, seed = NULL,
                                mode = c("smotenc", "round")) {
  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)
  ytab <- table(trainset$fenotipo)
  if (length(ytab) < 2L || min(ytab) < 2L) return(trainset)
  pred_cols <- setdiff(names(trainset), "fenotipo")
  if (mode == "smotenc") {
    train_fac <- trainset
    for (nm in pred_cols) {
      v  <- if (is.factor(trainset[[nm]])) as.numeric(as.character(trainset[[nm]])) else trainset[[nm]]
      lv <- sort(unique(stats::na.omit(v)))
      train_fac[[nm]] <- factor(v, levels = as.character(lv))
    }
    rec <- recipes::recipe(fenotipo ~ ., data = train_fac)
    rec <- themis::step_smotenc(rec, fenotipo, over_ratio = over_ratio)
    .err <- NULL
    out <- tryCatch(recipes::bake(recipes::prep(rec, training = train_fac), new_data = NULL),
                    error = function(e) { .err <<- conditionMessage(e); NULL })
    if (is.null(out)) { warning("SMOTE-NC falhou: ", .err, call. = FALSE); return(trainset) }
    out <- as.data.frame(out)
    for (nm in pred_cols) out[[nm]] <- as.integer(as.character(out[[nm]]))
    n_orig <- nrow(trainset)
    out$.smote_sintetico <- c(rep(FALSE, n_orig), rep(TRUE, nrow(out) - n_orig))
  } else {
    rec <- recipes::recipe(fenotipo ~ ., data = trainset)
    rec <- themis::step_smote(rec, fenotipo, over_ratio = over_ratio)
    .err <- NULL
    out <- tryCatch(recipes::bake(recipes::prep(rec, training = trainset), new_data = NULL),
                    error = function(e) { .err <<- conditionMessage(e); NULL })
    if (is.null(out)) { warning("SMOTE (round) falhou: ", .err, call. = FALSE); return(trainset) }
    out  <- as.data.frame(out)
    n_orig <- nrow(trainset)
    out$.smote_sintetico <- c(rep(FALSE, n_orig), rep(TRUE, nrow(out) - n_orig))
    is_syn <- as.logical(out$.smote_sintetico)
    if (any(is_syn)) {
      for (nm in pred_cols) {
        vals_orig <- stats::na.omit(trainset[[nm]])
        uniq_orig <- sort(unique(vals_orig))
        if (length(uniq_orig) == 0L) next
        v_syn <- out[[nm]][is_syn]
        snap  <- vapply(v_syn, function(x) { if (is.na(x)) NA_real_ else uniq_orig[which.min(abs(uniq_orig - x))] }, numeric(1))
        out[[nm]][is_syn] <- snap
      }
    }
  }
  out$.smote_sintetico <- NULL
  sanitize_predictors(prepara_cls(out), target = "fenotipo")
}

validacao_cruzada_cls <- function(data, folds = 10, gamma = 0.01, cost = 1,
                                  kernel = c("radial", "linear"),
                                  use_smote = TRUE, smote_over_ratio = 1,
                                  smote_mode = c("smotenc", "round"),
                                  cv_seed = 123) {
  kernel     <- match.arg(kernel)
  smote_mode <- match.arg(smote_mode)
  data       <- sanitize_predictors(prepara_cls(data), target = "fenotipo")
  n          <- nrow(data)
  set.seed(cv_seed)
  idx     <- numeric(n)
  classes <- levels(as.factor(data$fenotipo))
  for (classe in classes) {
    ci  <- which(data$fenotipo == classe)
    ci2 <- sample(ci)
    fa  <- rep(1:folds, length.out = length(ci2))
    for (j in seq_along(ci2)) idx[ci2[j]] <- fa[j]
  }
  acc <- prec <- rec <- f1 <- prec0 <- rec0 <- f1_0 <- numeric(folds)
  library(e1071)
  for (i in seq_len(folds)) {
    trainset <- na.omit(data[idx != i, , drop = FALSE])
    testset  <- na.omit(data[idx == i, , drop = FALSE])
    trainset <- sanitize_predictors(trainset, target = "fenotipo")
    testset  <- sanitize_predictors(testset,  target = "fenotipo")
    if (isTRUE(use_smote))
      trainset <- aplica_smote_treino(trainset, over_ratio = smote_over_ratio,
                                     seed = cv_seed + 1000L * i, mode = smote_mode)
    X_train <- as.matrix(trainset[, setdiff(names(trainset), "fenotipo"), drop = FALSE])
    y_train <- trainset$fenotipo
    X_test  <- as.matrix(testset[, setdiff(names(testset), "fenotipo"), drop = FALSE])
    y_true  <- factor(testset$fenotipo, levels = c("0", "1"))
    svm_fit <- svm(x = X_train, y = y_train, kernel = kernel, gamma = gamma,
                   cost = cost, type = "C-classification", scale = TRUE)
    y_pred  <- factor(predict(svm_fit, X_test), levels = c("0", "1"))
    cm  <- table(Truth = y_true, Pred = y_pred)
    TP  <- ifelse("1" %in% rownames(cm) && "1" %in% colnames(cm), cm["1","1"], 0)
    TN  <- ifelse("0" %in% rownames(cm) && "0" %in% colnames(cm), cm["0","0"], 0)
    FP  <- ifelse("0" %in% rownames(cm) && "1" %in% colnames(cm), cm["0","1"], 0)
    FN  <- ifelse("1" %in% rownames(cm) && "0" %in% colnames(cm), cm["1","0"], 0)
    acc[i]   <- (TP + TN) / max(1, TP + TN + FP + FN)
    prec[i]  <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
    rec[i]   <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    f1[i]    <- ifelse((prec[i] + rec[i]) > 0, 2 * prec[i] * rec[i] / (prec[i] + rec[i]), 0)
    prec0[i] <- ifelse((TN + FN) > 0, TN / (TN + FN), 0)
    rec0[i]  <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
    f1_0[i]  <- ifelse((prec0[i] + rec0[i]) > 0,
                       2 * prec0[i] * rec0[i] / (prec0[i] + rec0[i]), 0)
  }
  c(acc_mean  = mean(acc),  acc_sd  = sd(acc),
    prec_mean = mean(prec), prec_sd = sd(prec),
    rec_mean  = mean(rec),  rec_sd  = sd(rec),
    f1_mean   = mean(f1),   f1_sd   = sd(f1),
    f1_0_mean = mean(f1_0), f1_0_sd = sd(f1_0))
}

metricas_snps <- function(snps_selecionados, snps_causais = snps_causais_verdadeiros,
                          verbose = TRUE) {
  V <- length(intersect(snps_selecionados, snps_causais))
  T <- length(snps_selecionados)
  C <- length(snps_causais)
  precisao      <- if (T > 0) V / T else 0
  sensibilidade <- if (C > 0) V / C else 0
  f1_var        <- if ((precisao + sensibilidade) > 0)
    2 * precisao * sensibilidade / (precisao + sensibilidade) else 0
  if (verbose) {
    cat(sprintf("  SNPs causais selecionados (V): %d / %d\n", V, C))
    cat(sprintf("  Total selecionados (T)       : %d\n", T))
    cat(sprintf("  Precisão var.    : %.1f%%\n", 100 * precisao))
    cat(sprintf("  Sensibilidade var: %.1f%%\n", 100 * sensibilidade))
    cat(sprintf("  F1 var.          : %.1f%%\n", 100 * f1_var))
  }
  invisible(list(V = V, T = T, C = C, precisao = precisao,
                 sensibilidade = sensibilidade, f1_var = f1_var,
                 snps_causais_selecionados = intersect(snps_selecionados, snps_causais),
                 snps_falso_positivo = setdiff(snps_selecionados, snps_causais)))
}

# ==============================================================================
# SEÇÃO 2 – PARÂMETROS GLOBAIS
# ==============================================================================

library(randomForest); library(e1071); library(GA)
library(ggplot2);       library(recipes); library(themis)

num_snp  <- 100
list.snp <- list(1, 2, 3, c(4, 5), c(6, 7, 8))
list.ia  <- list(1, 2, 3, c(-1, 3), c(1, 2, 3))
beta     <- c(2.0, 1.3, 0.9, 2.0, 3.0)
maf      <- c(0.1, 0.4)
snps_causais_verdadeiros <- paste0("SNP", 1:8)

n1_alvo <- 862
n0_alvo <- 138

kernels_config <- list(
  list(nome = "Linear",       kernel = "linear", gamma = 0.01),
  list(nome = "Radial_0.001", kernel = "radial", gamma = 0.001),
  list(nome = "Radial_0.01",  kernel = "radial", gamma = 0.01),
  list(nome = "Radial_0.1",   kernel = "radial", gamma = 0.1),
  list(nome = "Radial_1",     kernel = "radial", gamma = 1.0)
)
nomes_kernels <- sapply(kernels_config, `[[`, "nome")

run_ga <- 30; maxiter_ga <- 10; pcross <- 0.8; pmut <- 0.1
elitism <- 5; popSize <- 100

percentual_snps <- 0.95
passo           <- 10

n_exec     <- 10
# Sementes algorítmicas distintas por execução (RF, CV, GA)
# A base é fixa; apenas os algoritmos variam.
seeds_exec <- 1:n_exec

t_total_inicio <- proc.time()

# ==============================================================================
# SEÇÃO 3 – SIMULAR A BASE UMA ÚNICA VEZ (protocolo Oliveira 2015)
# ==============================================================================

cat("=== Simulando base da Simulação 6 (semente fixa = 42) ===\n")
cat(sprintf("  Alvo: %d casos / %d controles (razão 1:%.2f)\n",
            n1_alvo, n0_alvo, n1_alvo / n0_alvo))
set.seed(42)
t_sim_ini <- proc.time()
res_sim6 <- simulate_scrime_binary_with_ratio(
  n1 = n1_alvo, n0 = n0_alvo,
  n.snp    = num_snp,
  list.snp = list.snp, list.ia = list.ia,
  beta     = beta, maf = maf,
  rand     = 42,
  reps_beta0    = 5,
  exact_by_rank = TRUE
)
dados_base <- subset(res_sim6$data, select = c(-linpred, -prob))
t_sim <- round((proc.time() - t_sim_ini)[["elapsed"]])
cat(sprintf("  β₀ calibrado: %.4f\n", res_sim6$beta0))
cat("  Distribuição obtida:\n"); print(res_sim6$counts)
cat(sprintf("  [TIMER] Simulação: %.1f s\n\n", t_sim))
write.csv(dados_base, file = "dados_sim6_base_unica.csv", row.names = FALSE)

limite <- floor((num_snp / passo) * percentual_snps)

# ==============================================================================
# SEÇÃO 4 – ARMAZENAMENTO
# ==============================================================================

snps_GA_all   <- vector("list", n_exec)  # [[exec]][[kernel_nome]]
metricas_exec <- vector("list", n_exec)
corte_exec    <- vector("list", n_exec)

t_rf_vec    <- numeric(n_exec)
t_corte_mat <- matrix(0, nrow = n_exec, ncol = length(nomes_kernels),
                      dimnames = list(NULL, nomes_kernels))
t_ga_mat    <- matrix(0, nrow = n_exec, ncol = length(nomes_kernels),
                      dimnames = list(NULL, nomes_kernels))

# ==============================================================================
# SEÇÃO 5 – LOOP: 10 EXECUÇÕES DO PIPELINE SOBRE A MESMA BASE
# ==============================================================================

for (exec in seq_len(n_exec)) {
  seed_exec <- seeds_exec[exec]
  cat(sprintf("\n%s\n", strrep("=", 70)))
  cat(sprintf("  EXECUÇÃO %d / %d  (semente algorítmica: %d | base: fixa)\n",
              exec, n_exec, seed_exec))
  cat(sprintf("%s\n", strrep("=", 70)))

  # ---- 5.1 Random Forest (nova semente a cada execução) --------------------
  cat(sprintf("  Executando Random Forest (semente=%d)...\n", seed_exec))
  data_rf <- dados_base
  if (!is.factor(data_rf$fenotipo))
    data_rf$fenotipo <- factor(data_rf$fenotipo, levels = c(0, 1))
  set.seed(seed_exec)
  t_rf_ini <- proc.time()
  RF      <- randomForest(fenotipo ~ ., data = data_rf, ntree = 4000,
                          mtry = ncol(data_rf) - 1, importance = TRUE)
  imp     <- importance(RF)
  mdg_col <- if ("MeanDecreaseGini"       %in% colnames(imp)) "MeanDecreaseGini"       else colnames(imp)[ncol(imp)]
  mda_col <- if ("MeanDecreaseAccuracy"   %in% colnames(imp)) "MeanDecreaseAccuracy"   else colnames(imp)[1]
  rank_RF     <- sort(imp[, mdg_col], decreasing = TRUE)   # gVI
  rank_RF_pVI <- sort(imp[, mda_col], decreasing = TRUE)   # pVI

  # Salvar ambos os índices na primeira execução (base fixa — mesmos dados)
  if (exec == 1) {
    df_imp <- data.frame(
      SNP      = rownames(imp),
      gVI      = round(imp[, mdg_col], 3),
      pVI      = round(imp[, mda_col], 3),
      rank_gVI = rank(-imp[, mdg_col]),
      rank_pVI = rank(-imp[, mda_col]),
      causal   = rownames(imp) %in% snps_causais_verdadeiros,
      stringsAsFactors = FALSE
    )
    df_imp <- df_imp[order(df_imp$rank_gVI), ]
    write.csv(df_imp, file.path(pasta_saida, "comparacao_gVI_pVI.csv"), row.names = FALSE)
    cat(sprintf("  [INFO] Importâncias gVI/pVI salvas em comparacao_gVI_pVI.csv\n"))
    cat(sprintf("  [INFO] Top-3 pVI: %s\n", paste(names(rank_RF_pVI)[1:3], collapse = ", ")))
    cat(sprintf("  [INFO] Corr. Spearman gVI vs pVI: %.3f\n",
                cor(df_imp$rank_gVI, df_imp$rank_pVI, method = "spearman")))
  }

  t_rf_vec[exec] <- round((proc.time() - t_rf_ini)[["elapsed"]])
  cat(sprintf("  [TIMER] RF exec %d: %.1f s | Top-3 gVI: %s\n",
              exec, t_rf_vec[exec],
              paste(names(rank_RF)[1:3], collapse = ", ")))

  snps_GA_all[[exec]]   <- list()
  metricas_exec[[exec]] <- list()
  corte_exec[[exec]]    <- list()

  # ---- 5.2 Loop de kernels --------------------------------------------------
  for (ki in seq_along(kernels_config)) {
    cfg    <- kernels_config[[ki]]
    nome   <- cfg$nome
    kernel <- cfg$kernel
    gamma  <- cfg$gamma
    folds  <- 10
    cost   <- 1.0

    cat(sprintf("\n  --- Exec %d | Kernel %s (gamma=%.3f) ---\n",
                exec, nome, gamma))

    # -- Corte ----------------------------------------------------------------
    mean_f1_curve <- numeric()
    cv_seed       <- seed_exec * 1000 + ki * 100
    t_corte_ini   <- proc.time()
    for (cont in 1:(limite + 1)) {
      j       <- cont * passo
      var_sel <- head(names(rank_RF), j)
      base_cv <- dados_base[c(var_sel, "fenotipo")]
      svm_cv  <- validacao_cruzada_cls(base_cv, folds = folds, gamma = gamma,
                                       cost = cost, kernel = kernel,
                                       cv_seed = cv_seed)
      mean_f1_curve[cont] <- svm_cv["f1_0_mean"]
      cat(sprintf("    Corte grupo %d (%d SNPs): F1_min=%.4f\n",
                  cont, j, mean_f1_curve[cont]))
    }

    # Gráfico de corte (somente exec 1)
    if (exec == 1) {
      pdf(sprintf("Grafico_F1_minoria_SVM_%s_exec01.pdf", nome), height = 5, width = 9)
      xs <- seq(passo, (limite + 1) * passo, by = passo)
      plot(xs, mean_f1_curve, type = "o", lwd = 2,
           xlab = "Grupo de Marcadores",
           ylab = sprintf("F1 Minoria/Controles (SVM %s, 10-fold)", nome),
           main = sprintf("Oliveira 2015 – Sim6 (1 base, exec 1) | Kernel %s", nome))
      maximo_idx  <- which.max(mean_f1_curve)
      buffer_alvo <- (maximo_idx + 1) * passo
      pt_corte    <- min(buffer_alvo, length(rank_RF))
      abline(v = pt_corte, lty = 2, lwd = 1.5)
      dev.off()
    }

    maximo_idx  <- which.max(mean_f1_curve)
    buffer_alvo <- (maximo_idx + 1) * passo
    pt_corte    <- min(buffer_alvo, length(rank_RF))
    corte_exec[[exec]][[nome]] <- pt_corte
    t_corte_mat[exec, nome] <- round((proc.time() - t_corte_ini)[["elapsed"]])
    cat(sprintf("    [TIMER] Corte: %.1f s | Ponto de corte: %d SNPs\n",
                t_corte_mat[exec, nome], pt_corte))

    snps_corte     <- head(names(rank_RF), pt_corte)
    genotipo_corte <- dados_base[, snps_corte, drop = FALSE]
    dados2         <- cbind(genotipo_corte, fenotipo = dados_base$fenotipo)
    if (!is.factor(dados2$fenotipo))
      dados2$fenotipo <- factor(dados2$fenotipo, levels = c(0, 1))

    # -- GA -------------------------------------------------------------------
    ga_seed <- seed_exec * 100 + ki
    t_ga_ini <- proc.time()
    f <- local({
      .kernel <- kernel; .gamma <- gamma; .cost <- cost; .folds <- folds
      .dados2 <- dados2; .cv_seed <- cv_seed
      function(x) {
        inc <- which(x == 1)
        if (length(inc) == 0) return(0)
        dv <- .dados2[, c(inc, ncol(.dados2)), drop = FALSE]
        colnames(dv)[ncol(dv)] <- "fenotipo"
        res <- validacao_cruzada_cls(dv, folds = .folds, gamma = .gamma,
                                     cost = .cost, kernel = .kernel,
                                     cv_seed = .cv_seed)
        res["f1_0_mean"]
      }
    })

    set.seed(ga_seed)
    ga_result <- ga(
      type       = "binary",
      fitness    = f,
      nBits      = ncol(genotipo_corte),
      popSize    = popSize,
      names      = colnames(genotipo_corte),
      maxiter    = maxiter_ga,
      seed       = ga_seed,
      parallel   = TRUE,
      run        = run_ga,
      pcrossover = pcross,
      pmutation  = pmut,
      elitism    = elitism,
      suggestions = matrix(rep(1, ncol(genotipo_corte)), nrow = 1)
    )

    sol_mat  <- as.matrix(ga_result@solution)
    sel_cols <- which(sol_mat[1, ] == 1)
    snps_ga  <- colnames(sol_mat)[sel_cols]
    snps_GA_all[[exec]][[nome]] <- snps_ga
    t_ga_mat[exec, nome] <- round((proc.time() - t_ga_ini)[["elapsed"]])
    cat(sprintf("    [TIMER] GA: %.1f s | SNPs selecionados: %d\n",
                t_ga_mat[exec, nome], length(snps_ga)))

    m_corte <- metricas_snps(snps_corte, verbose = FALSE)
    m_ga    <- metricas_snps(snps_ga,    verbose = FALSE)
    metricas_exec[[exec]][[nome]] <- list(corte = m_corte, ga = m_ga)

    # Gráfico GA (somente exec 1)
    if (exec == 1) {
      pdf(sprintf("Grafico_GA_F1_%s_exec01.pdf", nome), height = 5, width = 9)
      geracao <- seq_len(ga_result@iter)
      df_ga   <- data.frame(
        Geracao     = rep(geracao, 3),
        Aptidao     = c(ga_result@summary[, 2],
                        ga_result@summary[, 4],
                        ga_result@summary[, 1]),
        Estatistica = rep(c("Média", "Mediana", "Melhor"), each = length(geracao))
      )
      print(
        ggplot(df_ga, aes(x = Geracao, y = Aptidao, group = Estatistica)) +
          geom_line(aes(colour = Estatistica, linetype = Estatistica), linewidth = 1.2) +
          geom_point(size = 1.5) +
          scale_x_continuous(breaks = geracao) +
          labs(title  = sprintf("GA – F1 Minoria | %s | Exec 1 | 1 base × 10 exec", nome),
               x = "Geração", y = "F1-Score classe 0 (fitness)") +
          theme_bw()
      )
      dev.off()
    }
  }  # fim loop kernels

  cat(sprintf("\n  Exec %d concluída.\n", exec))
}  # fim loop execuções

# ==============================================================================
# SEÇÃO 6 – AGREGAÇÃO: UNIÃO POR KERNEL (10 execuções, mesma base)
# ==============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("  RESUMO FINAL — UNIÃO POR KERNEL (10 execuções, 1 base fixa)\n")
cat("  Protocolo equivalente a Oliveira (2015), Seção 8.8\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("  Base: dados_sim6_base_unica.csv (semente 42)\n"))
cat(sprintf("  SNPs causais verdadeiros: %s\n\n",
            paste(snps_causais_verdadeiros, collapse = ", ")))

resultados_uniao <- data.frame(
  Kernel     = character(),
  N_uniao    = integer(),
  V_uniao    = integer(),
  Prec_uniao = numeric(),
  Sens_uniao = numeric(),
  F1_uniao   = numeric(),
  stringsAsFactors = FALSE
)

snps_uniao_por_kernel    <- list()
snps_intersec_por_kernel <- list()

for (nome in nomes_kernels) {
  snps_lista   <- lapply(seq_len(n_exec), function(e) snps_GA_all[[e]][[nome]])

  # Salvar listas individuais por execução (para cálculo de interseção)
  snps_exec_df <- do.call(rbind, lapply(seq_along(snps_lista), function(e) {
    data.frame(Exec = e, Kernel = nome, SNP = snps_lista[[e]],
               stringsAsFactors = FALSE)
  }))
  nome_arq_exec <- file.path(pasta_saida,
    paste0("SNPs_por_exec_", gsub("\\.", "", nome), "_1base10exec.csv"))
  write.csv(snps_exec_df, nome_arq_exec, row.names = FALSE)

  # União e interseção por kernel
  snps_uniao_k    <- Reduce(union,     snps_lista)
  snps_intersec_k <- Reduce(intersect, snps_lista)
  snps_uniao_por_kernel[[nome]]    <- snps_uniao_k
  snps_intersec_por_kernel[[nome]] <- snps_intersec_k

  m <- metricas_snps(snps_uniao_k)
  cat(sprintf("  Kernel: %s\n", nome))
  resultados_uniao <- rbind(resultados_uniao, data.frame(
    Kernel     = nome,
    N_uniao    = m$T,
    V_uniao    = m$V,
    Prec_uniao = round(100 * m$precisao, 1),
    Sens_uniao = round(100 * m$sensibilidade, 1),
    F1_uniao   = round(100 * m$f1_var, 1),
    stringsAsFactors = FALSE
  ))

  # Log da interseção por kernel
  cat(sprintf("  Interseção (%s): %d SNPs\n", nome, length(snps_intersec_k)))
  if (length(snps_intersec_k) > 0) metricas_snps(snps_intersec_k)
  cat("\n")
}

cat("\n--- Tabela resumo (União 10 execuções por kernel, 1 base) ---\n")
print(resultados_uniao, row.names = FALSE)

# Salvar interseções por kernel em CSV
for (nome in nomes_kernels) {
  snps_int <- snps_intersec_por_kernel[[nome]]
  df_int <- data.frame(
    SNP    = if (length(snps_int) > 0) snps_int else character(0),
    causal = if (length(snps_int) > 0) snps_int %in% snps_causais_verdadeiros else logical(0),
    stringsAsFactors = FALSE
  )
  nome_arq_int <- file.path(pasta_saida,
    paste0("SNPs_intersecao_", gsub("\\.", "", nome), "_1base10exec.csv"))
  write.csv(df_int, nome_arq_int, row.names = FALSE)
}

snps_super_uniao    <- Reduce(union,     snps_uniao_por_kernel)
snps_super_intersec <- Reduce(intersect, snps_uniao_por_kernel)

cat("\n--- Super-União (todos os kernels) ---\n")
metricas_snps(snps_super_uniao)

cat("\n--- Super-Interseção (todos os kernels) ---\n")
if (length(snps_super_intersec) > 0) {
  metricas_snps(snps_super_intersec)
} else {
  cat("  (conjunto vazio)\n")
}

# ==============================================================================
# SEÇÃO 7 – DETALHAMENTO POR EXECUÇÃO
# ==============================================================================

cat("\n--- Detalhamento por execução (GA) ---\n")
resultados_por_exec <- do.call(rbind, lapply(seq_len(n_exec), function(exec) {
  do.call(rbind, lapply(nomes_kernels, function(nome) {
    m <- metricas_exec[[exec]][[nome]]$ga
    data.frame(Exec = exec, Kernel = nome,
               N_GA = m$T, V_GA = m$V,
               Prec_GA = round(100 * m$precisao, 1),
               Sens_GA = round(100 * m$sensibilidade, 1),
               F1_GA   = round(100 * m$f1_var, 1),
               stringsAsFactors = FALSE)
  }))
}))
print(resultados_por_exec, row.names = FALSE)

# ==============================================================================
# SEÇÃO 8 – SALVAR CSVs
# ==============================================================================

write.csv(resultados_uniao,
          file = "Resultados_Uniao_1base10exec.csv", row.names = FALSE)
write.csv(resultados_por_exec,
          file = "Resultados_Detalhados_1base10exec.csv", row.names = FALSE)
write.csv(
  data.frame(SNP    = snps_super_uniao,
             causal = snps_super_uniao %in% snps_causais_verdadeiros),
  file = "SNPs_super_uniao_1base10exec.csv", row.names = FALSE)
write.csv(
  data.frame(SNP    = snps_super_intersec,
             causal = snps_super_intersec %in% snps_causais_verdadeiros),
  file = "SNPs_super_intersecao_1base10exec.csv", row.names = FALSE)
for (nome in nomes_kernels) {
  snps_k <- snps_uniao_por_kernel[[nome]]
  write.csv(
    data.frame(SNP    = snps_k,
               causal = snps_k %in% snps_causais_verdadeiros),
    file = sprintf("SNPs_uniao_%s_1base10exec.csv", nome), row.names = FALSE)
}

# ==============================================================================
# SEÇÃO 9 – TEMPOS DE EXECUÇÃO
# ==============================================================================

t_total <- round((proc.time() - t_total_inicio)[["elapsed"]])

cat("\n===== Tempos de Execução (segundos) =====\n")
cat(sprintf("  %-32s: %6.1f s\n", "Simulação (1 vez)", t_sim))
cat(sprintf("  %-32s: %6.1f s (média: %.1f s)\n",
            "Random Forest (10 exec)", sum(t_rf_vec), mean(t_rf_vec)))
cat("\n  -- Corte por kernel (soma 10 exec) --\n")
for (nome in nomes_kernels)
  cat(sprintf("  Corte %-22s: %6.1f s (média: %.1f s)\n",
              nome, sum(t_corte_mat[, nome]), mean(t_corte_mat[, nome])))
cat("\n  -- GA por kernel (soma 10 exec) --\n")
for (nome in nomes_kernels)
  cat(sprintf("  GA    %-22s: %6.1f s (média: %.1f s)\n",
              nome, sum(t_ga_mat[, nome]), mean(t_ga_mat[, nome])))
cat(sprintf("\n  %-32s: %6.1f s (%.1f min)\n",
            "TOTAL GERAL", t_total, t_total / 60))
cat("Arquivos gerados na pasta:", getwd(), "\n")
