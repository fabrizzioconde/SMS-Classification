# ==============================================================================
# SMS aplicado à Simulação 6 de Oliveira (2015)
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
# INOVAÇÕES DESTE EXPERIMENTO em relação a Oliveira (2015):
#   1. SMOTE-NC aplicado exclusivamente no treino de cada fold (sem data leakage)
#   2. Validação cruzada k-fold ESTRATIFICADA (k=10), mantendo proporção de
#      classes em cada fold — crucial com 862/138 casos/controles.
#   3. F1-Score como métrica de corte E fitness do GA (no lugar de AUC-ROC)
#   4. SVM binário explícito (C-classification) em vez de SVR adaptado
#   5. Avaliação dual ao nível do indivíduo (F1 obs.) e da variável (precisão,
#      sensibilidade e F1 dos SNPs selecionados vs. causais verdadeiros)
# ==============================================================================

# --- Definir pasta de trabalho automaticamente --------------------------------
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

# --- Timers -------------------------------------------------------------------
t_total_inicio <- proc.time()
t_kernels      <- list()
t_corte        <- list()  # tempo da etapa de corte por kernel
t_ga           <- list()  # tempo do GA por kernel
t_sim          <- 0       # tempo da simulação da base
t_rf           <- 0       # tempo da Random Forest

# --- Pacotes ------------------------------------------------------------------
packages <- c("scrime", "e1071", "kernlab", "randomForest",
              "doParallel", "GA", "ggplot2", "recipes", "themis")
installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) install.packages(packages[!installed_packages])
library(scrime)

# ==============================================================================
# SEÇÃO 1 – SIMULAÇÃO DA BASE (replicando Oliveira 2015, Simulação 6)
# ==============================================================================

# --- Funções auxiliares de simulação ------------------------------------------
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

# --- Parâmetros da Simulação 6 de Oliveira (2015) ----------------------------
# SNPs causais: 1, 2, 3 (aditivos isolados) + {4,5} (interação ord.2) + {6,7,8} (ord.3)
num_snp  <- 100
list.snp <- list(1, 2, 3, c(4, 5), c(6, 7, 8))
list.ia  <- list(1, 2, 3, c(-1, 3), c(1, 2, 3))
beta     <- c(2.0, 1.3, 0.9, 2.0, 3.0)   # efeitos originais de Oliveira (2015)
maf      <- c(0.1, 0.4)
snps_causais_verdadeiros <- paste0("SNP", 1:8)  # nomes scrime: SNP1..SNP8

# Distribuição de classes: 862 casos / 138 controles (Oliveira 2015, Seção 8.8)
n1_alvo <- 862
n0_alvo <- 138

cat("=== Simulando base da Simulação 6 (Oliveira 2015) ===\n")
cat(sprintf("  Alvo: %d casos / %d controles  (razão 1:%.1f)\n",
            n1_alvo, n0_alvo, n0_alvo / n1_alvo))
set.seed(42)
t_sim_ini <- proc.time()
res_sim6 <- simulate_scrime_binary_with_ratio(
  n1 = n1_alvo, n0 = n0_alvo,
  n.snp    = num_snp,
  list.snp = list.snp, list.ia = list.ia,
  beta     = beta, maf = maf,
  rand     = 42,
  reps_beta0    = 5,
  exact_by_rank = TRUE     # força exatamente 862/138 por ranqueamento do logito
)
cat(sprintf("  β₀ calibrado: %.4f\n", res_sim6$beta0))
cat("  Distribuição de classes obtida:\n"); print(res_sim6$counts)

# Base de trabalho (sem linpred e prob)
dados      <- list()
dados[[1]] <- subset(res_sim6$data, select = c(-linpred, -prob))
write.csv(dados[[1]], file = "dados_sim6_oliveira2015.csv", row.names = FALSE)
t_sim <- round((proc.time() - t_sim_ini)[["elapsed"]])
cat(sprintf("[TIMER] Simulação     : %8.1f s\n", t_sim))
cat("  Base salva em dados_sim6_oliveira2015.csv\n\n")

# ==============================================================================
# SEÇÃO 2 – FUNÇÕES DO PIPELINE SMS
# ==============================================================================

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
                                mode = c("smotenc", "round"),
                                save_path = NULL) {
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
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(out, save_path, row.names = FALSE)
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
  # Nota: f1 = F1 da classe 1 (maioria); f1_0 = F1 da classe 0 (minoria = controles).
  # Na Sim6, classe 1 tem 86% dos dados, então F1_1 é inflado por trivialidade.
  # Usamos f1_0 como métrica de corte e fitness do GA para ter curva discriminante.
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
    # Métricas da classe 1 (casos – maioria)
    prec[i]  <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
    rec[i]   <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    f1[i]    <- ifelse((prec[i] + rec[i]) > 0, 2 * prec[i] * rec[i] / (prec[i] + rec[i]), 0)
    # Métricas da classe 0 (controles – MINORIA): TP_0=TN, FP_0=FN, FN_0=FP
    prec0[i] <- ifelse((TN + FN) > 0, TN / (TN + FN), 0)
    rec0[i]  <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
    f1_0[i]  <- ifelse((prec0[i] + rec0[i]) > 0,
                       2 * prec0[i] * rec0[i] / (prec0[i] + rec0[i]), 0)
  }
  c(acc_mean  = mean(acc),  acc_sd  = sd(acc),
    prec_mean = mean(prec), prec_sd = sd(prec),
    rec_mean  = mean(rec),  rec_sd  = sd(rec),
    f1_mean   = mean(f1),   f1_sd   = sd(f1),
    # F1 da classe minoritária (controles) – usado como métrica de corte e GA
    f1_0_mean = mean(f1_0), f1_0_sd = sd(f1_0))
}

# Calcula métricas ao nível da variável (SNPs)
metricas_snps <- function(snps_selecionados, snps_causais = snps_causais_verdadeiros) {
  V <- length(intersect(snps_selecionados, snps_causais))  # verdadeiros positivos
  T <- length(snps_selecionados)                           # total selecionados
  C <- length(snps_causais)                                # total causais verdadeiros
  precisao    <- if (T > 0) V / T else 0
  sensibilidade <- if (C > 0) V / C else 0
  f1_var      <- if ((precisao + sensibilidade) > 0)
    2 * precisao * sensibilidade / (precisao + sensibilidade) else 0
  cat(sprintf("  SNPs causais selecionados (V): %d / %d\n", V, C))
  cat(sprintf("  Total selecionados (T)       : %d\n", T))
  cat(sprintf("  Precisão var.    : %.1f%%\n", 100 * precisao))
  cat(sprintf("  Sensibilidade var: %.1f%%\n", 100 * sensibilidade))
  cat(sprintf("  F1 var.          : %.1f%%\n", 100 * f1_var))
  invisible(list(V = V, T = T, C = C, precisao = precisao,
                 sensibilidade = sensibilidade, f1_var = f1_var,
                 snps_causais_selecionados = intersect(snps_selecionados, snps_causais),
                 snps_falso_positivo = setdiff(snps_selecionados, snps_causais)))
}

# ==============================================================================
# SEÇÃO 3 – PIPELINE SMS (RF → Corte → GA, por kernel)
# ==============================================================================

library(randomForest); library(e1071); library(GA)
library(ggplot2);       library(recipes); library(themis)

mean_svm_RF_list <- list()
GA_list          <- list()
corte            <- list()
snps_selec_ref   <- list()
metricas_var     <- list()
percentual_snps  <- 0.95

# --- Random Forest (único, independe do kernel SVM) ---------------------------
cat("=== Etapa 1: Random Forest (ntree=4000) ===\n")
data_temp <- dados[[1]]
if (!is.factor(data_temp$fenotipo)) data_temp$fenotipo <- factor(data_temp$fenotipo, levels = c(0, 1))
set.seed(1)
t_rf_ini <- proc.time()
RF      <- randomForest(fenotipo ~ ., data = data_temp, ntree = 4000,
                        mtry = ncol(data_temp) - 1, importance = TRUE)
imp     <- importance(RF)
mdg_col <- if ("MeanDecreaseGini" %in% colnames(imp)) "MeanDecreaseGini" else colnames(imp)[ncol(imp)]
rank_RF <- sort(imp[, mdg_col], decreasing = TRUE)
t_rf <- round((proc.time() - t_rf_ini)[["elapsed"]])
cat(sprintf("[TIMER] Random Forest : %8.1f s\n", t_rf))
cat(sprintf("  SNPs no ranking: %d\n", length(rank_RF)))
cat("  Top-10 do rank RF:\n"); print(head(rank_RF, 10))

# Parâmetros GA (comuns a todos os kernels)
run_ga  <- 30; maxiter_ga <- 10; pcross <- 0.8; pmut <- 0.1
elitism <- 5;  popSize    <- 100

passo   <- 10
limite  <- floor((ncol(dados[[1]]) - 1) / passo * percentual_snps)

# ============================================================
# KERNEL DEFINITIONS – loop sobre 5 kernels
# ============================================================
kernels_config <- list(
  list(nome = "Linear",       kernel = "linear", gamma = 0.01),
  list(nome = "Radial_0.001", kernel = "radial", gamma = 0.001),
  list(nome = "Radial_0.01",  kernel = "radial", gamma = 0.01),
  list(nome = "Radial_0.1",   kernel = "radial", gamma = 0.1),
  list(nome = "Radial_1",     kernel = "radial", gamma = 1.0)
)

for (i in seq_along(kernels_config)) {
  cfg    <- kernels_config[[i]]
  nome   <- cfg$nome
  kernel <- cfg$kernel
  gamma  <- cfg$gamma
  folds  <- 10
  cost   <- 1.0

  cat(sprintf("\n=== Kernel %s (gamma=%.3f) ===\n", nome, gamma))
  t_k_ini <- proc.time()

  # --- Etapa 2: Corte (curva F1_minoria na trilha RF) --------------------------
  # Usamos f1_0_mean (F1 da classe 0 = controles, minoria) para ter curva
  # discriminante; F1_1 (maioria) é inflado e não mostra pico claro.
  mean_svm_RF_list[[i]] <- numeric()
  t_corte_ini <- proc.time()
  for (cont in 1:(limite + 1)) {
    j       <- cont * passo
    var_sel <- if (j == 1) names(rank_RF)[1] else head(names(rank_RF), j)
    base_cv <- dados[[1]][c(var_sel, "fenotipo")]
    svm_cv  <- validacao_cruzada_cls(base_cv, folds = folds, gamma = gamma,
                                     cost = cost, kernel = kernel)
    mean_svm_RF_list[[i]][cont] <- svm_cv["f1_0_mean"]
    cat(sprintf("  Corte grupo %d (%d SNPs): F1_minoria=%.4f\n",
                cont, j, mean_svm_RF_list[[i]][cont]))
  }

  pdf(sprintf("Grafico_F1_minoria_SVM_%s.pdf", nome), height = 5, width = 9)
  xs <- seq(passo, limite * passo + 10, by = passo)
  plot(xs, mean_svm_RF_list[[i]], type = "o", lwd = 2,
       xlab = "Grupo de Marcadores",
       ylab = sprintf("F1 Minoria/Controles (SVM %s, 10-fold)", nome),
       main = sprintf("Oliveira 2015 – Sim6 | Kernel %s | F1 classe 0", nome))
  maximo_idx  <- which.max(mean_svm_RF_list[[i]])
  buffer_alvo <- (maximo_idx + 1) * passo
  corte[[i]]  <- min(buffer_alvo, length(rank_RF))
  saturou     <- buffer_alvo > length(rank_RF)
  abline(v = corte[[i]], col = if (saturou) "red" else "black", lty = 2,
         lwd = if (saturou) 2 else 1)
  if (saturou) legend("topleft", sprintf("saturado em %d", length(rank_RF)),
                      bty = "n", text.col = "red", cex = 0.8)
  dev.off()

  t_corte[[nome]] <- round((proc.time() - t_corte_ini)[["elapsed"]])
  cat(sprintf("[TIMER] Corte %-12s: %8.1f s\n", nome, t_corte[[nome]]))
  cat(sprintf("  Ponto de corte: %d SNPs\n", corte[[i]]))
  snps_corte <- head(names(rank_RF), corte[[i]])

  genotipo_corte  <- dados[[1]][, snps_corte, drop = FALSE]
  dados[[2]]      <- cbind(genotipo_corte, fenotipo = dados[[1]]$fenotipo)
  if (!is.factor(dados[[2]]$fenotipo))
    dados[[2]]$fenotipo <- factor(dados[[2]]$fenotipo, levels = c(0, 1))

  # --- Etapa 3: GA (refinamento, fitness = F1) --------------------------------
  t_ga_ini <- proc.time()
  f <- local({
    .kernel <- kernel; .gamma <- gamma; .cost <- cost; .folds <- folds
    .dados2 <- dados[[2]]
    function(x) {
      inc <- which(x == 1)
      if (length(inc) == 0) return(0)
      dv <- .dados2[, c(inc, ncol(.dados2)), drop = FALSE]
      colnames(dv)[ncol(dv)] <- "fenotipo"
      res <- validacao_cruzada_cls(dv, folds = .folds, gamma = .gamma,
                                   cost = .cost, kernel = .kernel)
      # Fitness = F1 da classe minoritária (controles): discrimina melhor
      # quando a maioria (casos) domina a base (86/14%).
      res["f1_0_mean"]
    }
  })
  fitness <- function(x) f(x)

  set.seed(i)
  ga_result <- ga(
    type       = "binary",
    fitness    = fitness,
    nBits      = ncol(genotipo_corte),
    popSize    = popSize,
    names      = colnames(genotipo_corte),
    maxiter    = maxiter_ga,
    seed       = i,
    parallel   = TRUE,
    run        = run_ga,
    pcrossover = pcross,
    pmutation  = pmut,
    elitism    = elitism,
    suggestions = matrix(rep(1, ncol(genotipo_corte)), nrow = 1)
  )
  GA_list[[i]] <- ga_result

  sol_mat           <- as.matrix(ga_result@solution)
  sel_cols          <- which(sol_mat[1, ] == 1)
  snps_selec_ref[[i]] <- colnames(sol_mat)[sel_cols]
  t_ga[[nome]] <- round((proc.time() - t_ga_ini)[["elapsed"]])
  cat(sprintf("[TIMER] GA    %-12s: %8.1f s\n", nome, t_ga[[nome]]))

  # --- Métricas ao nível da variável -----------------------------------------
  cat(sprintf("\n--- Métricas variável | Kernel %s ---\n", nome))
  cat("  [União: corte]\n")
  m_corte <- metricas_snps(snps_corte)
  cat("  [Refinamento GA]\n")
  m_ref   <- metricas_snps(snps_selec_ref[[i]])
  metricas_var[[i]] <- list(corte = m_corte, refinamento = m_ref)

  # --- Gráfico GA ------------------------------------------------------------
  pdf(sprintf("Grafico_GA_F1_%s.pdf", nome), height = 5, width = 9)
  geracao   <- seq_len(ga_result@iter)
  df_ga     <- data.frame(
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
      labs(title  = sprintf("GA – Evolução F1 | %s | Oliveira 2015 Sim6", nome),
           x = "Geração", y = "F1 Score (fitness)") +
      theme_bw()
  )
  dev.off()

  t_kernels[[nome]] <- round((proc.time() - t_k_ini)[["elapsed"]])
  cat(sprintf("  Tempo kernel %s: %d s\n", nome, t_kernels[[nome]]))
}

# ==============================================================================
# SEÇÃO 4 – RESUMO FINAL COMPARATIVO
# ==============================================================================
cat("\n")
cat("==============================================================\n")
cat("  RESUMO FINAL – Simulação 6 de Oliveira (2015)\n")
cat("  SMS com SMOTE-NC + k-fold estratificado (k=10) + F1-Score\n")
cat("==============================================================\n")
cat(sprintf("  n total: 1000 | Casos (1): %d | Controles (0): %d\n",
            sum(dados[[1]]$fenotipo == 1), sum(dados[[1]]$fenotipo == 0)))
cat(sprintf("  SNPs causais verdadeiros: %s\n",
            paste(snps_causais_verdadeiros, collapse = ", ")))
cat("\n")

# Tabela de resultados
resultados <- data.frame(
  Kernel     = sapply(kernels_config, `[[`, "nome"),
  N_corte    = sapply(corte, identity),
  V_corte    = sapply(metricas_var, function(m) m$corte$V),
  Prec_corte = sapply(metricas_var, function(m) round(100 * m$corte$precisao, 1)),
  Sens_corte = sapply(metricas_var, function(m) round(100 * m$corte$sensibilidade, 1)),
  F1_corte   = sapply(metricas_var, function(m) round(100 * m$corte$f1_var, 1)),
  N_GA       = sapply(snps_selec_ref, length),
  V_GA       = sapply(metricas_var, function(m) m$refinamento$V),
  Prec_GA    = sapply(metricas_var, function(m) round(100 * m$refinamento$precisao, 1)),
  Sens_GA    = sapply(metricas_var, function(m) round(100 * m$refinamento$sensibilidade, 1)),
  F1_GA      = sapply(metricas_var, function(m) round(100 * m$refinamento$f1_var, 1)),
  stringsAsFactors = FALSE
)
print(resultados, row.names = FALSE)

# União e interseção das soluções GA
snps_uniao     <- Reduce(union,     snps_selec_ref)
snps_intersec  <- Reduce(intersect, snps_selec_ref)
cat("\n--- União de todos os kernels (GA) ---\n")
metricas_snps(snps_uniao)
cat("\n--- Interseção de todos os kernels (GA) ---\n")
if (length(snps_intersec) > 0) metricas_snps(snps_intersec) else cat("  (conjunto vazio)\n")

# Salvar resultados em CSV
write.csv(resultados,
          file = "Resultados_SMS_Oliveira2015_Sim6.csv",
          row.names = FALSE)
write.csv(data.frame(SNP = snps_uniao,
                     causal = snps_uniao %in% snps_causais_verdadeiros),
          file = "SNPs_uniao_GA.csv", row.names = FALSE)
write.csv(data.frame(SNP = snps_intersec,
                     causal = snps_intersec %in% snps_causais_verdadeiros),
          file = "SNPs_intersecao_GA.csv", row.names = FALSE)

t_total <- round((proc.time() - t_total_inicio)[["elapsed"]])
cat("\n===== Tempos de Execução Detalhados (segundos) =====\n")
cat(sprintf("  %-24s: %8.1f s\n", "Simulação",     t_sim))
cat(sprintf("  %-24s: %8.1f s\n", "Random Forest", t_rf))
cat("\n  -- Corte (curva F1 na trilha RF) --\n")
for (.k in sapply(kernels_config, `[[`, "nome"))
  cat(sprintf("  Corte %-18s: %8.1f s\n", .k, t_corte[[.k]]))
cat("\n  -- GA (refinamento) --\n")
for (.k in sapply(kernels_config, `[[`, "nome"))
  cat(sprintf("  GA    %-18s: %8.1f s\n", .k, t_ga[[.k]]))
cat("\n  -- Kernel total (Corte+GA) --\n")
for (.k in sapply(kernels_config, `[[`, "nome"))
  cat(sprintf("  %-24s: %8.1f s\n", .k, t_kernels[[.k]]))
cat(sprintf("  %-24s: %8.1f s\n", "TOTAL", t_total))
cat(sprintf("\nTempo total: %d s (%.1f min)\n", t_total, t_total / 60))
cat("Arquivos gerados na pasta:", getwd(), "\n")
