# Definindo a trilha de dados
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ---- Timer global ----
t_total_inicio <- proc.time()
t_kernels      <- list()    # lista para tempos por kernel
t_corte        <- list()  # tempo da etapa de corte por kernel
t_ga           <- list()  # tempo do GA por kernel
t_sim          <- 0       # tempo da simulação da base
t_rf           <- 0       # tempo da Random Forest

# Instalando e carregando pacotes
# Nome dos pacotes
packages <- c("scrime","e1071","kernlab","randomForest", 
              "doParallel", "GA", "ggplot2")

# Instalando os pacotes que ainda nao foram instalados
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

#------------------------------------------------------------
# Função auxiliar: dado um beta0, simula e mede proporção de 1s
# A classe é definida por linpred > 0 (equivalente a prob > 0.5)
#------------------------------------------------------------
.prop_ones_given_beta0 <- function(beta0,
                                   n.obs, n.snp,
                                   list.snp, list.ia,
                                   beta, maf,
                                   err.fun = rnorm,
                                   rand = 123,
                                   reps = 1) {
  # Média sobre 'reps' para reduzir ruído estocástico (opcional)
  acc <- 0
  for (i in 1:reps) {
    sim <- simulateSNPglm(
      n.obs   = n.obs,
      n.snp   = n.snp,
      list.snp= list.snp,
      list.ia = list.ia,
      beta0   = beta0,
      beta    = beta,
      maf     = maf,
      err.fun = err.fun,
      rand    = rand + i # muda a semente em cada réplica
    )
    linpred <- sim$y
    acc <- acc + mean(linpred > 0)
  }
  acc / reps
}

#------------------------------------------------------------
# Busca de beta0 para atingir uma proporção-alvo de 1s
# Usa uniroot em intervalo amplo e, se necessário, expande.
#------------------------------------------------------------
find_beta0_for_target <- function(target_prop,
                                  n.obs, n.snp,
                                  list.snp, list.ia,
                                  beta, maf,
                                  err.fun = rnorm,
                                  rand = 123,
                                  reps = 3,
                                  lower = -20, upper = 20,
                                  tol = 1e-3, max_expand = 5) {
  # Garante que target está em (0,1)
  if (target_prop <= 0 || target_prop >= 1) {
    stop("target_prop deve estar em (0,1).")
  }
  
  f <- function(b0) {
    .prop_ones_given_beta0(
      beta0 = b0,
      n.obs = n.obs, n.snp = n.snp,
      list.snp = list.snp, list.ia = list.ia,
      beta = beta, maf = maf,
      err.fun = err.fun, rand = rand,
      reps = reps
    ) - target_prop
  }
  
  L <- lower; U <- upper
  fL <- f(L); fU <- f(U)
  
  # Expande o intervalo se necessário, algumas vezes
  k <- 0
  while (sign(fL) == sign(fU) && k < max_expand) {
    L <- L - 10; U <- U + 10
    fL <- f(L);  fU <- f(U)
    k <- k + 1
  }
  
  if (sign(fL) == sign(fU)) {
    stop("Não foi possível bracketing: tente ajustar 'lower/upper' ou aumentar 'max_expand'.")
  }
  
  uniroot(f, interval = c(L, U), tol = tol)$root
}

#------------------------------------------------------------
# Função principal:
# - entra com n1 (qtde de 1s) e n0 (qtde de 0s) desejados
# - retorna base simulada com proporção próxima à desejada
# - se exact_by_rank=TRUE, força exatamente n1 e n0 por ranqueamento
#------------------------------------------------------------
simulate_scrime_binary_with_ratio <- function(n1, n0,
                                              n.snp,
                                              list.snp, list.ia,
                                              beta = rep(2, length(list.snp)),
                                              maf  = c(0.1, 0.4),
                                              err.fun = rnorm,
                                              rand = 123,
                                              reps_beta0 = 3,
                                              exact_by_rank = FALSE) {
  n.obs <- n1 + n0
  target <- n1 / n.obs
  
  # 1) Encontrar beta0 que atinge a proporção desejada (em média)
  beta0_star <- find_beta0_for_target(
    target_prop = target,
    n.obs = n.obs, n.snp = n.snp,
    list.snp = list.snp, list.ia = list.ia,
    beta = beta, maf = maf,
    err.fun = err.fun, rand = rand,
    reps = reps_beta0,
    lower = -20, upper = 20, tol = 1e-3
  )
  
  # 2) Simulação final com esse beta0
  sim <- simulateSNPglm(
    n.obs   = n.obs,
    n.snp   = n.snp,
    list.snp= list.snp,
    list.ia = list.ia,
    beta0   = beta0_star,
    beta    = beta,
    maf     = maf,
    err.fun = err.fun,
    rand    = rand
  )
  
  X <- as.data.frame(sim$x)
  linpred <- sim$y
  prob <- plogis(linpred)
  
  if (!exact_by_rank) {
    # Classificação "natural" do modelo (limiar 0 no logito)
    Y <- as.numeric(linpred > 0)
  } else {
    # Força exatamente n1 e n0 via ranqueamento do linpred
    ord <- order(linpred, decreasing = TRUE)
    Y <- integer(length(linpred))
    Y[ord[seq_len(n1)]] <- 1
  }
  
  base <- data.frame(X, linpred = linpred, prob = prob, fenotipo = Y)
  out <- list(
    data = base,
    beta0 = beta0_star,
    achieved_prop_1 = mean(Y == 1),
    counts = table(Y)
  )
  return(out)
}

#------------------------------------------------------------
# EXEMPLOS DE USO
#------------------------------------------------------------
# Parâmetros do seu caso:
num_individuos <- 1000
num_snp        <- 100
list.snp <- list(1, 2, 3, c(4,5), c(6,7,8)) # Indica quais são os marcadores causais
list.ia  <- list(1, 3, 2, c(-1,3), c(1,2,3)) # Contrução das iterações entre os 5 marcadores causais
beta     <- c(2, 2, 2, 3, 4)  # Efeitos
maf      <- c(0.1, 0.4)


# A) Proporção aproximada (natural do modelo), ex.: 60% de 1s e 40% de 0s
set.seed(123)
res_approx <- simulate_scrime_binary_with_ratio(
  n1 = 600, n0 = 400,
  n.snp = num_snp,
  list.snp = list.snp, list.ia = list.ia,
  beta = beta, maf = maf,
  rand = 123,
  reps_beta0 = 3,           # usar 3-5 p/ reduzir ruído da busca
  exact_by_rank = FALSE     # NÃO força exatamente; fica muito próximo
)
res_approx$beta0
res_approx$counts
mean(res_approx$data$fenotipo)  # proporção de 1s atingida

# B) Exatamente n1 e n0 (por ranqueamento). Útil p/ estudos controlados 50/50
set.seed(123)
t_sim_ini <- proc.time()
res_exact <- simulate_scrime_binary_with_ratio(
  n1 = 200, n0 = 800,
  n.snp = num_snp,
  list.snp = list.snp, list.ia = list.ia,
  beta = beta, maf = maf,
  rand = 123,
  reps_beta0 = 3,
  exact_by_rank = TRUE       # FORÇA exatamente 500/500
)
res_exact$beta0
res_exact$counts            # Deve mostrar 500 e 500
head(res_exact$data[, c(1:8, ncol(res_exact$data))])

# Salvando em CSV
write.csv(res_exact, file = "dados.csv", row.names = FALSE)


# Lista para os dados
dados <- list()

# Carregar a base de dados
#dados_inicial <- read.csv("dados_simulados.csv")

# Remover a coluna 'prob' e manter o fenótipo binário
dados[[1]] <- subset(res_exact[[1]], select = c(-linpred, -prob))
t_sim <- round((proc.time() - t_sim_ini)[["elapsed"]])
cat(sprintf("[TIMER] Simulação     : %8.1f s\n", t_sim))

# Conferindo as primeiras linhas
head(dados[[1]])

# Conferindo a distribuição do fenótipo
table(dados[[1]]$fenotipo)


#Codigo SMS



# --- util: converte colunas list/char/factor para numéricas atômicas ---
sanitize_predictors <- function(df, target = "fenotipo") {
  pcols <- setdiff(names(df), target)
  for (nm in pcols) {
    # Se for lista, tenta extrair 1 valor por célula
    if (is.list(df[[nm]])) {
      df[[nm]] <- vapply(df[[nm]], function(z) {
        if (length(z) == 1) {
          as.numeric(z)
        } else {
          NA_real_
        }
      }, numeric(1))
    }
    # Se for fator/char, torna numérico de forma estável
    if (is.factor(df[[nm]])) {
      df[[nm]] <- as.numeric(as.character(df[[nm]]))
    } else if (is.character(df[[nm]])) {
      suppressWarnings({
        tmp <- as.numeric(df[[nm]])
      })
      df[[nm]] <- tmp
    }
  }
  df
}

# --- prepara alvo e remove 'prob' se existir ---
prepara_cls <- function(df) {
  if ("prob" %in% names(df)) df$prob <- NULL
  if (!is.factor(df$fenotipo)) {
    df$fenotipo <- factor(df$fenotipo, levels = c(0,1))
  } else {
    df$fenotipo <- factor(df$fenotipo, levels = c(0,1))
  }
  df
}

# --- validação cruzada para classificação SVM (radial/linear) ---
validacao_cruzada_cls <- function(data, folds = 5, gamma = 0.01, cost = 1, kernel = c("radial","linear")) {
  kernel <- match.arg(kernel)
  data <- prepara_cls(data)
  data <- sanitize_predictors(data, target = "fenotipo")  # <--- SANITIZA AQUI
  
  n <- nrow(data)
  index <- 1:n
  
  set.seed(123)
  geral <- matrix(sample(index), ncol = folds, byrow = TRUE)
  
  acc  <- numeric(folds)
  prec <- numeric(folds)
  rec  <- numeric(folds)
  f1   <- numeric(folds)
  
  for (i in 1:folds) {
    test_idx  <- geral[, i]
    train_idx <- setdiff(index, test_idx)
    
    trainset <- na.omit(data[train_idx, , drop = FALSE])
    testset  <- na.omit(data[test_idx , , drop = FALSE])
    
    # Se algo voltou a virar list depois do subset/na.omit (raro), saneia de novo
    trainset <- sanitize_predictors(trainset, target = "fenotipo")
    testset  <- sanitize_predictors(testset , target = "fenotipo")
    
    # Cria matrizes numéricas diretamente (sem fórmula/model.frame)
    X_train <- as.matrix(trainset[, setdiff(names(trainset), "fenotipo"), drop = FALSE])
    y_train <- trainset$fenotipo
    
    X_test  <- as.matrix(testset[, setdiff(names(testset), "fenotipo"), drop = FALSE])
    y_true  <- factor(testset$fenotipo, levels = c("0","1"))
    
    # Treina SVM
    svm_fit <- svm(
      x = X_train, y = y_train,
      kernel = kernel,
      gamma  = gamma,
      cost   = cost,
      type   = "C-classification",
      scale  = TRUE
    )
    
    # Predição
    y_pred <- predict(svm_fit, X_test)
    y_pred <- factor(y_pred, levels = c("0","1"))
    
    # Métricas
    cm <- table(Truth = y_true, Pred = y_pred)
    TP <- ifelse("1" %in% rownames(cm) && "1" %in% colnames(cm), cm["1","1"], 0)
    TN <- ifelse("0" %in% rownames(cm) && "0" %in% colnames(cm), cm["0","0"], 0)
    FP <- ifelse("0" %in% rownames(cm) && "1" %in% colnames(cm), cm["0","1"], 0)
    FN <- ifelse("1" %in% rownames(cm) && "0" %in% colnames(cm), cm["1","0"], 0)
    
    acc[i]  <- (TP + TN) / max(1, TP + TN + FP + FN)
    prec[i] <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
    rec[i]  <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    f1[i]   <- ifelse((prec[i] + rec[i]) > 0, 2 * prec[i] * rec[i] / (prec[i] + rec[i]), 0)
  }
  
  c(
    acc_mean  = mean(acc),  acc_sd  = sd(acc),
    prec_mean = mean(prec), prec_sd = sd(prec),
    rec_mean  = mean(rec),  rec_sd  = sd(rec),
    f1_mean   = mean(f1),   f1_sd   = sd(f1)
  )
}



# Carregue seus dados (sem a coluna 'prob', com 'fenotipo' binário 0/1)
# dados <- read.csv("dados_simulados.csv"); dados$prob <- NULL

# Kernel radial
#res_radial <- validacao_cruzada_cls(dados, folds = 10, gamma = 0.01, cost = 1, kernel = "radial")

# Kernel linear
#res_linear <- validacao_cruzada_cls(dados, folds = 5, gamma = 0.01, cost = 1, kernel = "linear")

#res_radial
#res_linear



# Calcula p-valor por SNP (fenótipo binário) + Bonferroni
valor.p.binario <- function(genotipo_fenotipo) {
  # Última coluna = fenótipo (0/1 ou fator com 2 níveis)
  df <- as.data.frame(genotipo_fenotipo)
  df <- na.omit(df)
  
  y <- df[[ncol(df)]]
  # Garante binário (0/1)
  if (is.factor(y)) {
    if (length(levels(y)) != 2) stop("Fenótipo fator deve ter 2 níveis.")
    # Reordena níveis para "0","1" se possível
    if (all(levels(y) %in% c("0","1"))) y <- factor(y, levels = c("0","1"))
  } else {
    # Se vier como char -> num
    if (is.character(y)) y <- as.numeric(y)
    # Checa binariedade
    if (!all(unique(na.omit(y)) %in% c(0,1))) {
      stop("Fenótipo deve ser 0/1 (ou fator com níveis 0 e 1).")
    }
    y <- factor(y, levels = c(0,1))
  }
  
  p_raw <- numeric(ncol(df) - 1L)
  m <- length(p_raw)
  
  for (i in seq_len(ncol(df) - 1L)) {
    x <- df[[i]]
    
    # Sanitiza SNP (0/1/2 típico). Se for fator/char -> num.
    if (is.factor(x)) x <- as.numeric(as.character(x))
    if (is.character(x)) x <- suppressWarnings(as.numeric(x))
    
    # Coluna constante ou sem variância -> p = 1
    if (all(is.na(x)) || length(unique(na.omit(x))) <= 1L) {
      p_raw[i] <- 1
      next
    }
    
    # Ajusta modelo logístico simples y ~ x
    fit <- try(
      glm(y ~ x, family = binomial(), na.action = na.omit),
      silent = TRUE
    )
    
    if (inherits(fit, "try-error")) {
      # Em caso de separação completa/erro numérico
      p_raw[i] <- NA_real_
    } else {
      sm <- summary(fit)$coefficients
      # coeficiente de 'x' está na 2ª linha, 4ª coluna (Pr(>|z|))
      p_raw[i] <- if (!is.na(sm[2, 4])) sm[2, 4] else NA_real_
    }
  }
  
  p_adj <- pmin(1, m * p_raw)  # Bonferroni
  
  out <- data.frame(
    `Valor p bruto`    = p_raw,
    `Valor p ajustado` = p_adj,
    row.names = names(df)[seq_len(ncol(df) - 1L)]
  )
  return(out)
}

# Calcula p-valor
pvals <- valor.p.binario(dados[[1]])

# Ordenar do menor para o maior valor-p bruto
pvals_ordenado_bruto <- pvals[order(pvals$`Valor.p.bruto`), ]

# Ordenar do menor para o maior valor-p ajustado
pvals_ordenado_ajustado <- pvals[order(pvals$`Valor.p.ajustado`), ]

head(pvals_ordenado_bruto)
head(pvals_ordenado_ajustado)


################ INÍCIO (CLASSIFICAÇÃO) ################

library(randomForest)
library(e1071)
library(GA)
library(ggplot2)


## ---------------- Parâmetros gerais ----------------
mean_svm_RF_list <- list()   # agora guarda a média de acurácia do SVM na trilha RF
GA                <- list()
minimo            <- list()  # (mantido para compatibilidade, mas usaremos máximo de acurácia)
corte             <- list()
snps_selec_corte  <- list()
snps_selec_ref    <- list()
percentual_snps   <- 0.95


########## KERNEL LINEAR #########
i <- 1                      # contador (aqui só 1 trilha)


## ---------------- Random Forest (CLASSIFICAÇÃO) ----------------
ntree <- 4000

data_temp <- as.data.frame(dados[[1]])
# garante fenotipo fator 0/1
if (!is.factor(data_temp$fenotipo)) {
  data_temp$fenotipo <- factor(data_temp$fenotipo, levels = c(0,1))
}

set.seed(1)
t_rf_ini <- proc.time()
RF <- randomForest(
  fenotipo ~ ., data = data_temp,
  ntree = ntree,
  mtry = ncol(data_temp) - 1,
  importance = TRUE
)

# Ranking por importância: MeanDecreaseGini (poderia usar "MeanDecreaseAccuracy")
imp <- importance(RF)
mdg_col <- if ("MeanDecreaseGini" %in% colnames(imp)) "MeanDecreaseGini" else colnames(imp)[ncol(imp)]
rank_RF <- sort(imp[, mdg_col], decreasing = TRUE)
t_rf <- round((proc.time() - t_rf_ini)[["elapsed"]])
cat(sprintf("[TIMER] Random Forest : %8.1f s\n", t_rf))

# ---- Timer kernel: Linear ----
t_k_Linear <- proc.time()

## ---------------- SVM (classificação) e k-folds na trilha do rank RF ----------------
gamma  <- 0.01
cost   <- 1.0
folds  <- 10
kernel <- "linear"

mean_svm_RF_list[[i]] <- numeric()
passo   <- 10
limite  <- floor((length(names(dados[[1]][, -ncol(dados[[1]])])) / passo) * percentual_snps)

t_corte_ini <- proc.time()  # inicio corte [Linear]
for (cont in 1:(limite + 1)) {
  cat("Grupo:", cont, "\n")
  j <- cont * passo
  
  # Ajuste quando j == 1
  if (j == 1) {
    var_sel <- names(rank_RF)[1]
  } else {
    var_sel <- names(rank_RF)[1:j]
  }
  
  # Base somente com as variáveis selecionadas + fenótipo
  base_cv <- dados[[1]][c(var_sel, "fenotipo")]
  
  # Avalia SVM classificação em k-folds -> pega acurácia média (primeiro elemento)
  svm_cv <- validacao_cruzada_cls(
    data   = base_cv,
    folds  = folds,
    gamma  = gamma,
    cost   = cost,
    kernel = kernel
  )
  
  acc_mean <- svm_cv[1]  # acc_mean
  mean_svm_RF_list[[i]][cont] <- acc_mean
}

## ---------------- Gráfico (Acurácia do SVM na trilha da RF) ----------------
pdf(file = "Grafico_ACCURACY_SVM_linear.pdf", height = 5, width = 9)
plot(seq(passo, limite * passo + 10, by = passo),
     mean_svm_RF_list[[i]],
     type = "o", lwd = 2,
     xlab = "Grupo de Marcadores",
     ylab = "Acurácia (SVM Linear, 10-fold)")
# Escolha do ponto de corte: agora é o MÁXIMO de acurácia
maximo_idx <- which.max(mean_svm_RF_list[[i]])
corte[[i]] <- (maximo_idx + 1) * passo
abline(v = corte[[i]], col = "black", lty = 2)
dev.off()

snps_selec_corte <- names(rank_RF[1:corte[[i]]])
snps_selec_corte

t_corte[["Linear"]] <- round((proc.time() - t_corte_ini)[["elapsed"]])
cat(sprintf("[TIMER] Corte Linear      : %8.1f s\n", t_corte[["Linear"]]))
## ---------------- Seleciona genótipo e monta dados[[2]] após corte ----------------
genotipo <- list()
genotipo[[2]] <- dados[[1]][, names(rank_RF[1:corte[[i]]]), drop = FALSE]

dados[[2]] <- cbind(genotipo[[2]], fenotipo = dados[[1]]$fenotipo)
if (!is.factor(dados[[2]]$fenotipo)) dados[[2]]$fenotipo <- factor(dados[[2]]$fenotipo, levels = c(0,1))

## ---------------- AG para REFINAMENTO (fitness = acurácia do SVM) ----------------
f <- function(x) {
  inc <- which(x == 1)
  # pelo menos 1 variável
  if (length(inc) == 0) return(0)
  
  dados_validacao <- dados[[2]][, c(inc, ncol(dados[[2]])), drop = FALSE]
  
  # se só uma coluna preditora, garante nomeações corretas
  colnames(dados_validacao)[ncol(dados_validacao)] <- "fenotipo"
  
  res <- validacao_cruzada_cls(
    data   = dados_validacao,
    folds  = folds,
    gamma  = gamma,
    cost   = cost,
    kernel = kernel
  )
  acc_mean <- res[1]  # acc_mean
  return(acc_mean)     # GA maximiza por padrão
}
fitness <- function(x) f(x)

# Parâmetros do GA
run     <- 30
maxiter <- 10
pcross  <- 0.8
pmut    <- 0.1
elitism <- 5
popSize <- 100

t_ga_ini <- proc.time()  # inicio GA [Linear]
set.seed(i)
GA[[i]] <- ga(
  type     = "binary",
  fitness  = fitness,
  nBits    = ncol(genotipo[[2]]),
  popSize  = popSize,
  names    = colnames(genotipo[[2]]),
  maxiter  = maxiter,
  seed     = i,
  parallel = TRUE,
  run      = run,
  pcrossover = pcross,
  pmutation  = pmut,
  elitism    = elitism,
  suggestions = matrix(rep(1, ncol(genotipo[[2]])),
                       ncol = ncol(genotipo[[2]]))
)

sol_mat   <- as.matrix(GA[[i]]@solution)      # garante formato matriz
best_row  <- 1                                 # escolha a 1ª solução
sel_cols  <- which(sol_mat[best_row, ] == 1)   # índices de colunas
snps_selec_ref[[i]] <- colnames(sol_mat)[sel_cols]
t_ga[["Linear"]] <- round((proc.time() - t_ga_ini)[["elapsed"]])
cat(sprintf("[TIMER] GA    Linear      : %8.1f s\n", t_ga[["Linear"]]))


snps_selec_SVM_linear_GA_acc <- snps_selec_ref[[i]]
t_kernels[["Linear"]] <- round((proc.time() - t_k_Linear)[["elapsed"]])
snps_selec_SVM_linear_GA_acc

# Relatório e gráfico do GA
summary(GA[[i]])
plot(GA[[i]])

# Salva gráfico do GA com ggplot2
pdf(file = "Grafico_GA_SVM_linear_ACC.pdf", height = 5, width = 9)
geracao        <- seq_len(GA[[i]]@iter)
mean_fitness   <- GA[[i]]@summary[, 2]
median_fitness <- GA[[i]]@summary[, 4]
best_fitness   <- GA[[i]]@summary[, 1]
Estatisticas   <- c(rep("Mediana", length(geracao)),
                    rep("Média",   length(geracao)),
                    rep("Melhor",  length(geracao)))
data_grafico_1 <- data.frame(
  Geracao     = rep(geracao, 3),
  Aptidao     = c(mean_fitness, median_fitness, best_fitness),
  Estatisticas = Estatisticas
)

ggplot(data_grafico_1, aes(x = Geracao, y = Aptidao, group = Estatisticas)) +
  geom_line(aes(colour = Estatisticas, linetype = Estatisticas), size = 2) +
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Geracao),
                                  max(data_grafico_1$Geracao), by = 1))
dev.off()



########## KERNEL RADIAL GAMMA = 0.001 #########
i <- 2                      # contador (aqui só 1 trilha)


# ---- Timer kernel: Radial_0.001 ----
t_k_Radial_0_001 <- proc.time()

## ---------------- SVM (classificação) e k-folds na trilha do rank RF ----------------
gamma  <- 0.001
cost   <- 1.0
folds  <- 10
kernel <- "radial"

mean_svm_RF_list[[i]] <- numeric()
passo   <- 10
limite  <- floor((length(names(dados[[1]][, -ncol(dados[[1]])])) / passo) * percentual_snps)

t_corte_ini <- proc.time()  # inicio corte [Radial_0.001]
for (cont in 1:(limite + 1)) {
  cat("Grupo:", cont, "\n")
  j <- cont * passo
  
  # Ajuste quando j == 1
  if (j == 1) {
    var_sel <- names(rank_RF)[1]
  } else {
    var_sel <- names(rank_RF)[1:j]
  }
  
  # Base somente com as variáveis selecionadas + fenótipo
  base_cv <- dados[[1]][c(var_sel, "fenotipo")]
  
  # Avalia SVM classificação em k-folds -> pega acurácia média (primeiro elemento)
  svm_cv <- validacao_cruzada_cls(
    data   = base_cv,
    folds  = folds,
    gamma  = gamma,
    cost   = cost,
    kernel = kernel
  )
  
  acc_mean <- svm_cv[1]  # acc_mean
  mean_svm_RF_list[[i]][cont] <- acc_mean
}

## ---------------- Gráfico (Acurácia do SVM na trilha da RF) ----------------
pdf(file = "Grafico_ACCURACY_SVM_radial_0001.pdf", height = 5, width = 9)
plot(seq(passo, limite * passo + 10, by = passo),
     mean_svm_RF_list[[i]],
     type = "o", lwd = 2,
     xlab = "Grupo de Marcadores",
     ylab = "Acurácia (SVM Radial 0.001, k-fold)")
# Escolha do ponto de corte: agora é o MÁXIMO de acurácia
maximo_idx <- which.max(mean_svm_RF_list[[i]])
corte[[i]] <- (maximo_idx + 1) * passo
abline(v = corte[[i]], col = "black", lty = 2)
dev.off()

snps_selec_corte <- names(rank_RF[1:corte[[i]]])
snps_selec_corte

t_corte[["Radial_0.001"]] <- round((proc.time() - t_corte_ini)[["elapsed"]])
cat(sprintf("[TIMER] Corte Radial_0.001: %8.1f s\n", t_corte[["Radial_0.001"]]))
## ---------------- Seleciona genótipo e monta dados[[2]] após corte ----------------
genotipo <- list()
genotipo[[2]] <- dados[[1]][, names(rank_RF[1:corte[[i]]]), drop = FALSE]

dados[[2]] <- cbind(genotipo[[2]], fenotipo = dados[[1]]$fenotipo)
if (!is.factor(dados[[2]]$fenotipo)) dados[[2]]$fenotipo <- factor(dados[[2]]$fenotipo, levels = c(0,1))

## ---------------- AG para REFINAMENTO (fitness = acurácia do SVM) ----------------
f <- function(x) {
  inc <- which(x == 1)
  # pelo menos 1 variável
  if (length(inc) == 0) return(0)
  
  dados_validacao <- dados[[2]][, c(inc, ncol(dados[[2]])), drop = FALSE]
  
  # se só uma coluna preditora, garante nomeações corretas
  colnames(dados_validacao)[ncol(dados_validacao)] <- "fenotipo"
  
  res <- validacao_cruzada_cls(
    data   = dados_validacao,
    folds  = folds,
    gamma  = gamma,
    cost   = cost,
    kernel = kernel
  )
  acc_mean <- res[1]  # acc_mean
  return(acc_mean)     # GA maximiza por padrão
}
fitness <- function(x) f(x)

# Parâmetros do GA
run     <- 30
maxiter <- 10
pcross  <- 0.8
pmut    <- 0.1
elitism <- 5
popSize <- 100

t_ga_ini <- proc.time()  # inicio GA [Radial_0.001]
set.seed(i)
GA[[i]] <- ga(
  type     = "binary",
  fitness  = fitness,
  nBits    = ncol(genotipo[[2]]),
  popSize  = popSize,
  names    = colnames(genotipo[[2]]),
  maxiter  = maxiter,
  seed     = i,
  parallel = TRUE,
  run      = run,
  pcrossover = pcross,
  pmutation  = pmut,
  elitism    = elitism,
  suggestions = matrix(rep(1, ncol(genotipo[[2]])),
                       ncol = ncol(genotipo[[2]]))
)

sol_mat   <- as.matrix(GA[[i]]@solution)      # garante formato matriz
best_row  <- 1                                 # escolha a 1ª solução
sel_cols  <- which(sol_mat[best_row, ] == 1)   # índices de colunas
snps_selec_ref[[i]] <- colnames(sol_mat)[sel_cols]
t_ga[["Radial_0.001"]] <- round((proc.time() - t_ga_ini)[["elapsed"]])
cat(sprintf("[TIMER] GA    Radial_0.001: %8.1f s\n", t_ga[["Radial_0.001"]]))


snps_selec_SVM_radial_0001_GA_acc <- snps_selec_ref[[i]]
t_kernels[["Radial_0.001"]] <- round((proc.time() - t_k_Radial_0_001)[["elapsed"]])
snps_selec_SVM_radial_0001_GA_acc

# Relatório e gráfico do GA
summary(GA[[i]])
plot(GA[[i]])

# Salva gráfico do GA com ggplot2
pdf(file = "Grafico_GA_SVM_radial_0001_ACC.pdf", height = 5, width = 9)
geracao        <- seq_len(GA[[i]]@iter)
mean_fitness   <- GA[[i]]@summary[, 2]
median_fitness <- GA[[i]]@summary[, 4]
best_fitness   <- GA[[i]]@summary[, 1]
Estatisticas   <- c(rep("Mediana", length(geracao)),
                    rep("Média",   length(geracao)),
                    rep("Melhor",  length(geracao)))
data_grafico_1 <- data.frame(
  Geracao     = rep(geracao, 3),
  Aptidao     = c(mean_fitness, median_fitness, best_fitness),
  Estatisticas = Estatisticas
)

ggplot(data_grafico_1, aes(x = Geracao, y = Aptidao, group = Estatisticas)) +
  geom_line(aes(colour = Estatisticas, linetype = Estatisticas), size = 2) +
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Geracao),
                                  max(data_grafico_1$Geracao), by = 1))
dev.off()


########## KERNEL RADIAL GAMMA = 0.01 #########
i <- 3                      # contador (aqui só 1 trilha)


# ---- Timer kernel: Radial_0.01 ----
t_k_Radial_0_01 <- proc.time()

## ---------------- SVM (classificação) e k-folds na trilha do rank RF ----------------
gamma  <- 0.01
cost   <- 1.0
folds  <- 10
kernel <- "radial"

mean_svm_RF_list[[i]] <- numeric()
passo   <- 10
limite  <- floor((length(names(dados[[1]][, -ncol(dados[[1]])])) / passo) * percentual_snps)

t_corte_ini <- proc.time()  # inicio corte [Radial_0.01]
for (cont in 1:(limite + 1)) {
  cat("Grupo:", cont, "\n")
  j <- cont * passo
  
  # Ajuste quando j == 1
  if (j == 1) {
    var_sel <- names(rank_RF)[1]
  } else {
    var_sel <- names(rank_RF)[1:j]
  }
  
  # Base somente com as variáveis selecionadas + fenótipo
  base_cv <- dados[[1]][c(var_sel, "fenotipo")]
  
  # Avalia SVM classificação em k-folds -> pega acurácia média (primeiro elemento)
  svm_cv <- validacao_cruzada_cls(
    data   = base_cv,
    folds  = folds,
    gamma  = gamma,
    cost   = cost,
    kernel = kernel
  )
  
  acc_mean <- svm_cv[1]  # acc_mean
  mean_svm_RF_list[[i]][cont] <- acc_mean
}

## ---------------- Gráfico (Acurácia do SVM na trilha da RF) ----------------
pdf(file = "Grafico_ACCURACY_SVM_radial_001.pdf", height = 5, width = 9)
plot(seq(passo, limite * passo + 10, by = passo),
     mean_svm_RF_list[[i]],
     type = "o", lwd = 2,
     xlab = "Grupo de Marcadores",
     ylab = "Acurácia (SVM Radial 0.01, 10-fold)")
# Escolha do ponto de corte: agora é o MÁXIMO de acurácia
maximo_idx <- which.max(mean_svm_RF_list[[i]])
corte[[i]] <- (maximo_idx + 1) * passo
abline(v = corte[[i]], col = "black", lty = 2)
dev.off()

snps_selec_corte <- names(rank_RF[1:corte[[i]]])
snps_selec_corte

t_corte[["Radial_0.01"]] <- round((proc.time() - t_corte_ini)[["elapsed"]])
cat(sprintf("[TIMER] Corte Radial_0.01 : %8.1f s\n", t_corte[["Radial_0.01"]]))
## ---------------- Seleciona genótipo e monta dados[[2]] após corte ----------------
genotipo <- list()
genotipo[[2]] <- dados[[1]][, names(rank_RF[1:corte[[i]]]), drop = FALSE]

dados[[2]] <- cbind(genotipo[[2]], fenotipo = dados[[1]]$fenotipo)
if (!is.factor(dados[[2]]$fenotipo)) dados[[2]]$fenotipo <- factor(dados[[2]]$fenotipo, levels = c(0,1))

## ---------------- AG para REFINAMENTO (fitness = acurácia do SVM) ----------------
f <- function(x) {
  inc <- which(x == 1)
  # pelo menos 1 variável
  if (length(inc) == 0) return(0)
  
  dados_validacao <- dados[[2]][, c(inc, ncol(dados[[2]])), drop = FALSE]
  
  # se só uma coluna preditora, garante nomeações corretas
  colnames(dados_validacao)[ncol(dados_validacao)] <- "fenotipo"
  
  res <- validacao_cruzada_cls(
    data   = dados_validacao,
    folds  = folds,
    gamma  = gamma,
    cost   = cost,
    kernel = kernel
  )
  acc_mean <- res[1]  # acc_mean
  return(acc_mean)     # GA maximiza por padrão
}
fitness <- function(x) f(x)

# Parâmetros do GA
run     <- 30
maxiter <- 10
pcross  <- 0.8
pmut    <- 0.1
elitism <- 5
popSize <- 100

t_ga_ini <- proc.time()  # inicio GA [Radial_0.01]
set.seed(i)
GA[[i]] <- ga(
  type     = "binary",
  fitness  = fitness,
  nBits    = ncol(genotipo[[2]]),
  popSize  = popSize,
  names    = colnames(genotipo[[2]]),
  maxiter  = maxiter,
  seed     = i,
  parallel = TRUE,
  run      = run,
  pcrossover = pcross,
  pmutation  = pmut,
  elitism    = elitism,
  suggestions = matrix(rep(1, ncol(genotipo[[2]])),
                       ncol = ncol(genotipo[[2]]))
)

sol_mat   <- as.matrix(GA[[i]]@solution)      # garante formato matriz
best_row  <- 1                                 # escolha a 1ª solução
sel_cols  <- which(sol_mat[best_row, ] == 1)   # índices de colunas
snps_selec_ref[[i]] <- colnames(sol_mat)[sel_cols]
t_ga[["Radial_0.01"]] <- round((proc.time() - t_ga_ini)[["elapsed"]])
cat(sprintf("[TIMER] GA    Radial_0.01 : %8.1f s\n", t_ga[["Radial_0.01"]]))


snps_selec_SVM_radial_001_GA_acc <- snps_selec_ref[[i]]
t_kernels[["Radial_0.01"]] <- round((proc.time() - t_k_Radial_0_01)[["elapsed"]])
snps_selec_SVM_radial_001_GA_acc

# Relatório e gráfico do GA
summary(GA[[i]])
plot(GA[[i]])

# Salva gráfico do GA com ggplot2
pdf(file = "Grafico_GA_SVM_radial_001_ACC.pdf", height = 5, width = 9)
geracao        <- seq_len(GA[[i]]@iter)
mean_fitness   <- GA[[i]]@summary[, 2]
median_fitness <- GA[[i]]@summary[, 4]
best_fitness   <- GA[[i]]@summary[, 1]
Estatisticas   <- c(rep("Mediana", length(geracao)),
                    rep("Média",   length(geracao)),
                    rep("Melhor",  length(geracao)))
data_grafico_1 <- data.frame(
  Geracao     = rep(geracao, 3),
  Aptidao     = c(mean_fitness, median_fitness, best_fitness),
  Estatisticas = Estatisticas
)

ggplot(data_grafico_1, aes(x = Geracao, y = Aptidao, group = Estatisticas)) +
  geom_line(aes(colour = Estatisticas, linetype = Estatisticas), size = 2) +
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Geracao),
                                  max(data_grafico_1$Geracao), by = 1))
dev.off()


########## KERNEL RADIAL GAMMA = 0.1 #########
i <- 4                      # contador (aqui só 1 trilha)


# ---- Timer kernel: Radial_0.1 ----
t_k_Radial_0_1 <- proc.time()

## ---------------- SVM (classificação) e k-folds na trilha do rank RF ----------------
gamma  <- 0.1
cost   <- 1.0
folds  <- 10
kernel <- "radial"

mean_svm_RF_list[[i]] <- numeric()
passo   <- 10
limite  <- floor((length(names(dados[[1]][, -ncol(dados[[1]])])) / passo) * percentual_snps)

t_corte_ini <- proc.time()  # inicio corte [Radial_0.1]
for (cont in 1:(limite + 1)) {
  cat("Grupo:", cont, "\n")
  j <- cont * passo
  
  # Ajuste quando j == 1
  if (j == 1) {
    var_sel <- names(rank_RF)[1]
  } else {
    var_sel <- names(rank_RF)[1:j]
  }
  
  # Base somente com as variáveis selecionadas + fenótipo
  base_cv <- dados[[1]][c(var_sel, "fenotipo")]
  
  # Avalia SVM classificação em k-folds -> pega acurácia média (primeiro elemento)
  svm_cv <- validacao_cruzada_cls(
    data   = base_cv,
    folds  = folds,
    gamma  = gamma,
    cost   = cost,
    kernel = kernel
  )
  
  acc_mean <- svm_cv[1]  # acc_mean
  mean_svm_RF_list[[i]][cont] <- acc_mean
}

## ---------------- Gráfico (Acurácia do SVM na trilha da RF) ----------------
pdf(file = "Grafico_ACCURACY_SVM_radial_01.pdf", height = 5, width = 9)
plot(seq(passo, limite * passo + 10, by = passo),
     mean_svm_RF_list[[i]],
     type = "o", lwd = 2,
     xlab = "Grupo de Marcadores",
     ylab = "Acurácia (SVM Radial 0.1, 10-fold)")
# Escolha do ponto de corte: agora é o MÁXIMO de acurácia
maximo_idx <- which.max(mean_svm_RF_list[[i]])
corte[[i]] <- (maximo_idx + 1) * passo
abline(v = corte[[i]], col = "black", lty = 2)
dev.off()

snps_selec_corte <- names(rank_RF[1:corte[[i]]])
snps_selec_corte

t_corte[["Radial_0.1"]] <- round((proc.time() - t_corte_ini)[["elapsed"]])
cat(sprintf("[TIMER] Corte Radial_0.1  : %8.1f s\n", t_corte[["Radial_0.1"]]))
## ---------------- Seleciona genótipo e monta dados[[2]] após corte ----------------
genotipo <- list()
genotipo[[2]] <- dados[[1]][, names(rank_RF[1:corte[[i]]]), drop = FALSE]

dados[[2]] <- cbind(genotipo[[2]], fenotipo = dados[[1]]$fenotipo)
if (!is.factor(dados[[2]]$fenotipo)) dados[[2]]$fenotipo <- factor(dados[[2]]$fenotipo, levels = c(0,1))

## ---------------- AG para REFINAMENTO (fitness = acurácia do SVM) ----------------
f <- function(x) {
  inc <- which(x == 1)
  # pelo menos 1 variável
  if (length(inc) == 0) return(0)
  
  dados_validacao <- dados[[2]][, c(inc, ncol(dados[[2]])), drop = FALSE]
  
  # se só uma coluna preditora, garante nomeações corretas
  colnames(dados_validacao)[ncol(dados_validacao)] <- "fenotipo"
  
  res <- validacao_cruzada_cls(
    data   = dados_validacao,
    folds  = folds,
    gamma  = gamma,
    cost   = cost,
    kernel = kernel
  )
  acc_mean <- res[1]  # acc_mean
  return(acc_mean)     # GA maximiza por padrão
}
fitness <- function(x) f(x)

# Parâmetros do GA
run     <- 30
maxiter <- 10
pcross  <- 0.8
pmut    <- 0.1
elitism <- 5
popSize <- 100

t_ga_ini <- proc.time()  # inicio GA [Radial_0.1]
set.seed(i)
GA[[i]] <- ga(
  type     = "binary",
  fitness  = fitness,
  nBits    = ncol(genotipo[[2]]),
  popSize  = popSize,
  names    = colnames(genotipo[[2]]),
  maxiter  = maxiter,
  seed     = i,
  parallel = TRUE,
  run      = run,
  pcrossover = pcross,
  pmutation  = pmut,
  elitism    = elitism,
  suggestions = matrix(rep(1, ncol(genotipo[[2]])),
                       ncol = ncol(genotipo[[2]]))
)

sol_mat   <- as.matrix(GA[[i]]@solution)      # garante formato matriz
best_row  <- 1                                 # escolha a 1ª solução
sel_cols  <- which(sol_mat[best_row, ] == 1)   # índices de colunas
snps_selec_ref[[i]] <- colnames(sol_mat)[sel_cols]
t_ga[["Radial_0.1"]] <- round((proc.time() - t_ga_ini)[["elapsed"]])
cat(sprintf("[TIMER] GA    Radial_0.1  : %8.1f s\n", t_ga[["Radial_0.1"]]))


snps_selec_SVM_radial_01_GA_acc <- snps_selec_ref[[i]]
t_kernels[["Radial_0.1"]] <- round((proc.time() - t_k_Radial_0_1)[["elapsed"]])
snps_selec_SVM_radial_01_GA_acc

# Relatório e gráfico do GA
summary(GA[[i]])
plot(GA[[i]])

# Salva gráfico do GA com ggplot2
pdf(file = "Grafico_GA_SVM_radial_01_ACC.pdf", height = 5, width = 9)
geracao        <- seq_len(GA[[i]]@iter)
mean_fitness   <- GA[[i]]@summary[, 2]
median_fitness <- GA[[i]]@summary[, 4]
best_fitness   <- GA[[i]]@summary[, 1]
Estatisticas   <- c(rep("Mediana", length(geracao)),
                    rep("Média",   length(geracao)),
                    rep("Melhor",  length(geracao)))
data_grafico_1 <- data.frame(
  Geracao     = rep(geracao, 3),
  Aptidao     = c(mean_fitness, median_fitness, best_fitness),
  Estatisticas = Estatisticas
)

ggplot(data_grafico_1, aes(x = Geracao, y = Aptidao, group = Estatisticas)) +
  geom_line(aes(colour = Estatisticas, linetype = Estatisticas), size = 2) +
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Geracao),
                                  max(data_grafico_1$Geracao), by = 1))
dev.off()


########## KERNEL RADIAL GAMMA = 1 #########
i <- 5                      # contador (aqui só 1 trilha)


# ---- Timer kernel: Radial_1.0 ----
t_k_Radial_1_0 <- proc.time()

## ---------------- SVM (classificação) e k-folds na trilha do rank RF ----------------
gamma  <- 1.0
cost   <- 1.0
folds  <- 10
kernel <- "radial"

mean_svm_RF_list[[i]] <- numeric()
passo   <- 10
limite  <- floor((length(names(dados[[1]][, -ncol(dados[[1]])])) / passo) * percentual_snps)

t_corte_ini <- proc.time()  # inicio corte [Radial_1.0]
for (cont in 1:(limite + 1)) {
  cat("Grupo:", cont, "\n")
  j <- cont * passo
  
  # Ajuste quando j == 1
  if (j == 1) {
    var_sel <- names(rank_RF)[1]
  } else {
    var_sel <- names(rank_RF)[1:j]
  }
  
  # Base somente com as variáveis selecionadas + fenótipo
  base_cv <- dados[[1]][c(var_sel, "fenotipo")]
  
  # Avalia SVM classificação em k-folds -> pega acurácia média (primeiro elemento)
  svm_cv <- validacao_cruzada_cls(
    data   = base_cv,
    folds  = folds,
    gamma  = gamma,
    cost   = cost,
    kernel = kernel
  )
  
  acc_mean <- svm_cv[1]  # acc_mean
  mean_svm_RF_list[[i]][cont] <- acc_mean
}

## ---------------- Gráfico (Acurácia do SVM na trilha da RF) ----------------
pdf(file = "Grafico_ACCURACY_SVM_radial_1.pdf", height = 5, width = 9)
plot(seq(passo, limite * passo + 10, by = passo),
     mean_svm_RF_list[[i]],
     type = "o", lwd = 2,
     xlab = "Grupo de Marcadores",
     ylab = "Acurácia (SVM Radial 1.0, 10-fold)")
# Escolha do ponto de corte: agora é o MÁXIMO de acurácia
maximo_idx <- which.max(mean_svm_RF_list[[i]])
corte[[i]] <- (maximo_idx + 1) * passo
abline(v = corte[[i]], col = "black", lty = 2)
dev.off()

snps_selec_corte <- names(rank_RF[1:corte[[i]]])
snps_selec_corte

t_corte[["Radial_1.0"]] <- round((proc.time() - t_corte_ini)[["elapsed"]])
cat(sprintf("[TIMER] Corte Radial_1.0  : %8.1f s\n", t_corte[["Radial_1.0"]]))
## ---------------- Seleciona genótipo e monta dados[[2]] após corte ----------------
genotipo <- list()
genotipo[[2]] <- dados[[1]][, names(rank_RF[1:corte[[i]]]), drop = FALSE]

dados[[2]] <- cbind(genotipo[[2]], fenotipo = dados[[1]]$fenotipo)
if (!is.factor(dados[[2]]$fenotipo)) dados[[2]]$fenotipo <- factor(dados[[2]]$fenotipo, levels = c(0,1))

## ---------------- AG para REFINAMENTO (fitness = acurácia do SVM) ----------------
f <- function(x) {
  inc <- which(x == 1)
  # pelo menos 1 variável
  if (length(inc) == 0) return(0)
  
  dados_validacao <- dados[[2]][, c(inc, ncol(dados[[2]])), drop = FALSE]
  
  # se só uma coluna preditora, garante nomeações corretas
  colnames(dados_validacao)[ncol(dados_validacao)] <- "fenotipo"
  
  res <- validacao_cruzada_cls(
    data   = dados_validacao,
    folds  = folds,
    gamma  = gamma,
    cost   = cost,
    kernel = kernel
  )
  acc_mean <- res[1]  # acc_mean
  return(acc_mean)     # GA maximiza por padrão
}
fitness <- function(x) f(x)

# Parâmetros do GA
run     <- 30
maxiter <- 10
pcross  <- 0.8
pmut    <- 0.1
elitism <- 5
popSize <- 100

t_ga_ini <- proc.time()  # inicio GA [Radial_1.0]
set.seed(i)
GA[[i]] <- ga(
  type     = "binary",
  fitness  = fitness,
  nBits    = ncol(genotipo[[2]]),
  popSize  = popSize,
  names    = colnames(genotipo[[2]]),
  maxiter  = maxiter,
  seed     = i,
  parallel = TRUE,
  run      = run,
  pcrossover = pcross,
  pmutation  = pmut,
  elitism    = elitism,
  suggestions = matrix(rep(1, ncol(genotipo[[2]])),
                       ncol = ncol(genotipo[[2]]))
)

sol_mat   <- as.matrix(GA[[i]]@solution)      # garante formato matriz
best_row  <- 1                                 # escolha a 1ª solução
sel_cols  <- which(sol_mat[best_row, ] == 1)   # índices de colunas
snps_selec_ref[[i]] <- colnames(sol_mat)[sel_cols]
t_ga[["Radial_1.0"]] <- round((proc.time() - t_ga_ini)[["elapsed"]])
cat(sprintf("[TIMER] GA    Radial_1.0  : %8.1f s\n", t_ga[["Radial_1.0"]]))


snps_selec_SVM_radial_1_GA_acc <- snps_selec_ref[[i]]
t_kernels[["Radial_1.0"]] <- round((proc.time() - t_k_Radial_1_0)[["elapsed"]])
snps_selec_SVM_radial_1_GA_acc

# Relatório e gráfico do GA
summary(GA[[i]])
plot(GA[[i]])

# Salva gráfico do GA com ggplot2
pdf(file = "Grafico_GA_SVM_radial_1_ACC.pdf", height = 5, width = 9)
geracao        <- seq_len(GA[[i]]@iter)
mean_fitness   <- GA[[i]]@summary[, 2]
median_fitness <- GA[[i]]@summary[, 4]
best_fitness   <- GA[[i]]@summary[, 1]
Estatisticas   <- c(rep("Mediana", length(geracao)),
                    rep("Média",   length(geracao)),
                    rep("Melhor",  length(geracao)))
data_grafico_1 <- data.frame(
  Geracao     = rep(geracao, 3),
  Aptidao     = c(mean_fitness, median_fitness, best_fitness),
  Estatisticas = Estatisticas
)

ggplot(data_grafico_1, aes(x = Geracao, y = Aptidao, group = Estatisticas)) +
  geom_line(aes(colour = Estatisticas, linetype = Estatisticas), size = 2) +
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Geracao),
                                  max(data_grafico_1$Geracao), by = 1))
dev.off()

############ Termino do SMS #############


####Uniao dos SNPs selecionados####
#O argumento lista_SNPs tem que estar no formato lista com pelo menos 2 elementos.
uniao_snps<-function(lista_SNPs){
  uniao <- list()
  uniao[[1]]<-union(lista_SNPs[[1]],lista_SNPs[[2]])
  for (a in 1:(length(lista_SNPs)-2)){
    uniao[[a+1]]<-union(uniao[[a]],lista_SNPs[[a+2]])
  }
  uniao_final <- uniao[[a+1]]
  return(uniao_final)
}

uniao_snps(snps_selec_ref)

uniao_final<-uniao_snps(snps_selec_ref)

####Intersecao dos SNPs selecionados####
#O argumento lista_SNPs tem que estar no formato lista com pelo menos 2 elementos.
intersecao_snps<-function(lista_SNPs){
  intersecao <- list()
  intersecao[[1]]<-intersect(lista_SNPs[[1]],lista_SNPs[[2]])
  for (a in 1:(length(lista_SNPs)-2)){
    intersecao[[a+1]]<-intersect(intersecao[[a]],lista_SNPs[[a+2]])
  }
  intersecao_final <- intersecao[[a+1]]
  if (length(intersecao_final) == 0L) {
    return("Conjunto vazio")
  } else
    return(intersecao_final)
}

intersecao_snps(snps_selec_ref)

intersecao_final<-intersecao_snps(snps_selec_ref)

###############################################
## ACRESCIMOS: Padronização de saídas (CLASSIFICAÇÃO)
## NÃO altera nada do SMS original; apenas acrescenta
###############################################

## 1) Preparos: importância da RF como data.frame + valor-p binário
# 'imp' e 'mdg_col' já foram definidos no seu SMS de classificação
rank_RF_df <- data.frame(Importancia_RF = imp[, mdg_col, drop = TRUE])
# Salva o rank completo (equivalente ao "rank_global" do Bessel)
write.csv(rank_RF_df, "rank_Random_Forest_CLASSIF.txt", row.names = TRUE)

# Recalcula/garante p-valores binários com a função original
# (se você já executou 'pvals <- valor.p.binario(dados[[1]])' acima, isto reaproveita)
if (!exists("pvals")) {
  pvals <- valor.p.binario(dados[[1]])
}

# Seleções por valor-p (bruto e Bonferroni) para imprimir no relatório
pvals_bruto_selec    <- rownames(pvals)[which(pvals[["Valor.p.bruto"]]    <= 0.05)]
pvals_ajustado_selec <- rownames(pvals)[which(pvals[["Valor.p.ajustado"]] <= 0.05)]

## 2) Consolida as listas de SNPs selecionados por kernel (já produzidas no seu SMS)
# Nomes conforme você já criou nas trilhas:
#   snps_selec_SVM_linear_GA_acc
#   snps_selec_SVM_radial_0001_GA_acc
#   snps_selec_SVM_radial_001_GA_acc
#   snps_selec_SVM_radial_01_GA_acc
#   snps_selec_SVM_radial_1_GA_acc

selec_por_kernel <- list(
  "Linear"    = if (exists("snps_selec_SVM_linear_GA_acc"))        snps_selec_SVM_linear_GA_acc        else character(0),
  "Radial_0.001" = if (exists("snps_selec_SVM_radial_0001_GA_acc")) snps_selec_SVM_radial_0001_GA_acc else character(0),
  "Radial_0.01"  = if (exists("snps_selec_SVM_radial_001_GA_acc"))  snps_selec_SVM_radial_001_GA_acc  else character(0),
  "Radial_0.1"   = if (exists("snps_selec_SVM_radial_01_GA_acc"))   snps_selec_SVM_radial_01_GA_acc   else character(0),
  "Radial_1"     = if (exists("snps_selec_SVM_radial_1_GA_acc"))    snps_selec_SVM_radial_1_GA_acc    else character(0)
)

## 3) Cria dataframes "rank_global_selec" (RF + valor-p) para cada seleção e salva
dir.create("rank_global_classif", showWarnings = FALSE)
rank_global_selec_CLASSIF <- list()

for (lbl in names(selec_por_kernel)) {
  snps_sel <- selec_por_kernel[[lbl]]
  if (length(snps_sel) > 0) {
    rf_part  <- rank_RF_df[snps_sel, , drop = FALSE]
    pv_part  <- pvals[snps_sel, , drop = FALSE]
    df_out   <- cbind(rf_part, pv_part)
    rank_global_selec_CLASSIF[[lbl]] <- df_out
    write.csv(df_out,
              file = file.path("rank_global_classif",
                               paste0("rank_global_CLASSIF_", lbl, ".txt")),
              row.names = TRUE)
  }
}

## 4) União e Interseção já computadas no seu SMS:
#   uniao_final, intersecao_final
# Aqui apenas garantimos colapsar para impressão bonita.
collapse_or_label <- function(x) {
  if (length(x) == 0L) return("Conjunto vazio")
  if (is.character(x) && length(x) == 1L && x == "Conjunto vazio") return(x)
  paste(x, collapse = ", ")
}

## 5) “Medidas de Otimalidade” dos SNPs causais (CLASSIFICAÇÃO)
##    — análogo ao bloco do SVR/Bessel, mas com métricas de classificação.
##    Ajuste 'var_sel_causais' conforme seu experimento:
var_sel_causais <- c("SNP1","SNP2","SNP3","SNP4","SNP5","SNP6","SNP7","SNP8")

# Fun auxiliar para extrair métricas nomeadas do validacao_cruzada_cls
pega_metricas_cls <- function(res_vec) {
  # Ordem retornada pelo seu validacao_cruzada_cls():
  # acc_mean, acc_sd, prec_mean, prec_sd, rec_mean, rec_sd, f1_mean, f1_sd
  out <- c(
    ACC   = unname(res_vec["acc_mean"]),
    PREC  = unname(res_vec["prec_mean"]),
    REC   = unname(res_vec["rec_mean"]),
    F1    = unname(res_vec["f1_mean"])
  )
  return(out)
}

# Avalia métricas com kernels: linear e radial (gammas = 0.001, 0.01, 0.1, 1)
metricas_causais <- list()

if (all(var_sel_causais %in% names(dados[[1]]))) {
  base_causais <- dados[[1]][c(var_sel_causais, "fenotipo")]
  
  # Linear
  res_lin <- validacao_cruzada_cls(
    data = base_causais, folds = 10, gamma = 0.01, cost = 1, kernel = "linear"
  )
  metricas_causais[["Linear"]] <- pega_metricas_cls(res_lin)
  
  # Radial 0.001
  res_r1 <- validacao_cruzada_cls(
    data = base_causais, folds = 10, gamma = 0.001, cost = 1, kernel = "radial"
  )
  metricas_causais[["Radial_0.001"]] <- pega_metricas_cls(res_r1)
  
  # Radial 0.01
  res_r2 <- validacao_cruzada_cls(
    data = base_causais, folds = 10, gamma = 0.01, cost = 1, kernel = "radial"
  )
  metricas_causais[["Radial_0.01"]] <- pega_metricas_cls(res_r2)
  
  # Radial 0.1
  res_r3 <- validacao_cruzada_cls(
    data = base_causais, folds = 10, gamma = 0.1, cost = 1, kernel = "radial"
  )
  metricas_causais[["Radial_0.1"]] <- pega_metricas_cls(res_r3)
  
  # Radial 1
  res_r4 <- validacao_cruzada_cls(
    data = base_causais, folds = 10, gamma = 1.0, cost = 1, kernel = "radial"
  )
  metricas_causais[["Radial_1"]] <- pega_metricas_cls(res_r4)
}

## 6) Relatório consolidado (formato semelhante ao Bessel)
# ---- Calcula tempo total ----
t_total_elapsed <- round((proc.time() - t_total_inicio)[["elapsed"]])

sink("my_output_SMS_CLASSIF_ACC_2.txt")

cat("===== Seleções por Kernel (GA, Classificação) =====\n")
cat("SMS Linear = ",          collapse_or_label(selec_por_kernel[["Linear"]]),        "\n")
cat("SMS Radial 0.001 = ",    collapse_or_label(selec_por_kernel[["Radial_0.001"]]), "\n")
cat("SMS Radial 0.01 = ",     collapse_or_label(selec_por_kernel[["Radial_0.01"]]),  "\n")
cat("SMS Radial 0.1 = ",      collapse_or_label(selec_por_kernel[["Radial_0.1"]]),   "\n")
cat("SMS Radial 1 = ",        collapse_or_label(selec_por_kernel[["Radial_1"]]),     "\n\n")

cat("===== União / Interseção =====\n")
cat("União = ",       collapse_or_label(uniao_final),       "\n")
cat("Interseção = ",  collapse_or_label(intersecao_final),  "\n\n")

cat("===== Seleção por Valor-p (binário) =====\n")
cat("Valor-p bruto (<= 0.05) = ",    collapse_or_label(pvals_bruto_selec),    "\n")
cat("Valor-p ajustado (<= 0.05) = ", collapse_or_label(pvals_ajustado_selec), "\n\n")

cat("===== Medidas de Otimalidade (SNPs causais, Classificação) =====\n")
if (length(metricas_causais) == 0L) {
  cat("Conjunto 'var_sel_causais' não encontrado nas colunas. Ajuste os nomes e reexecute.\n")
} else {
  for (lbl in names(metricas_causais)) {
    m <- metricas_causais[[lbl]]
    cat(lbl, " ->  ACC =", sprintf("%.4f", m["ACC"]),
        " | PREC =", sprintf("%.4f", m["PREC"]),
        " | REC  =", sprintf("%.4f", m["REC"]),
        " | F1   =", sprintf("%.4f", m["F1"]), "\n")
  }
}


cat("\n===== Tempos de Execução (segundos) =====\n")
cat("\n===== Tempos de Execução Detalhados (segundos) =====\n")
cat(sprintf("  %-24s: %8.1f s\n", "Simulação",     t_sim))
cat(sprintf("  %-24s: %8.1f s\n", "Random Forest", t_rf))
cat("\n  -- Corte (curva na trilha RF) --\n")
cat(sprintf("  Corte %-18s: %8.1f s\n", "Linear", t_corte[["Linear"]]))
cat(sprintf("  Corte %-18s: %8.1f s\n", "Radial_0.001", t_corte[["Radial_0.001"]]))
cat(sprintf("  Corte %-18s: %8.1f s\n", "Radial_0.01", t_corte[["Radial_0.01"]]))
cat(sprintf("  Corte %-18s: %8.1f s\n", "Radial_0.1", t_corte[["Radial_0.1"]]))
cat(sprintf("  Corte %-18s: %8.1f s\n", "Radial_1.0", t_corte[["Radial_1.0"]]))
cat("\n  -- GA (refinamento) --\n")
cat(sprintf("  GA    %-18s: %8.1f s\n", "Linear", t_ga[["Linear"]]))
cat(sprintf("  GA    %-18s: %8.1f s\n", "Radial_0.001", t_ga[["Radial_0.001"]]))
cat(sprintf("  GA    %-18s: %8.1f s\n", "Radial_0.01", t_ga[["Radial_0.01"]]))
cat(sprintf("  GA    %-18s: %8.1f s\n", "Radial_0.1", t_ga[["Radial_0.1"]]))
cat(sprintf("  GA    %-18s: %8.1f s\n", "Radial_1.0", t_ga[["Radial_1.0"]]))
cat("\n  -- Kernel total (Corte+GA) --\n")
cat(sprintf("  %-24s: %8.1f s\n", "Linear", t_kernels[["Linear"]]))
cat(sprintf("  %-24s: %8.1f s\n", "Radial_0.001", t_kernels[["Radial_0.001"]]))
cat(sprintf("  %-24s: %8.1f s\n", "Radial_0.01", t_kernels[["Radial_0.01"]]))
cat(sprintf("  %-24s: %8.1f s\n", "Radial_0.1", t_kernels[["Radial_0.1"]]))
cat(sprintf("  %-24s: %8.1f s\n", "Radial_1.0", t_kernels[["Radial_1.0"]]))
cat(sprintf("  %-24s: %8.1f s\n", "TOTAL", t_total_elapsed))
