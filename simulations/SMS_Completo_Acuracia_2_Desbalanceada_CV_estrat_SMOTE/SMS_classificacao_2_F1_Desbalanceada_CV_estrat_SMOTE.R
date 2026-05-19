# ----------------------------------------------------------------------------
# Bloco de configuração inicial: pasta de trabalho e instalação de pacotes.
# ----------------------------------------------------------------------------

# Define a pasta de trabalho (working directory) como a pasta do próprio script.
# Funciona em três cenários:
#   1) RStudio: usa rstudioapi para localizar o .R atualmente aberto.
#   2) Rscript: lê o argumento --file= passado pelo Rscript ao processo.
#   3) Fallback: mantém o getwd() atual.
.script_dir <- tryCatch({
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path)) {
    dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    .args <- commandArgs(trailingOnly = FALSE)
    .file_arg <- sub("^--file=", "", .args[grep("^--file=", .args)])
    if (length(.file_arg) > 0) dirname(normalizePath(.file_arg)) else getwd()
  }
}, error = function(e) getwd())
setwd(.script_dir)

# ---- Timer global ----
t_total_inicio <- proc.time()
t_kernels      <- list()    # lista para tempos por kernel
t_corte        <- list()  # tempo da etapa de corte por kernel
t_ga           <- list()  # tempo do GA por kernel
t_sim          <- 0       # tempo da simulação da base
t_rf           <- 0       # tempo da Random Forest
cat("[setwd] pasta de trabalho:", getwd(), "\n")

# Lista de pacotes necessários:
#   scrime       -> simulateSNPglm (simulação de SNPs/fenótipo)
#   e1071        -> SVM (Support Vector Machine)
#   kernlab      -> kernels alternativos (não usado direto, mas mantido)
#   randomForest -> RF + medidas de importância (MeanDecreaseGini)
#   doParallel   -> paralelismo (usado pelo GA com parallel = TRUE)
#   GA           -> algoritmo genético para refinar a seleção de SNPs
#   ggplot2      -> gráficos da evolução do GA
#   recipes      -> pipeline tidymodels (usado pelo SMOTE sem vazamento)
#   themis       -> step_smote / step_smotenc (oversampling no treino)
packages <- c("scrime","e1071","kernlab","randomForest",
              "doParallel", "GA", "ggplot2", "recipes", "themis")

# Vetor lógico TRUE/FALSE indicando, para cada nome em 'packages', se já está instalado.
installed_packages <- packages %in% rownames(installed.packages())
# Se houver pelo menos um FALSE, instala apenas os que faltam (não reinstala os existentes).
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Carrega 'scrime' agora porque suas funções (simulateSNPglm) são usadas logo abaixo.
# Os demais pacotes são carregados perto do bloco onde são usados (na seção CLASSIFICAÇÃO).
library(scrime)

#------------------------------------------------------------
# Função auxiliar interna (prefixo "." sinaliza uso privado):
# Dado um intercepto 'beta0', simula uma base via simulateSNPglm e devolve
# a proporção de indivíduos classificados como 1 (i.e., linpred > 0).
# A regra "linpred > 0" equivale a "prob > 0.5" no modelo logístico.
#------------------------------------------------------------
.prop_ones_given_beta0 <- function(beta0,            # intercepto candidato
                                   n.obs, n.snp,     # nº de indivíduos e nº de SNPs
                                   list.snp, list.ia,# SNPs causais e suas interações (ver scrime)
                                   beta, maf,        # efeitos lineares e faixa de MAF (Minor Allele Frequency)
                                   err.fun = rnorm,  # distribuição do ruído (gaussiano por padrão)
                                   rand = 123,       # semente base
                                   reps = 1) {       # nº de réplicas p/ reduzir ruído estocástico
  # Acumulador da proporção média de 1s ao longo das réplicas.
  acc <- 0
  for (i in 1:reps) {
    # Simula uma base com o pacote scrime; muda a semente em cada réplica
    # somando 'i' a 'rand' (assim 'reps' bases são diferentes mas reproduzíveis).
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
    # sim$y aqui é o preditor linear (logit), não o rótulo binário.
    linpred <- sim$y
    # mean(linpred > 0) -> fração de indivíduos com probabilidade > 0.5.
    acc <- acc + mean(linpred > 0)
  }
  # Devolve a média das proporções (estimativa mais estável).
  acc / reps
}

#------------------------------------------------------------
# Procura o valor de 'beta0' (intercepto do logito) que faz a proporção
# de 1s simulada chegar próxima de 'target_prop' (ex.: 0.20, 0.50, 0.80).
# Estratégia: usa 'uniroot' (busca de raiz) sobre f(b0) = prop(b0) - target_prop.
# Se as bordas iniciais (lower/upper) não cercarem a raiz, alarga o intervalo
# até 'max_expand' vezes (em incrementos de ±10) antes de desistir.
#------------------------------------------------------------
find_beta0_for_target <- function(target_prop,        # proporção desejada de 1s (entre 0 e 1)
                                  n.obs, n.snp,       # dimensões da base simulada
                                  list.snp, list.ia,  # SNPs causais + interações
                                  beta, maf,          # efeitos e MAFs
                                  err.fun = rnorm,    # distribuição do ruído
                                  rand = 123,         # semente base
                                  reps = 3,           # nº de réplicas usadas em cada avaliação
                                  lower = -20, upper = 20,  # intervalo inicial de busca em beta0
                                  tol = 1e-3, max_expand = 5) { # tolerância e nº máx. de expansões
  # Sanidade: target_prop tem que ser uma proporção válida (0 < p < 1).
  if (target_prop <= 0 || target_prop >= 1) {
    stop("target_prop deve estar em (0,1).")
  }

  # Função objetivo: queremos a raiz de f(b0) = prop(b0) - target_prop.
  # Quando f(b0) = 0, a proporção simulada bate com o alvo.
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

  # Inicializa os limites do bracket (L, U) e avalia f nas bordas.
  L <- lower; U <- upper
  fL <- f(L); fU <- f(U)

  # 'uniroot' só funciona se f(L) e f(U) têm sinais opostos.
  # Se não tem (mesmo sinal), alarga o intervalo simetricamente ±10 a cada tentativa,
  # repetindo até 'max_expand' vezes.
  k <- 0
  while (sign(fL) == sign(fU) && k < max_expand) {
    L <- L - 10; U <- U + 10
    fL <- f(L);  fU <- f(U)
    k <- k + 1
  }

  # Se mesmo depois das expansões os sinais coincidem, não há como localizar a raiz.
  if (sign(fL) == sign(fU)) {
    stop("Não foi possível bracketing: tente ajustar 'lower/upper' ou aumentar 'max_expand'.")
  }

  # 'uniroot' devolve uma lista; '$root' é o valor de beta0 que zera f (até a tolerância).
  uniroot(f, interval = c(L, U), tol = tol)$root
}

#------------------------------------------------------------
# Função "de alto nível": gera uma base simulada com:
#   - n1 indivíduos de classe 1 e n0 de classe 0 (proporção alvo = n1/(n1+n0));
#   - se exact_by_rank=TRUE, força EXATAMENTE essas contagens via ranking do logito;
#   - se exact_by_rank=FALSE, aceita a classificação natural (linpred > 0),
#     ficando próximo da proporção desejada por ajuste de beta0.
# Retorna lista com a base + beta0 estimado + diagnósticos.
#------------------------------------------------------------
simulate_scrime_binary_with_ratio <- function(n1, n0,                # quantidades alvo por classe
                                              n.snp,                 # nº de SNPs simulados
                                              list.snp, list.ia,     # SNPs causais e suas interações
                                              beta = rep(2, length(list.snp)), # efeitos lineares
                                              maf  = c(0.1, 0.4),    # faixa de MAF para sortear cada SNP
                                              err.fun = rnorm,       # distribuição do ruído
                                              rand = 123,            # semente base
                                              reps_beta0 = 3,        # réplicas usadas na busca de beta0
                                              exact_by_rank = FALSE) {  # se TRUE, força exatamente n1/n0
  # Total de observações e proporção-alvo p = n1/(n1+n0).
  n.obs <- n1 + n0
  target <- n1 / n.obs

  # 1) Acha o intercepto beta0 que produz, em média, a proporção desejada de 1s.
  #    'find_beta0_for_target' usa busca de raiz (uniroot) sobre simulações repetidas.
  beta0_star <- find_beta0_for_target(
    target_prop = target,
    n.obs = n.obs, n.snp = n.snp,
    list.snp = list.snp, list.ia = list.ia,
    beta = beta, maf = maf,
    err.fun = err.fun, rand = rand,
    reps = reps_beta0,
    lower = -20, upper = 20, tol = 1e-3
  )

  # 2) Simulação final, agora com o beta0 ajustado (uma única chamada, semente 'rand').
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

  # sim$x -> matriz dos SNPs (n.obs x n.snp); convertemos para data.frame.
  X <- as.data.frame(sim$x)
  # sim$y -> preditor linear (logito) por indivíduo (não é probabilidade nem rótulo).
  linpred <- sim$y
  # plogis(linpred) = 1 / (1 + exp(-linpred)) -> probabilidade da classe 1.
  prob <- plogis(linpred)

  if (!exact_by_rank) {
    # "Classificação natural": rótulo 1 se prob > 0.5, equivalente a linpred > 0.
    # Aqui o nº exato de 1s pode variar um pouco em torno de n1 (mas perto, pois ajustamos beta0).
    Y <- as.numeric(linpred > 0)
  } else {
    # 'order(linpred, decreasing=TRUE)' devolve os índices em ordem decrescente de logito.
    ord <- order(linpred, decreasing = TRUE)
    # Vetor inicial todo zero, do tamanho da base.
    Y <- integer(length(linpred))
    # Marca como 1 exatamente os 'n1' indivíduos com maior logito (top-n1 = classe positiva).
    Y[ord[seq_len(n1)]] <- 1
  }

  # Monta o data.frame final: X + colunas auxiliares (linpred, prob) + fenotipo.
  # Mantemos linpred/prob para auditoria; serão removidos antes do SVM.
  base <- data.frame(X, linpred = linpred, prob = prob, fenotipo = Y)
  # Lista de saída com a base + parâmetros diagnósticos da simulação.
  out <- list(
    data            = base,                # base simulada (data.frame)
    beta0           = beta0_star,          # intercepto efetivamente usado
    achieved_prop_1 = mean(Y == 1),        # proporção realmente atingida de classe 1
    counts          = table(Y)             # contagens por classe
  )
  return(out)
}

#------------------------------------------------------------
# EXEMPLOS DE USO
# Gera duas bases (A: ~60/40 aproximado; B: exatamente 200/800) e
# usa a base B (desbalanceada) como entrada do pipeline SMS.
#------------------------------------------------------------

# Parâmetros do experimento:
num_individuos <- 1000               # tamanho amostral total (referência)
num_snp        <- 100                # total de SNPs simulados por indivíduo
list.snp <- list(1, 2, 3, c(4,5), c(6,7,8))    # SNPs CAUSAIS (índices das colunas)
list.ia  <- list(1, 3, 2, c(-1,3), c(1,2,3)) # padrão de interações entre causais (ver scrime::simulateSNPglm)
beta     <- c(2, 2, 2, 3, 4)                # tamanho do efeito (igual para os 8 causais)
maf      <- c(0.1, 0.4)              # MAF é sorteado uniformemente nesse intervalo por SNP

# ------- A) Caso aproximado: ~60% classe 1 / 40% classe 0 -------
set.seed(123)  # semente externa para reprodutibilidade (a interna do scrime é 'rand')
res_approx <- simulate_scrime_binary_with_ratio(
  n1 = 600, n0 = 400,                # alvos aproximados
  n.snp = num_snp,
  list.snp = list.snp, list.ia = list.ia,
  beta = beta, maf = maf,
  rand = 123,
  reps_beta0 = 3,           # 3 a 5 réplicas reduzem ruído na busca de beta0
  exact_by_rank = FALSE     # NÃO força exatidão: a proporção fica próxima do alvo
)
res_approx$beta0                # intercepto efetivamente usado
res_approx$counts               # contagens reais por classe (próximas de 600/400)
mean(res_approx$data$fenotipo)  # proporção de 1s atingida (próxima de 0.6)

# ------- B) Caso exato (desbalanceado 20/80) -------
set.seed(123)
t_sim_ini <- proc.time()
res_exact <- simulate_scrime_binary_with_ratio(
  n1 = 200, n0 = 800,                # alvos exatos: 200 classe 1, 800 classe 0
  n.snp = num_snp,
  list.snp = list.snp, list.ia = list.ia,
  beta = beta, maf = maf,
  rand = 123,
  reps_beta0 = 3,
  exact_by_rank = TRUE       # FORÇA exatamente n1/n0 via top-n do logito
)
res_exact$beta0
res_exact$counts            # deve mostrar 800/200 (na ordem 0,1)
# Olhada rápida: 8 primeiros SNPs (os causais) + última coluna (fenotipo).
head(res_exact$data[, c(1:8, ncol(res_exact$data))])

# Persiste em CSV para auditoria/eventual reuso fora do R.
# Atenção: 'res_exact' é uma LISTA; write.csv vai serializar parte dela como colunas.
write.csv(res_exact, file = "dados.csv", row.names = FALSE)


# ----------------------------------------------------------------------------
# Empacotamento dos dados para o pipeline SMS.
# ----------------------------------------------------------------------------
# 'dados' será uma lista; dados[[1]] = base completa (SNPs + fenotipo),
# e dados[[2]] (mais adiante) = subconjunto após o corte do rank da RF.
dados <- list()

# (Alternativa: ler de CSV; mantido comentado.)
#dados_inicial <- read.csv("dados_simulados.csv")

# 'res_exact$data' é o data.frame retornado pela simulação.
# 'subset(..., select = c(-linpred, -prob))' remove essas duas colunas auxiliares
# (eram úteis apenas para auditoria do logito/probabilidade na simulação).
dados[[1]] <- subset(res_exact[[1]], select = c(-linpred, -prob))
t_sim <- round((proc.time() - t_sim_ini)[["elapsed"]])
cat(sprintf("[TIMER] Simulação     : %8.1f s\n", t_sim))

# Sanidade visual: olha as primeiras linhas (espera-se SNPs + fenotipo).
head(dados[[1]])

# Sanidade numérica: distribuição da classe (esperado 800 zeros e 200 uns).
table(dados[[1]]$fenotipo)


#Codigo SMS  (Início do pipeline: helpers, CV estratificada, SMOTE, RF + GA + SVM)



# --- util: converte colunas list/char/factor para numéricas atômicas ---
# Necessário porque o e1071::svm precisa de matriz numérica e o pacote scrime
# pode devolver colunas em formatos não-numéricos (ex.: lista de comprimento 1
# por célula). Esta função padroniza tudo, sem alterar a coluna alvo.
sanitize_predictors <- function(df, target = "fenotipo") {
  # Nomes das colunas preditoras (tudo menos o alvo).
  pcols <- setdiff(names(df), target)
  for (nm in pcols) {
    # Caso 1: a coluna é uma LISTA (cada célula é um vetor possivelmente unitário).
    if (is.list(df[[nm]])) {
      # 'vapply' garante saída numérica do comprimento 1; segura erros de tipos heterogêneos.
      df[[nm]] <- vapply(df[[nm]], function(z) {
        if (length(z) == 1) {
          # Caso normal: 1 valor por célula -> converte para numeric.
          as.numeric(z)
        } else {
          # Caso anômalo: nenhum ou mais de um valor -> marca como NA numérico.
          NA_real_
        }
      }, numeric(1))
    }
    # Caso 2: a coluna é FATOR (ex.: "0","1","2") -> passa por character p/ não pegar o código interno.
    if (is.factor(df[[nm]])) {
      df[[nm]] <- as.numeric(as.character(df[[nm]]))
    } else if (is.character(df[[nm]])) {
      # Caso 3: a coluna é CHARACTER -> coerce direto; 'suppressWarnings' ignora "NAs introduced".
      suppressWarnings({
        tmp <- as.numeric(df[[nm]])
      })
      df[[nm]] <- tmp
    }
  }
  # Devolve o data.frame com preditores já numéricos; alvo permanece intacto.
  df
}

# --- prepara alvo e remove 'prob' se existir ---
# Garante que o data.frame está pronto para classificação binária:
#   - sem coluna 'prob' (resíduo da simulação);
#   - 'fenotipo' como fator com níveis fixos c("0","1") na ordem certa.
prepara_cls <- function(df) {
  # Remove a coluna 'prob' (probabilidade simulada) se ainda estiver presente.
  if ("prob" %in% names(df)) df$prob <- NULL
  # Em ambos os ramos forçamos níveis c("0","1") — isso garante que predict()
  # e a matriz de confusão tenham SEMPRE as duas classes na mesma ordem.
  if (!is.factor(df$fenotipo)) {
    df$fenotipo <- factor(df$fenotipo, levels = c(0,1))
  } else {
    df$fenotipo <- factor(df$fenotipo, levels = c(0,1))
  }
  df
}

# --- SMOTE apenas no treino (chamar dentro de cada fold após o split; evita vazamento) ---
# Mantém os SNPs em valores discretos (mesmo conjunto {1,2,3,...} dos originais),
# evitando o viés de interpolação contínua do SMOTE clássico.
# 'mode':
#   "smotenc" -> themis::step_smotenc tratando cada SNP como FATOR;
#                para cada SNP do sintético, atribui a MODA entre os K vizinhos.
#                (Padrão; mais coerente para dados categóricos/ordinais.)
#   "round"   -> themis::step_smote (interpolação contínua) + "snap" de cada
#                valor sintético para o inteiro mais próximo observado na
#                respectiva coluna do treino daquele fold.
# 'save_path' (opcional): grava o treino pós-SMOTE em CSV com a coluna
# '.smote_sintetico' (FALSE = original, TRUE = sintético). A coluna NÃO é
# repassada ao SVM.
aplica_smote_treino <- function(trainset, over_ratio = 1, seed = NULL,
                                mode = c("smotenc", "round"),
                                save_path = NULL) {
  # 'match.arg' valida que 'mode' é uma das opções permitidas e seleciona a 1ª se faltar.
  mode <- match.arg(mode)
  # Semente para reprodutibilidade da geração de sintéticos dentro deste fold.
  if (!is.null(seed)) set.seed(seed)

  # ---- utilitário interno: grava o treino em CSV (apenas se 'save_path' definido) ----
  # 'marca_origem_FALSE=TRUE' indica que o data.frame ainda não tem '.smote_sintetico'
  # (caso em que estamos devolvendo o treino original sem SMOTE).
  .salvar_se_preciso <- function(df, marca_origem_FALSE = TRUE) {
    # Sem caminho de saída -> nada a fazer.
    if (is.null(save_path)) return(invisible(NULL))
    # Cópia para não mutar o argumento original.
    df_save <- df
    # Adiciona '.smote_sintetico = FALSE' (todas as linhas são originais).
    if (marca_origem_FALSE && !".smote_sintetico" %in% names(df_save)) {
      df_save$.smote_sintetico <- FALSE
    }
    # Cria a pasta de destino se ainda não existir.
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    # Grava o CSV (sem nomes de linha).
    utils::write.csv(df_save, save_path, row.names = FALSE)
  }

  # Distribuição de classes no treino do fold atual.
  ytab <- table(trainset$fenotipo)
  # Se houver < 2 classes ou a minoritária tiver < 2 amostras, SMOTE não é viável
  # (precisa de pelo menos 'k+1' = 6 vizinhos, mas evitamos erro adicional aqui).
  if (length(ytab) < 2L || min(ytab) < 2L) {
    # Apenas salva o treino original (marcando todas as linhas como originais) e retorna.
    .salvar_se_preciso(trainset)
    return(trainset)
  }

  # Nomes das colunas preditoras (todas, menos o alvo 'fenotipo').
  pred_cols <- setdiff(names(trainset), "fenotipo")

  if (mode == "smotenc") {
    # ---------- MODO SMOTE-NC: trata SNPs como fatores; mantém valores discretos ----------
    # Cópia do treino que será convertida coluna a coluna em fatores.
    train_fac <- trainset
    # Estrutura auxiliar para guardar quais níveis cada SNP tem (não usada depois,
    # mas útil para depuração).
    levs_por_col <- vector("list", length(pred_cols))
    names(levs_por_col) <- pred_cols
    for (nm in pred_cols) {
      v <- trainset[[nm]]
      # Se já vier como fator, converte para número para padronizar.
      if (is.factor(v)) v <- as.numeric(as.character(v))
      # Níveis distintos observados na coluna, em ordem crescente (ex.: 1,2,3).
      lv <- sort(unique(stats::na.omit(v)))
      levs_por_col[[nm]] <- lv
      # Converte a coluna em fator com esses níveis (rotulados como "1","2","3").
      train_fac[[nm]] <- factor(v, levels = as.character(lv))
    }

    # Receita do tidymodels: fórmula 'fenotipo ~ .' sobre os dados em fator.
    rec <- recipes::recipe(fenotipo ~ ., data = train_fac)
    # Acrescenta o passo SMOTE-NC: distância de Gower + MODA dos vizinhos para
    # variáveis categóricas. O 'themis' (>= 0.2 e <= 1.0.3) NÃO aceita o
    # argumento 'indicator_column', então não tentamos passá-lo. Em vez disso,
    # marcamos sintéticos por posição mais abaixo (themis sempre devolve os
    # originais nas primeiras nrow(train_fac) linhas e os sintéticos depois).
    rec <- themis::step_smotenc(rec, fenotipo,
                                over_ratio = over_ratio)

    # 'prep' aprende os parâmetros do passo só sobre o treino;
    # 'bake(new_data = NULL)' devolve o próprio treino balanceado.
    # Capturamos a mensagem de erro para que falhas reais não passem em silêncio.
    .smote_err <- NULL
    out <- tryCatch(
      recipes::bake(recipes::prep(rec, training = train_fac), new_data = NULL),
      error = function(e) { .smote_err <<- conditionMessage(e); NULL }
    )
    # Falha no SMOTE -> avisa, mantém o treino original e sai (sem balancear).
    if (is.null(out)) {
      warning("SMOTE-NC falhou neste fold: ", .smote_err, call. = FALSE)
      .salvar_se_preciso(trainset)
      return(trainset)
    }
    # Converte o tibble devolvido em data.frame para manipulação posterior.
    out <- as.data.frame(out)

    # SMOTE-NC devolve as colunas como fator; revertemos para inteiro
    # (mesmo tipo das colunas originais — facilita o SVM e o CSV).
    for (nm in pred_cols) {
      out[[nm]] <- as.integer(as.character(out[[nm]]))
    }

    # Marca origem das linhas: as primeiras nrow(trainset) são originais
    # (themis::smotenc faz rbind(originais, sintéticos) preservando a ordem)
    # e o restante são sintéticas geradas pelo SMOTE-NC.
    n_orig <- nrow(trainset)
    out$.smote_sintetico <- c(rep(FALSE, n_orig),
                              rep(TRUE,  nrow(out) - n_orig))
  } else {
    # ---------- MODO ROUND: SMOTE clássico (interpolação) + snap aos inteiros ----------
    # Receita igual à do modo SMOTE-NC, mas usando 'step_smote' (interpolação contínua).
    # Também aqui o themis (<= 1.0.3) NÃO aceita 'indicator_column'; marcamos
    # sintéticos por posição (originais nas primeiras nrow(trainset) linhas).
    rec <- recipes::recipe(fenotipo ~ ., data = trainset)
    rec <- themis::step_smote(rec, fenotipo,
                              over_ratio = over_ratio)
    # Prep + bake sobre o treino, com proteção contra erro (capturando a mensagem).
    .smote_err <- NULL
    out <- tryCatch(
      recipes::bake(recipes::prep(rec, training = trainset), new_data = NULL),
      error = function(e) { .smote_err <<- conditionMessage(e); NULL }
    )
    if (is.null(out)) {
      warning("SMOTE (round) falhou neste fold: ", .smote_err, call. = FALSE)
      .salvar_se_preciso(trainset)
      return(trainset)
    }
    out <- as.data.frame(out)

    # Marca origem das linhas: as primeiras nrow(trainset) são originais e o
    # restante são sintéticas. themis::smote() faz rbind(originais, sintéticos).
    n_orig <- nrow(trainset)
    out$.smote_sintetico <- c(rep(FALSE, n_orig),
                              rep(TRUE,  nrow(out) - n_orig))

    # Vetor lógico marcando quais linhas do resultado são sintéticas.
    is_syn <- as.logical(out$.smote_sintetico)
    # Só faz sentido aplicar o snap se houver pelo menos um sintético.
    if (any(is_syn)) {
      for (nm in pred_cols) {
        # Valores originais da coluna no treino daquele fold (sem NAs).
        vals_orig <- stats::na.omit(trainset[[nm]])
        # Conjunto único de valores observados (ex.: c(1,2,3)).
        uniq_orig <- sort(unique(vals_orig))
        # Coluna sem nenhum valor utilizável -> pula.
        if (length(uniq_orig) == 0L) next
        # Valores interpolados (contínuos) das linhas sintéticas nesta coluna.
        v_syn <- out[[nm]][is_syn]
        # Para cada valor sintético, encontra o valor original mais próximo
        # (em distância absoluta) — preserva o suporte discreto observado.
        snap <- vapply(v_syn, function(x) {
          if (is.na(x)) return(NA_real_)
          uniq_orig[which.min(abs(uniq_orig - x))]
        }, numeric(1))
        # Substitui apenas as linhas sintéticas pelos valores "snapped".
        out[[nm]][is_syn] <- snap
      }
    }
  }

  # Grava o CSV do treino pós-SMOTE com a coluna .smote_sintetico (se 'save_path').
  .salvar_se_preciso(out, marca_origem_FALSE = FALSE)

  # Remove a coluna indicadora antes de devolver para o SVM (não é preditor).
  out$.smote_sintetico <- NULL
  # Padroniza tipos (numéricos atômicos) e formato do alvo (fator c("0","1")).
  sanitize_predictors(prepara_cls(out), target = "fenotipo")
}

# --- validação cruzada estratificada + opcional SMOTE só no treino de cada fold ---
# Parâmetros relevantes:
#   smote_mode   : "smotenc" (padrão) ou "round" — ver 'aplica_smote_treino'.
#   save_smote   : se TRUE, grava o treino pós-SMOTE de cada fold em CSV.
#   smote_outdir : pasta de saída dos CSVs por fold.
#   smote_label  : rótulo a ser usado no nome dos arquivos.
validacao_cruzada_cls <- function(data, folds = 5, gamma = 0.01, cost = 1,
                                  kernel = c("radial", "linear"),
                                  use_smote = TRUE, smote_over_ratio = 1,
                                  smote_mode = c("smotenc", "round"),
                                  cv_seed = 123,
                                  save_smote = FALSE,
                                  smote_outdir = "smote_folds",
                                  smote_label = "run") {
  # Valida argumentos categóricos e fixa um único valor (1º se omitido).
  kernel     <- match.arg(kernel)
  smote_mode <- match.arg(smote_mode)
  # 'prepara_cls' garante que 'fenotipo' é fator com níveis c("0","1") e remove 'prob'.
  data <- prepara_cls(data)
  # 'sanitize_predictors' converte preditores list/char/factor em numéricos atômicos.
  data <- sanitize_predictors(data, target = "fenotipo")

  # Número de observações e índice 1..n usado em todo o particionamento.
  n     <- nrow(data)
  index <- 1:n

  # Semente única para garantir partição reprodutível em qualquer chamada.
  set.seed(cv_seed)
  # 'idx[i]' guardará o número do fold (1..folds) ao qual a linha i pertence.
  idx <- numeric(n)
  # Níveis de classe (ex.: c("0","1")) — itera por classe para estratificar.
  classes <- levels(as.factor(data$fenotipo))
  for (classe in classes) {
    # Índices das linhas que pertencem à classe atual.
    classe_idx          <- which(data$fenotipo == classe)
    n_classe            <- length(classe_idx)
    # Embaralha os índices dessa classe (estratifica DENTRO da classe).
    classe_idx_shuffled <- sample(classe_idx)
    # Vetor 1,2,...,folds,1,2,... com n_classe elementos: distribui amostras
    # da classe entre os folds em rodízio (proporção ~ igual em cada fold).
    fold_assignments    <- rep(1:folds, length.out = n_classe)
    # Atribui o fold a cada linha da classe (na ordem embaralhada).
    for (j in seq_len(n_classe)) {
      idx[classe_idx_shuffled[j]] <- fold_assignments[j]
    }
  }
  
  # Vetores que guardarão a métrica de cada fold (preenchidos no loop).
  acc  <- numeric(folds)   # acurácia (TP+TN)/(TP+TN+FP+FN)
  prec <- numeric(folds)   # precisão  TP/(TP+FP)
  rec  <- numeric(folds)   # recall    TP/(TP+FN)
  f1   <- numeric(folds)   # F1        2*prec*rec/(prec+rec)

  for (i in seq_len(folds)) {
    # Índices das linhas que vão para teste (i = fold atual) e para treino (todas as outras).
    test_idx  <- which(idx == i)
    train_idx <- which(idx != i)

    # Subconjuntos correspondentes; 'na.omit' descarta linhas com qualquer NA.
    # 'drop = FALSE' garante data.frame mesmo se restar 1 coluna (caso o GA selecione 1 SNP).
    trainset <- na.omit(data[train_idx, , drop = FALSE])
    testset  <- na.omit(data[test_idx , , drop = FALSE])

    # Reaplica o sanitizador (algumas operações como na.omit/subset podem reverter tipos).
    trainset <- sanitize_predictors(trainset, target = "fenotipo")
    testset  <- sanitize_predictors(testset , target = "fenotipo")
    
    if (isTRUE(use_smote)) {
      # Caminho de gravação do CSV pós-SMOTE deste fold (NULL = não grava).
      sp_fold <- NULL
      if (isTRUE(save_smote)) {
        # Nome do arquivo: smote_<label>_fold01.csv ... fold10.csv (zero à esquerda).
        sp_fold <- file.path(smote_outdir,
                             sprintf("smote_%s_fold%02d.csv", smote_label, i))
      }
      # SMOTE aplicado SOMENTE ao 'trainset'; 'testset' fica intacto (sem vazamento).
      # Semente derivada do fold para reprodutibilidade independente entre folds.
      trainset <- aplica_smote_treino(trainset, over_ratio = smote_over_ratio,
                                      seed      = cv_seed + 1000L * i,
                                      mode      = smote_mode,
                                      save_path = sp_fold)
    }
    
    # Constrói matrizes numéricas X (preditores) e vetor y (alvo) para o SVM.
    # Usar a interface (x = ..., y = ...) é mais robusta que a fórmula, especialmente
    # quando o número de colunas muda dinamicamente (como no GA).
    X_train <- as.matrix(trainset[, setdiff(names(trainset), "fenotipo"), drop = FALSE])
    y_train <- trainset$fenotipo

    X_test  <- as.matrix(testset[, setdiff(names(testset), "fenotipo"), drop = FALSE])
    # Garante factor c("0","1") para alinhar com 'y_pred' (necessário para a matriz de confusão).
    y_true  <- factor(testset$fenotipo, levels = c("0","1"))

    # Treina SVM no conjunto de TREINO (já balanceado por SMOTE se 'use_smote=TRUE').
    svm_fit <- svm(
      x = X_train, y = y_train,
      kernel = kernel,        # "linear" ou "radial"
      gamma  = gamma,          # parâmetro do kernel radial (ignorado se kernel='linear')
      cost   = cost,           # custo de classificação errada (parâmetro C)
      type   = "C-classification",  # tarefa de classificação binária
      scale  = TRUE            # normaliza features internamente (média 0, var 1)
    )

    # Predição no conjunto de teste e alinhamento de níveis para a tabela de confusão.
    y_pred <- predict(svm_fit, X_test)
    y_pred <- factor(y_pred, levels = c("0","1"))

    # Matriz de confusão: linhas = verdade, colunas = predição.
    cm <- table(Truth = y_true, Pred = y_pred)
    # Extração defensiva das células TP/TN/FP/FN: se uma classe não apareceu, devolve 0.
    # Convenção: classe "1" = positiva; classe "0" = negativa.
    TP <- ifelse("1" %in% rownames(cm) && "1" %in% colnames(cm), cm["1","1"], 0) # predisse 1 e era 1
    TN <- ifelse("0" %in% rownames(cm) && "0" %in% colnames(cm), cm["0","0"], 0) # predisse 0 e era 0
    FP <- ifelse("0" %in% rownames(cm) && "1" %in% colnames(cm), cm["0","1"], 0) # predisse 1 mas era 0
    FN <- ifelse("1" %in% rownames(cm) && "0" %in% colnames(cm), cm["1","0"], 0) # predisse 0 mas era 1

    # 'max(1, ...)' protege contra divisão por zero (caso patológico de fold vazio).
    acc[i]  <- (TP + TN) / max(1, TP + TN + FP + FN)
    # Em precisão/recall/F1: se denominador é 0, retornamos 0 em vez de NaN.
    prec[i] <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
    rec[i]  <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    f1[i]   <- ifelse((prec[i] + rec[i]) > 0, 2 * prec[i] * rec[i] / (prec[i] + rec[i]), 0)
  }

  # Retorna VETOR NOMEADO com média e desvio padrão de cada métrica (ordem fixa).
  # O restante do pipeline acessa por nome (ex.: res["acc_mean"]).
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



# ---------------------------------------------------------------------------
# Calcula um valor-p de associação por SNP (fenótipo binário) + correção Bonferroni.
# Para cada SNP: ajusta uma regressão logística simples 'fenotipo ~ SNP' e
# extrai Pr(>|z|) do coeficiente do SNP. Retorna data.frame com p_raw e p_adj.
# ---------------------------------------------------------------------------
valor.p.binario <- function(genotipo_fenotipo) {
  # Convenção: a ÚLTIMA coluna do data.frame de entrada é o fenótipo (0/1 ou fator binário).
  df <- as.data.frame(genotipo_fenotipo)
  # Descarta linhas com qualquer NA antes de modelar (evita erros do glm).
  df <- na.omit(df)

  # Vetor do fenótipo (alvo). 'df[[ncol(df)]]' pega a última coluna.
  y <- df[[ncol(df)]]
  if (is.factor(y)) {
    # Fator: exige exatamente 2 níveis.
    if (length(levels(y)) != 2) stop("Fenótipo fator deve ter 2 níveis.")
    # Reordena para c("0","1") quando os rótulos forem 0/1 (mantém referência consistente).
    if (all(levels(y) %in% c("0","1"))) y <- factor(y, levels = c("0","1"))
  } else {
    # Não é fator: converte char->num se preciso, valida que só tem 0/1 e refatora.
    if (is.character(y)) y <- as.numeric(y)
    if (!all(unique(na.omit(y)) %in% c(0,1))) {
      stop("Fenótipo deve ser 0/1 (ou fator com níveis 0 e 1).")
    }
    y <- factor(y, levels = c(0,1))
  }

  # Vetor de p-valores brutos (1 por SNP). 'ncol(df)-1L' = nº de SNPs (todas menos o alvo).
  p_raw <- numeric(ncol(df) - 1L)
  # Número total de testes (usado na correção de Bonferroni adiante).
  m <- length(p_raw)

  for (i in seq_len(ncol(df) - 1L)) {
    # Vetor do i-ésimo SNP (preditor candidato).
    x <- df[[i]]

    # Garante numérico mesmo se vier como fator/char (0/1/2 são valores típicos).
    if (is.factor(x)) x <- as.numeric(as.character(x))
    if (is.character(x)) x <- suppressWarnings(as.numeric(x))

    # SNP totalmente NA ou monomórfico (1 único valor) -> não há associação testável.
    # Atribuímos p = 1 (sinal de "nada a rejeitar").
    if (all(is.na(x)) || length(unique(na.omit(x))) <= 1L) {
      p_raw[i] <- 1
      next
    }

    # Ajusta regressão logística simples y ~ x via glm com link logit.
    # 'try' isola erros numéricos (ex.: separação completa) sem parar o loop.
    fit <- try(
      glm(y ~ x, family = binomial(), na.action = na.omit),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      # Erro no ajuste -> registra NA (será tratado depois pela correção/ordenamento).
      p_raw[i] <- NA_real_
    } else {
      # 'summary(fit)$coefficients' é uma matriz: linhas = (Intercept, x), colunas
      # = (Estimate, Std. Error, z value, Pr(>|z|)). Pegamos linha 2 (x), coluna 4 (p-valor).
      sm <- summary(fit)$coefficients
      p_raw[i] <- if (!is.na(sm[2, 4])) sm[2, 4] else NA_real_
    }
  }

  # Correção de Bonferroni: p_adj = min(1, m * p_raw). 'pmin' atua elemento a elemento.
  p_adj <- pmin(1, m * p_raw)

  # Empacota em data.frame; 'row.names' usa o nome do SNP (coluna do df de entrada).
  # Os nomes ficam com pontos por causa da regra do R (data.frame não aceita espaços
  # como nome de coluna por padrão): "Valor p bruto" -> "Valor.p.bruto".
  out <- data.frame(
    `Valor p bruto`    = p_raw,
    `Valor p ajustado` = p_adj,
    row.names = names(df)[seq_len(ncol(df) - 1L)]
  )
  return(out)
}

# Executa a função sobre a base completa (dados[[1]]); armazena o resultado em 'pvals'.
pvals <- valor.p.binario(dados[[1]])

# 'order' devolve índices na ordem crescente do vetor — reordena o data.frame por p_raw.
pvals_ordenado_bruto <- pvals[order(pvals$`Valor.p.bruto`), ]

# Mesma ideia, mas pelo p ajustado (Bonferroni).
pvals_ordenado_ajustado <- pvals[order(pvals$`Valor.p.ajustado`), ]

# 'head' exibe as 6 primeiras linhas — usado só para inspeção interativa no console.
head(pvals_ordenado_bruto)
head(pvals_ordenado_ajustado)


################ INÍCIO (CLASSIFICAÇÃO) ################

# Carregamento dos pacotes usados desta seção até o fim do pipeline.
library(randomForest)  # randomForest() para o ranking de importância de SNPs
library(e1071)         # svm() para o classificador
library(GA)            # ga() algoritmo genético (refinamento da seleção)
library(ggplot2)       # gráficos da evolução do GA
library(recipes)       # pipeline tidymodels (usado pelo SMOTE)
library(themis)        # step_smote / step_smotenc (oversampling)


## ---------------- Parâmetros gerais ----------------
# Estruturas que serão preenchidas durante o pipeline (uma posição por kernel SVM):
mean_svm_RF_list <- list()   # média da acurácia do SVM por tamanho de grupo (trilha do rank RF)
GA                <- list()  # objetos GA por kernel (resultado do algoritmo genético)
minimo            <- list()  # legado (mantido por compatibilidade)
corte             <- list()  # ponto de corte ótimo (em nº de SNPs) por kernel
snps_selec_corte  <- list()  # SNPs selecionados pelo corte (pré-GA)
snps_selec_ref    <- list()  # SNPs selecionados pelo GA (refinamento)
percentual_snps   <- 0.95    # fração máxima de SNPs varrida ao construir a curva de acurácia


########## KERNEL LINEAR #########
i <- 1   # índice da "trilha" atual (1=linear, 2=radial 0.001, 3=radial 0.01, ...)


## ---------------- Random Forest (CLASSIFICAÇÃO) ----------------
# Número de árvores da RF; quanto maior, mais estável o ranking (custa mais CPU).
ntree <- 4000

# Cópia do data.frame de entrada (RF aceita data.frame com fator no alvo).
data_temp <- as.data.frame(dados[[1]])
# Garante que 'fenotipo' é fator c("0","1") para a RF tratar como classificação.
if (!is.factor(data_temp$fenotipo)) {
  data_temp$fenotipo <- factor(data_temp$fenotipo, levels = c(0,1))
}

# Semente para reprodutibilidade do bootstrap interno da RF.
set.seed(1)
t_rf_ini <- proc.time()
RF <- randomForest(
  fenotipo ~ ., data = data_temp,    # fórmula: alvo ~ todos os outros (SNPs)
  ntree = ntree,
  mtry = ncol(data_temp) - 1,        # mtry = nº total de preditores (todos sorteados em cada split)
  importance = TRUE                  # ativa cálculo das medidas de importância
)

# 'importance(RF)' devolve matriz com várias medidas (MeanDecreaseAccuracy/Gini, etc).
# Aqui escolhemos MeanDecreaseGini quando disponível; caso contrário, a última coluna.
imp <- importance(RF)
mdg_col <- if ("MeanDecreaseGini" %in% colnames(imp)) "MeanDecreaseGini" else colnames(imp)[ncol(imp)]
# Vetor nomeado de importâncias ORDENADAS do mais ao menos importante.
rank_RF <- sort(imp[, mdg_col], decreasing = TRUE)
t_rf <- round((proc.time() - t_rf_ini)[["elapsed"]])
cat(sprintf("[TIMER] Random Forest : %8.1f s\n", t_rf))
cat("[debug RF] length(rank_RF) =", length(rank_RF),
    "| any NA names =", any(is.na(names(rank_RF))),
    "| any NA vals =", any(is.na(rank_RF)), "\n")

# ---- Timer kernel: Linear ----
t_k_Linear <- proc.time()

## ---------------- SVM (classificação) e k-folds na trilha do rank RF ----------------
# Hiperparâmetros do SVM avaliado neste kernel.
gamma  <- 0.01      # ignorado pelo kernel linear (mantido por consistência)
cost   <- 1.0       # parâmetro C
folds  <- 10        # nº de folds na CV estratificada
kernel <- "linear"  # tipo de kernel desta trilha

# Vetor que receberá a acurácia média para cada tamanho de grupo (10, 20, 30...).
mean_svm_RF_list[[i]] <- numeric()
# Tamanho do passo (de quantos em quantos SNPs aumentamos o grupo a cada iteração).
passo   <- 10
# Quantos passos varremos no máximo: ~95% das colunas / 10.
# 'dados[[1]][, -ncol(dados[[1]])]' = todas as colunas menos a última (fenotipo).
limite  <- floor((length(names(dados[[1]][, -ncol(dados[[1]])])) / passo) * percentual_snps)

# Loop sobre os "grupos" 1, 2, ..., limite+1 (cada grupo é uma curva-ponto).
t_corte_ini <- proc.time()  # inicio corte [Linear]
for (cont in 1:(limite + 1)) {
  cat("Grupo:", cont, "\n")  # log de progresso no console
  j <- cont * passo          # nº de SNPs usados neste ponto da curva (10, 20, 30, ...)

  # Caso especial: 1 único SNP -> 'names(rank_RF)[1]' (não usar :1 que daria seq vazia).
  if (j == 1) {
    var_sel <- names(rank_RF)[1]
  } else {
    var_sel <- head(names(rank_RF), j)  # top-j SNPs do ranking RF (head evita NAs se j > length)
  }

  # Sub-base apenas com os top-j SNPs + fenotipo (preserva nome 'fenotipo' na última col).
  base_cv <- dados[[1]][c(var_sel, "fenotipo")]

  # Roda a CV estratificada + SMOTE no treino e devolve métricas médias.
  # Aqui só consumimos a 1ª métrica (acc_mean) para alimentar a curva.
  svm_cv <- validacao_cruzada_cls(
    data   = base_cv,
    folds  = folds,
    gamma  = gamma,
    cost   = cost,
    kernel = kernel
  )

  f1_mean <- svm_cv["f1_mean"]                   # F1 médio da CV (posição nomeada no vetor)
  mean_svm_RF_list[[i]][cont] <- f1_mean          # armazena na posição 'cont' da trilha 'i'
}

## ---------------- Gráfico (Acurácia do SVM na trilha da RF) ----------------
# Abre o dispositivo PDF; tudo entre pdf(...) e dev.off() é gravado nesse arquivo.
pdf(file = "Grafico_F1_SVM_linear.pdf", height = 5, width = 9)
plot(seq(passo, limite * passo + 10, by = passo),  # eixo X: nº de SNPs em cada ponto
     mean_svm_RF_list[[i]],                         # eixo Y: F1 médio da CV
     type = "o", lwd = 2,                           # 'o' = linha+pontos; espessura 2
     xlab = "Grupo de Marcadores",
     ylab = "F1 (SVM Linear, 10-fold)")
# Ponto de corte = índice do MÁXIMO de acurácia ao longo da trilha.
maximo_idx  <- which.max(mean_svm_RF_list[[i]])
# Converte o índice (1..limite+1) para nº de SNPs (passo * (idx+1)) + buffer de 1 passo.
buffer_alvo <- (maximo_idx + 1) * passo
corte[[i]]  <- min(buffer_alvo, length(rank_RF))
saturou     <- buffer_alvo > length(rank_RF)
# Linha vertical pontilhada marcando o corte (vermelha se buffer saturado em length(rank_RF)).
abline(v   = corte[[i]],
       col = if (saturou) "red" else "black",
       lty = 2,
       lwd = if (saturou) 2     else 1)
if (saturou) {
  legend("topleft",
         legend  = sprintf("corte saturado em %d (buffer pediu %d)",
                           length(rank_RF), buffer_alvo),
         bty = "n", text.col = "red", cex = 0.8)
}
dev.off()  # fecha o PDF

# SNPs selecionados pelo ponto de corte (os top-N do rank RF, com N=corte[[i]]).
snps_selec_corte <- head(names(rank_RF), corte[[i]])
t_corte[["Linear"]] <- round((proc.time() - t_corte_ini)[["elapsed"]])
cat(sprintf("[TIMER] Corte Linear      : %8.1f s\n", t_corte[["Linear"]]))
snps_selec_corte  # imprime no console para inspeção

## ---------------- Seleciona genótipo e monta dados[[2]] após corte ----------------
# 'genotipo[[2]]' = matriz com somente os SNPs do corte (sem o fenotipo).
genotipo <- list()
genotipo[[2]] <- dados[[1]][, head(names(rank_RF), corte[[i]]), drop = FALSE]

# 'dados[[2]]' = data.frame com SNPs do corte + coluna 'fenotipo' do dados[[1]].
dados[[2]] <- cbind(genotipo[[2]], fenotipo = dados[[1]]$fenotipo)
# Reaplica o fator (cbind pode ter quebrado).
if (!is.factor(dados[[2]]$fenotipo)) dados[[2]]$fenotipo <- factor(dados[[2]]$fenotipo, levels = c(0,1))

## ---------------- AG para REFINAMENTO (fitness = acurácia do SVM) ----------------
# Função objetivo do GA: recebe vetor binário 'x' (1 bit por SNP do corte) e
# devolve a acurácia média da CV usando APENAS os SNPs com bit=1.
f <- function(x) {
  inc <- which(x == 1)                     # índices das colunas selecionadas pelo GA
  if (length(inc) == 0) return(0)          # vetor todo zero -> fitness = 0 (penaliza)

  # Sub-base com as colunas escolhidas + última coluna (fenotipo).
  dados_validacao <- dados[[2]][, c(inc, ncol(dados[[2]])), drop = FALSE]

  # Renomeia a última coluna para 'fenotipo' (caso o slicing tenha mudado o nome).
  colnames(dados_validacao)[ncol(dados_validacao)] <- "fenotipo"

  res <- validacao_cruzada_cls(
    data   = dados_validacao,
    folds  = folds,
    gamma  = gamma,
    cost   = cost,
    kernel = kernel
  )
  f1_mean <- res["f1_mean"]                 # F1 médio do conjunto candidato
  return(f1_mean)                           # GA maximiza F1 em vez de acurácia
}
# Alias usado pelo pacote GA (assinatura padrão).
fitness <- function(x) f(x)

# Hiperparâmetros do GA:
run     <- 30    # nº máximo de gerações SEM melhora antes de parar (parada antecipada)
maxiter <- 10    # nº máximo de gerações
pcross  <- 0.8   # probabilidade de cruzamento
pmut    <- 0.1   # probabilidade de mutação
elitism <- 5     # nº de melhores indivíduos copiados intactos para a próxima geração
popSize <- 100   # tamanho da população

# Semente derivada do índice da trilha (i=1 aqui).
t_ga_ini <- proc.time()  # inicio GA [Linear]
set.seed(i)
GA[[i]] <- ga(
  type     = "binary",                        # cromossomos binários (1 bit por SNP)
  fitness  = fitness,
  nBits    = ncol(genotipo[[2]]),             # comprimento do cromossomo = nº de SNPs do corte
  popSize  = popSize,
  names    = colnames(genotipo[[2]]),         # nomes dos bits = nomes dos SNPs (facilita leitura)
  maxiter  = maxiter,
  seed     = i,                               # semente interna do GA (para reprodutibilidade)
  parallel = TRUE,                            # avalia fitness em paralelo (doParallel)
  run      = run,
  pcrossover = pcross,
  pmutation  = pmut,
  elitism    = elitism,
  # Sugestão inicial: 1 cromossomo com TODOS os bits = 1 (usar todos os SNPs do corte).
  # Isso semeia a população com a "solução do corte" e acelera a convergência.
  suggestions = matrix(rep(1, ncol(genotipo[[2]])),
                       ncol = ncol(genotipo[[2]]))
)

# Extrai a 1ª solução da matriz GA@solution (pode ter empate; pegamos a primeira).
sol_mat   <- as.matrix(GA[[i]]@solution)
best_row  <- 1                                 # índice da solução escolhida
sel_cols  <- which(sol_mat[best_row, ] == 1)   # bits ativos -> SNPs selecionados
snps_selec_ref[[i]] <- colnames(sol_mat)[sel_cols]  # nomes dos SNPs do refinamento
t_ga[["Linear"]] <- round((proc.time() - t_ga_ini)[["elapsed"]])
cat(sprintf("[TIMER] GA    Linear      : %8.1f s\n", t_ga[["Linear"]]))


# Alias amigável para o relatório final (1 nome por trilha SVM).
snps_selec_SVM_linear_GA_acc <- snps_selec_ref[[i]]
t_kernels[["Linear"]] <- round((proc.time() - t_k_Linear)[["elapsed"]])
snps_selec_SVM_linear_GA_acc

# Resumo textual + gráfico básico do GA (auxiliar para acompanhar no console).
summary(GA[[i]])
plot(GA[[i]])

# Versão "bonita" do gráfico de evolução do GA via ggplot2.
pdf(file = "Grafico_GA_SVM_linear_F1.pdf", height = 5, width = 9)
# Eixo x = nº de gerações realmente executadas (pode ser menor que maxiter por parada antecipada).
geracao        <- seq_len(GA[[i]]@iter)
# GA[[i]]@summary é uma matriz: col1=best, col2=mean, col3=q1, col4=median, etc.
mean_fitness   <- GA[[i]]@summary[, 2]
median_fitness <- GA[[i]]@summary[, 4]
best_fitness   <- GA[[i]]@summary[, 1]
# Vetor categórico com 3 níveis (mediana/média/melhor), repetido nº_gerações vezes cada.
Estatisticas   <- c(rep("Mediana", length(geracao)),
                    rep("Média",   length(geracao)),
                    rep("Melhor",  length(geracao)))
# Data.frame "long" para ggplot: 3 séries empilhadas em uma coluna 'Aptidao'.
data_grafico_1 <- data.frame(
  Geracao     = rep(geracao, 3),
  Aptidao     = c(mean_fitness, median_fitness, best_fitness),
  Estatisticas = Estatisticas
)

# 3 linhas (uma por estatística), com pontos sobre cada geração; eixo X com breaks inteiros.
ggplot(data_grafico_1, aes(x = Geracao, y = Aptidao, group = Estatisticas)) +
  geom_line(aes(colour = Estatisticas, linetype = Estatisticas), size = 2) +
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Geracao),
                                  max(data_grafico_1$Geracao), by = 1))
dev.off()



########## KERNEL RADIAL GAMMA = 0.001 #########
# Mesma estrutura da trilha anterior (LINEAR), agora com kernel RBF e gamma=0.001.
# Reaproveita 'rank_RF' já computado uma vez (RF não depende do kernel SVM).
# Veja o bloco LINEAR para comentários linha-a-linha; aqui anotamos só o essencial.
i <- 2                      # índice da trilha


# ---- Timer kernel: Radial_0.001 ----
t_k_Radial_0_001 <- proc.time()

## ---------------- SVM (classificação) e k-folds na trilha do rank RF ----------------
gamma  <- 0.001             # largura inversa do kernel gaussiano (quanto menor, mais "suave")
cost   <- 1.0
folds  <- 10
kernel <- "radial"          # RBF: K(x,y) = exp(-gamma * ||x - y||^2)

mean_svm_RF_list[[i]] <- numeric()
passo   <- 10
limite  <- floor((length(names(dados[[1]][, -ncol(dados[[1]])])) / passo) * percentual_snps)

t_corte_ini <- proc.time()  # inicio corte [Radial_0.001]
for (cont in 1:(limite + 1)) {
  cat("Grupo:", cont, "\n")
  j <- cont * passo
  if (j == 1) {
    var_sel <- names(rank_RF)[1]                # caso especial: 1 SNP só
  } else {
    var_sel <- head(names(rank_RF), j)          # top-j SNPs (head evita NAs se j > length)
  }
  base_cv <- dados[[1]][c(var_sel, "fenotipo")]
  # CV estratificada + SMOTE no treino com o kernel/gamma desta trilha.
  svm_cv <- validacao_cruzada_cls(
    data   = base_cv,
    folds  = folds,
    gamma  = gamma,
    cost   = cost,
    kernel = kernel
  )
  f1_mean <- svm_cv["f1_mean"]                   # F1 médio da CV
  mean_svm_RF_list[[i]][cont] <- f1_mean
}

## ---------------- Gráfico (Acurácia do SVM na trilha da RF) ----------------
pdf(file = "Grafico_F1_SVM_radial_0001.pdf", height = 5, width = 9)
plot(seq(passo, limite * passo + 10, by = passo),
     mean_svm_RF_list[[i]],
     type = "o", lwd = 2,
     xlab = "Grupo de Marcadores",
     ylab = "F1 (SVM Radial 0.001, 10-fold)")
maximo_idx  <- which.max(mean_svm_RF_list[[i]])   # nº ótimo de SNPs nesta trilha
buffer_alvo <- (maximo_idx + 1) * passo            # +1 passo de buffer p/ o AG refinar
corte[[i]]  <- min(buffer_alvo, length(rank_RF))
saturou     <- buffer_alvo > length(rank_RF)
abline(v   = corte[[i]],
       col = if (saturou) "red" else "black",
       lty = 2,
       lwd = if (saturou) 2     else 1)
if (saturou) {
  legend("topleft",
         legend  = sprintf("corte saturado em %d (buffer pediu %d)",
                           length(rank_RF), buffer_alvo),
         bty = "n", text.col = "red", cex = 0.8)
}
dev.off()

snps_selec_corte <- head(names(rank_RF), corte[[i]]) # top-N SNPs definidos pelo corte
snps_selec_corte

t_corte[["Radial_0.001"]] <- round((proc.time() - t_corte_ini)[["elapsed"]])
cat(sprintf("[TIMER] Corte Radial_0.001: %8.1f s\n", t_corte[["Radial_0.001"]]))
## ---------------- Seleciona genótipo e monta dados[[2]] após corte ----------------
genotipo <- list()
genotipo[[2]] <- dados[[1]][, head(names(rank_RF), corte[[i]]), drop = FALSE]  # só SNPs do corte
dados[[2]] <- cbind(genotipo[[2]], fenotipo = dados[[1]]$fenotipo)         # + alvo
if (!is.factor(dados[[2]]$fenotipo)) dados[[2]]$fenotipo <- factor(dados[[2]]$fenotipo, levels = c(0,1))

## ---------------- AG para REFINAMENTO (fitness = acurácia do SVM) ----------------
# 'f' é redefinida aqui propositalmente para capturar o kernel/gamma/cost ATUAIS via
# closure — sem essa redefinição o GA usaria os parâmetros da trilha anterior.
f <- function(x) {
  inc <- which(x == 1)                     # bits ativos = SNPs selecionados pelo cromossomo
  if (length(inc) == 0) return(0)
  dados_validacao <- dados[[2]][, c(inc, ncol(dados[[2]])), drop = FALSE]
  colnames(dados_validacao)[ncol(dados_validacao)] <- "fenotipo"
  res <- validacao_cruzada_cls(
    data   = dados_validacao,
    folds  = folds,
    gamma  = gamma,
    cost   = cost,
    kernel = kernel
  )
  f1_mean <- res["f1_mean"]
  return(f1_mean)
}
fitness <- function(x) f(x)

# Hiperparâmetros do GA (mesmos da trilha linear).
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

sol_mat   <- as.matrix(GA[[i]]@solution)
best_row  <- 1                                # 1ª solução em caso de empate
sel_cols  <- which(sol_mat[best_row, ] == 1)
snps_selec_ref[[i]] <- colnames(sol_mat)[sel_cols]
t_ga[["Radial_0.001"]] <- round((proc.time() - t_ga_ini)[["elapsed"]])
cat(sprintf("[TIMER] GA    Radial_0.001: %8.1f s\n", t_ga[["Radial_0.001"]]))


# Alias amigável para o relatório final desta trilha.
snps_selec_SVM_radial_0001_GA_acc <- snps_selec_ref[[i]]
t_kernels[["Radial_0.001"]] <- round((proc.time() - t_k_Radial_0_001)[["elapsed"]])
snps_selec_SVM_radial_0001_GA_acc

# Resumo textual + gráfico básico do GA.
summary(GA[[i]])
plot(GA[[i]])

# Gráfico ggplot2 da evolução do fitness do GA (mediana/média/melhor por geração).
pdf(file = "Grafico_GA_SVM_radial_0001_F1.pdf", height = 5, width = 9)
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
# Idêntico à trilha anterior (gamma=0.001), com gamma=0.01 (RBF mais "agudo").
# Comentários linha-a-linha estão no bloco LINEAR; aqui anotamos só o essencial.
i <- 3                      # índice da trilha


# ---- Timer kernel: Radial_0.01 ----
t_k_Radial_0_01 <- proc.time()

## ---------------- SVM (classificação) e k-folds na trilha do rank RF ----------------
gamma  <- 0.01              # gamma maior -> kernel mais sensível a pontos próximos
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
    var_sel <- head(names(rank_RF), j)
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
  
  f1_mean <- svm_cv["f1_mean"]  # F1 médio da CV
  mean_svm_RF_list[[i]][cont] <- f1_mean
}

## ---------------- Gráfico (Acurácia do SVM na trilha da RF) ----------------
pdf(file = "Grafico_F1_SVM_radial_001.pdf", height = 5, width = 9)
plot(seq(passo, limite * passo + 10, by = passo),
     mean_svm_RF_list[[i]],
     type = "o", lwd = 2,
     xlab = "Grupo de Marcadores",
     ylab = "F1 (SVM Radial 0.01, 10-fold)")
# Escolha do ponto de corte: agora é o MÁXIMO de acurácia
maximo_idx  <- which.max(mean_svm_RF_list[[i]])
buffer_alvo <- (maximo_idx + 1) * passo            # +1 passo de buffer p/ o AG refinar
corte[[i]]  <- min(buffer_alvo, length(rank_RF))
saturou     <- buffer_alvo > length(rank_RF)
abline(v   = corte[[i]],
       col = if (saturou) "red" else "black",
       lty = 2,
       lwd = if (saturou) 2     else 1)
if (saturou) {
  legend("topleft",
         legend  = sprintf("corte saturado em %d (buffer pediu %d)",
                           length(rank_RF), buffer_alvo),
         bty = "n", text.col = "red", cex = 0.8)
}
dev.off()

snps_selec_corte <- head(names(rank_RF), corte[[i]])
snps_selec_corte

t_corte[["Radial_0.01"]] <- round((proc.time() - t_corte_ini)[["elapsed"]])
cat(sprintf("[TIMER] Corte Radial_0.01 : %8.1f s\n", t_corte[["Radial_0.01"]]))
## ---------------- Seleciona genótipo e monta dados[[2]] após corte ----------------
genotipo <- list()
genotipo[[2]] <- dados[[1]][, head(names(rank_RF), corte[[i]]), drop = FALSE]

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
  f1_mean <- res["f1_mean"]  # F1 médio do conjunto candidato
  return(f1_mean)             # GA maximiza F1 em vez de acurácia
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
pdf(file = "Grafico_GA_SVM_radial_001_F1.pdf", height = 5, width = 9)
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
# Idêntico às trilhas anteriores, agora com gamma=0.1. Veja LINEAR para detalhes.
i <- 4                      # índice da trilha


# ---- Timer kernel: Radial_0.1 ----
t_k_Radial_0_1 <- proc.time()

## ---------------- SVM (classificação) e k-folds na trilha do rank RF ----------------
gamma  <- 0.1               # gamma intermediário (entre 0.01 e 1)
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
    var_sel <- head(names(rank_RF), j)
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
  
  f1_mean <- svm_cv["f1_mean"]  # F1 médio da CV
  mean_svm_RF_list[[i]][cont] <- f1_mean
}

## ---------------- Gráfico (Acurácia do SVM na trilha da RF) ----------------
pdf(file = "Grafico_F1_SVM_radial_01.pdf", height = 5, width = 9)
plot(seq(passo, limite * passo + 10, by = passo),
     mean_svm_RF_list[[i]],
     type = "o", lwd = 2,
     xlab = "Grupo de Marcadores",
     ylab = "F1 (SVM Radial 0.1, 10-fold)")
# Escolha do ponto de corte: MÁXIMO do F1
maximo_idx  <- which.max(mean_svm_RF_list[[i]])
buffer_alvo <- (maximo_idx + 1) * passo            # +1 passo de buffer p/ o AG refinar
corte[[i]]  <- min(buffer_alvo, length(rank_RF))
saturou     <- buffer_alvo > length(rank_RF)
abline(v   = corte[[i]],
       col = if (saturou) "red" else "black",
       lty = 2,
       lwd = if (saturou) 2     else 1)
if (saturou) {
  legend("topleft",
         legend  = sprintf("corte saturado em %d (buffer pediu %d)",
                           length(rank_RF), buffer_alvo),
         bty = "n", text.col = "red", cex = 0.8)
}
dev.off()

snps_selec_corte <- head(names(rank_RF), corte[[i]])
snps_selec_corte

t_corte[["Radial_0.1"]] <- round((proc.time() - t_corte_ini)[["elapsed"]])
cat(sprintf("[TIMER] Corte Radial_0.1  : %8.1f s\n", t_corte[["Radial_0.1"]]))
## ---------------- Seleciona genótipo e monta dados[[2]] após corte ----------------
genotipo <- list()
genotipo[[2]] <- dados[[1]][, head(names(rank_RF), corte[[i]]), drop = FALSE]

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
  f1_mean <- res["f1_mean"]  # F1 médio do conjunto candidato
  return(f1_mean)             # GA maximiza F1 em vez de acurácia
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
pdf(file = "Grafico_GA_SVM_radial_01_F1.pdf", height = 5, width = 9)
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
# Última trilha SVM (gamma=1.0). Veja LINEAR para os comentários linha-a-linha.
i <- 5                      # índice da trilha


# ---- Timer kernel: Radial_1.0 ----
t_k_Radial_1_0 <- proc.time()

## ---------------- SVM (classificação) e k-folds na trilha do rank RF ----------------
gamma  <- 1.0               # kernel bastante "agudo" — ajuste mais local, alto risco de overfit
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
    var_sel <- head(names(rank_RF), j)
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
  
  f1_mean <- svm_cv["f1_mean"]  # F1 médio da CV
  mean_svm_RF_list[[i]][cont] <- f1_mean
}

## ---------------- Gráfico (Acurácia do SVM na trilha da RF) ----------------
pdf(file = "Grafico_F1_SVM_radial_1.pdf", height = 5, width = 9)
plot(seq(passo, limite * passo + 10, by = passo),
     mean_svm_RF_list[[i]],
     type = "o", lwd = 2,
     xlab = "Grupo de Marcadores",
     ylab = "F1 (SVM Radial 1.0, 10-fold)")
# Escolha do ponto de corte: MÁXIMO do F1
maximo_idx  <- which.max(mean_svm_RF_list[[i]])
buffer_alvo <- (maximo_idx + 1) * passo            # +1 passo de buffer p/ o AG refinar
corte[[i]]  <- min(buffer_alvo, length(rank_RF))
saturou     <- buffer_alvo > length(rank_RF)
abline(v   = corte[[i]],
       col = if (saturou) "red" else "black",
       lty = 2,
       lwd = if (saturou) 2     else 1)
if (saturou) {
  legend("topleft",
         legend  = sprintf("corte saturado em %d (buffer pediu %d)",
                           length(rank_RF), buffer_alvo),
         bty = "n", text.col = "red", cex = 0.8)
}
dev.off()

snps_selec_corte <- head(names(rank_RF), corte[[i]])
snps_selec_corte

t_corte[["Radial_1.0"]] <- round((proc.time() - t_corte_ini)[["elapsed"]])
cat(sprintf("[TIMER] Corte Radial_1.0  : %8.1f s\n", t_corte[["Radial_1.0"]]))
## ---------------- Seleciona genótipo e monta dados[[2]] após corte ----------------
genotipo <- list()
genotipo[[2]] <- dados[[1]][, head(names(rank_RF), corte[[i]]), drop = FALSE]

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
  f1_mean <- res["f1_mean"]  # F1 médio do conjunto candidato
  return(f1_mean)             # GA maximiza F1 em vez de acurácia
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
pdf(file = "Grafico_GA_SVM_radial_1_F1.pdf", height = 5, width = 9)
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
# Fim das 5 trilhas SVM (linear + 4 radiais). A partir daqui consolidamos
# as seleções (união/interseção) e geramos relatório.


####Uniao dos SNPs selecionados####
# 'lista_SNPs' precisa ser uma LISTA com >= 2 elementos (cada elemento = vetor
# de SNPs selecionados por uma trilha). Faz a união cumulativa em cadeia:
# union(L[[1]], L[[2]]) -> depois com L[[3]] -> depois com L[[4]] -> ...
uniao_snps<-function(lista_SNPs){
  uniao <- list()
  # Caso base: união dos dois primeiros elementos.
  uniao[[1]]<-union(lista_SNPs[[1]],lista_SNPs[[2]])
  # Passos seguintes: vai acumulando a união com o próximo elemento da lista.
  for (a in 1:(length(lista_SNPs)-2)){
    uniao[[a+1]]<-union(uniao[[a]],lista_SNPs[[a+2]])
  }
  # 'a' permanece após o loop (R não tem escopo de bloco); aqui é o último índice usado.
  uniao_final <- uniao[[a+1]]
  return(uniao_final)
}

# Chamada inicial só para imprimir no console (não captura).
uniao_snps(snps_selec_ref)

# União final entre as seleções das 5 trilhas (1 vetor com SNPs únicos).
uniao_final<-uniao_snps(snps_selec_ref)

####Intersecao dos SNPs selecionados####
# Mesma ideia da união, mas com 'intersect' (interseção cumulativa).
# Se a interseção for vazia, devolve a string "Conjunto vazio" para o relatório.
intersecao_snps<-function(lista_SNPs){
  intersecao <- list()
  intersecao[[1]]<-intersect(lista_SNPs[[1]],lista_SNPs[[2]])
  for (a in 1:(length(lista_SNPs)-2)){
    intersecao[[a+1]]<-intersect(intersecao[[a]],lista_SNPs[[a+2]])
  }
  intersecao_final <- intersecao[[a+1]]
  # Mensagem amigável quando não há nenhum SNP comum a todas as trilhas.
  if (length(intersecao_final) == 0L) {
    return("Conjunto vazio")
  } else
    return(intersecao_final)
}

intersecao_snps(snps_selec_ref)

# Interseção final (SNPs presentes em TODAS as 5 trilhas).
intersecao_final<-intersecao_snps(snps_selec_ref)

###############################################
## ACRESCIMOS: Padronização de saídas (CLASSIFICAÇÃO)
## NÃO altera nada do SMS original; apenas acrescenta.
## Gera CSVs/TXT consolidados (rank RF, valor-p, seleção por kernel, métricas).
###############################################

## 1) Preparos: importância da RF como data.frame + valor-p binário
# 'imp' e 'mdg_col' foram definidos lá em cima (uma única RF compartilhada).
# Convertemos a coluna escolhida (MeanDecreaseGini) em data.frame para reuso/escrita.
rank_RF_df <- data.frame(Importancia_RF = imp[, mdg_col, drop = TRUE])
# Salva ranking completo (todos os SNPs com sua importância) em arquivo de texto.
write.csv(rank_RF_df, "rank_Random_Forest_CLASSIF_F1.txt", row.names = TRUE)

# Reaproveita 'pvals' já calculado, se existir; do contrário, recalcula da base completa.
if (!exists("pvals")) {
  pvals <- valor.p.binario(dados[[1]])
}

# Seleções por valor-p para o relatório: bruto e ajustado por Bonferroni.
# 'which(... <= 0.05)' devolve índices significativos; 'rownames' converte em nomes de SNP.
pvals_bruto_selec    <- rownames(pvals)[which(pvals[["Valor.p.bruto"]]    <= 0.05)]
pvals_ajustado_selec <- rownames(pvals)[which(pvals[["Valor.p.ajustado"]] <= 0.05)]

## 2) Consolida as listas de SNPs selecionados por kernel (já produzidas nas 5 trilhas)
# Para cada trilha, usamos 'if(exists(...)) ... else character(0)' como salvaguarda:
# se a trilha não tiver sido executada (ex.: rodando o script em pedaços), entra vetor vazio.
selec_por_kernel <- list(
  "Linear"       = if (exists("snps_selec_SVM_linear_GA_acc"))        snps_selec_SVM_linear_GA_acc        else character(0),
  "Radial_0.001" = if (exists("snps_selec_SVM_radial_0001_GA_acc")) snps_selec_SVM_radial_0001_GA_acc else character(0),
  "Radial_0.01"  = if (exists("snps_selec_SVM_radial_001_GA_acc"))  snps_selec_SVM_radial_001_GA_acc  else character(0),
  "Radial_0.1"   = if (exists("snps_selec_SVM_radial_01_GA_acc"))   snps_selec_SVM_radial_01_GA_acc   else character(0),
  "Radial_1"     = if (exists("snps_selec_SVM_radial_1_GA_acc"))    snps_selec_SVM_radial_1_GA_acc    else character(0)
)

## 3) Cria dataframes "rank_global_selec" (RF + valor-p) para cada seleção e salva
# Cria pasta de saída (silencia warning se ela já existir).
  dir.create("rank_global_classif_F1", showWarnings = FALSE)
# Lista que receberá o data.frame consolidado por kernel.
rank_global_selec_CLASSIF <- list()

# Para cada trilha (kernel), filtra rank_RF_df + pvals nos SNPs daquela seleção e grava.
for (lbl in names(selec_por_kernel)) {
  snps_sel <- selec_por_kernel[[lbl]]
  if (length(snps_sel) > 0) {                          # pula trilhas sem seleção
    rf_part  <- rank_RF_df[snps_sel, , drop = FALSE]   # importância dos SNPs selecionados
    pv_part  <- pvals[snps_sel, , drop = FALSE]        # p-valores dos SNPs selecionados
    df_out   <- cbind(rf_part, pv_part)                # une lado a lado (mesma ordem de linhas)
    rank_global_selec_CLASSIF[[lbl]] <- df_out
    write.csv(df_out,
              file = file.path("rank_global_classif_F1",
                               paste0("rank_global_CLASSIF_F1_", lbl, ".txt")),
              row.names = TRUE)
  }
}

## 4) Helpers de impressão para o relatório
# Devolve "Conjunto vazio" para vetor vazio ou marca textual, e
# "snp1, snp2, ..." (separado por vírgula) para o caso normal.
collapse_or_label <- function(x) {
  if (length(x) == 0L) return("Conjunto vazio")
  if (is.character(x) && length(x) == 1L && x == "Conjunto vazio") return(x)
  paste(x, collapse = ", ")
}

## 5) "Medidas de Otimalidade" dos SNPs causais (CLASSIFICAÇÃO)
##    Mede o desempenho ideal: rodando o SVM SÓ com os 8 SNPs causais conhecidos
##    da simulação (SNP1..SNP8). Serve como teto de referência por kernel.
##    Em dados reais, basta ajustar 'var_sel_causais' para o conjunto candidato.
var_sel_causais <- c("SNP1","SNP2","SNP3","SNP4","SNP5","SNP6","SNP7","SNP8")

# Extrai as 4 métricas médias (ACC/PREC/REC/F1) do vetor nomeado devolvido pela CV.
# 'unname' tira os names (deixa o vetor com 4 nomes próprios).
pega_metricas_cls <- function(res_vec) {
  # Ordem das posições no vetor de saída de validacao_cruzada_cls():
  # 1=acc_mean, 2=acc_sd, 3=prec_mean, 4=prec_sd, 5=rec_mean, 6=rec_sd, 7=f1_mean, 8=f1_sd
  out <- c(
    ACC   = unname(res_vec["acc_mean"]),
    PREC  = unname(res_vec["prec_mean"]),
    REC   = unname(res_vec["rec_mean"]),
    F1    = unname(res_vec["f1_mean"])
  )
  return(out)
}

# Lista que receberá as métricas (por kernel) calculadas só com os SNPs causais.
metricas_causais <- list()

# 'all(var_sel_causais %in% names(dados[[1]]))' = TRUE quando todos os SNPs
# nomeados existem como coluna do data.frame original.
if (all(var_sel_causais %in% names(dados[[1]]))) {
  # Sub-base só com os 8 causais + alvo.
  base_causais <- dados[[1]][c(var_sel_causais, "fenotipo")]

  # ----- Kernel Linear -----
  res_lin <- validacao_cruzada_cls(
    data = base_causais, folds = 10, gamma = 0.01, cost = 1, kernel = "linear"
  )
  metricas_causais[["Linear"]] <- pega_metricas_cls(res_lin)

  # ----- Radial gamma=0.001 -----
  res_r1 <- validacao_cruzada_cls(
    data = base_causais, folds = 10, gamma = 0.001, cost = 1, kernel = "radial"
  )
  metricas_causais[["Radial_0.001"]] <- pega_metricas_cls(res_r1)

  # ----- Radial gamma=0.01 -----
  res_r2 <- validacao_cruzada_cls(
    data = base_causais, folds = 10, gamma = 0.01, cost = 1, kernel = "radial"
  )
  metricas_causais[["Radial_0.01"]] <- pega_metricas_cls(res_r2)

  # ----- Radial gamma=0.1 -----
  res_r3 <- validacao_cruzada_cls(
    data = base_causais, folds = 10, gamma = 0.1, cost = 1, kernel = "radial"
  )
  metricas_causais[["Radial_0.1"]] <- pega_metricas_cls(res_r3)

  # ----- Radial gamma=1 -----
  res_r4 <- validacao_cruzada_cls(
    data = base_causais, folds = 10, gamma = 1.0, cost = 1, kernel = "radial"
  )
  metricas_causais[["Radial_1"]] <- pega_metricas_cls(res_r4)
}

## 6) Relatório consolidado (texto único)
# 'sink' redireciona TODA a saída de 'cat' para o arquivo até o próximo 'sink()'.
# ---- Calcula tempo total ----
t_total_elapsed <- round((proc.time() - t_total_inicio)[["elapsed"]])

sink("my_output_SMS_CLASSIF_F1_2.txt")

# Seção 1: seleções do GA por kernel.
cat("===== Seleções por Kernel (GA, Classificação) =====\n")
cat("SMS Linear = ",          collapse_or_label(selec_por_kernel[["Linear"]]),        "\n")
cat("SMS Radial 0.001 = ",    collapse_or_label(selec_por_kernel[["Radial_0.001"]]), "\n")
cat("SMS Radial 0.01 = ",     collapse_or_label(selec_por_kernel[["Radial_0.01"]]),  "\n")
cat("SMS Radial 0.1 = ",      collapse_or_label(selec_por_kernel[["Radial_0.1"]]),   "\n")
cat("SMS Radial 1 = ",        collapse_or_label(selec_por_kernel[["Radial_1"]]),     "\n\n")

# Seção 2: união e interseção das seleções (computadas em uniao_snps/intersecao_snps).
cat("===== União / Interseção =====\n")
cat("União = ",       collapse_or_label(uniao_final),       "\n")
cat("Interseção = ",  collapse_or_label(intersecao_final),  "\n\n")

# Seção 3: seleções por significância estatística direta (sem ML).
cat("===== Seleção por Valor-p (binário) =====\n")
cat("Valor-p bruto (<= 0.05) = ",    collapse_or_label(pvals_bruto_selec),    "\n")
cat("Valor-p ajustado (<= 0.05) = ", collapse_or_label(pvals_ajustado_selec), "\n\n")

# Seção 4: métricas com os SNPs causais "verdadeiros" (referência de teto).
cat("===== Medidas de Otimalidade (SNPs causais, Classificação) =====\n")
if (length(metricas_causais) == 0L) {
  cat("Conjunto 'var_sel_causais' não encontrado nas colunas. Ajuste os nomes e reexecute.\n")
} else {
  for (lbl in names(metricas_causais)) {
    m <- metricas_causais[[lbl]]
    # sprintf("%.4f", ...) -> formata com 4 casas decimais para alinhar a impressão.
    cat(lbl, " ->  ACC =", sprintf("%.4f", m["ACC"]),
        " | PREC =", sprintf("%.4f", m["PREC"]),
        " | REC  =", sprintf("%.4f", m["REC"]),
        " | F1   =", sprintf("%.4f", m["F1"]), "\n")
  }
}

# Fecha o redirecionamento (a próxima saída volta para o console).

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
