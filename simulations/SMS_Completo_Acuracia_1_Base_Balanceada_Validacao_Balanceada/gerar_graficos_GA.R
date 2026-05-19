# ============================================
# Script para gerar apenas os gráficos GA
# a partir dos dados salvos em .RData
# ============================================

# Definindo a trilha de dados
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Carregando pacotes necessários
library(ggplot2)

# Carregando os dados salvos
load("SMS_acuracia_1.RData")

# Verificando se os objetos GA existem
if (!exists("GA")) {
  stop("Objeto 'GA' não encontrado no arquivo .RData")
}

cat("Gerando gráficos GA a partir dos dados salvos...\n")

# Função auxiliar para gerar gráfico GA
gerar_grafico_GA <- function(ga_obj, nome_arquivo, i) {
  cat(sprintf("Gerando %s...\n", nome_arquivo))
  
  pdf(file = nome_arquivo, height = 5, width = 9)
  geracao        <- seq_len(ga_obj@iter)
  mean_fitness   <- ga_obj@summary[, 2]
  median_fitness <- ga_obj@summary[, 4]
  best_fitness   <- ga_obj@summary[, 1]
  Estatisticas   <- c(rep("Mediana", length(geracao)),
                      rep("Média",   length(geracao)),
                      rep("Melhor",  length(geracao)))
  data_grafico_1 <- data.frame(
    Geracao     = rep(geracao, 3),
    Aptidao     = c(mean_fitness, median_fitness, best_fitness),
    Estatisticas = Estatisticas
  )
  
  p <- ggplot(data_grafico_1, aes(x = Geracao, y = Aptidao, group = Estatisticas)) +
    geom_line(aes(colour = Estatisticas, linetype = Estatisticas), size = 2) +
    geom_point() +
    scale_x_continuous(breaks = seq(min(data_grafico_1$Geracao),
                                    max(data_grafico_1$Geracao), by = 1))
  print(p)
  dev.off()
  
  cat(sprintf("  ✓ %s gerado com sucesso!\n", nome_arquivo))
}

# Gerando os 5 gráficos GA
# i = 1: Linear
if (length(GA) >= 1 && !is.null(GA[[1]])) {
  gerar_grafico_GA(GA[[1]], "Grafico_GA_SVM_linear_ACC.pdf", 1)
} else {
  cat("Aviso: GA[[1]] não encontrado. Pulando gráfico linear.\n")
}

# i = 2: Radial gamma = 0.001
if (length(GA) >= 2 && !is.null(GA[[2]])) {
  gerar_grafico_GA(GA[[2]], "Grafico_GA_SVM_radial_0001_ACC.pdf", 2)
} else {
  cat("Aviso: GA[[2]] não encontrado. Pulando gráfico radial 0.001.\n")
}

# i = 3: Radial gamma = 0.01
if (length(GA) >= 3 && !is.null(GA[[3]])) {
  gerar_grafico_GA(GA[[3]], "Grafico_GA_SVM_radial_001_ACC.pdf", 3)
} else {
  cat("Aviso: GA[[3]] não encontrado. Pulando gráfico radial 0.01.\n")
}

# i = 4: Radial gamma = 0.1
if (length(GA) >= 4 && !is.null(GA[[4]])) {
  gerar_grafico_GA(GA[[4]], "Grafico_GA_SVM_radial_01_ACC.pdf", 4)
} else {
  cat("Aviso: GA[[4]] não encontrado. Pulando gráfico radial 0.1.\n")
}

# i = 5: Radial gamma = 1.0
if (length(GA) >= 5 && !is.null(GA[[5]])) {
  gerar_grafico_GA(GA[[5]], "Grafico_GA_SVM_radial_1_ACC.pdf", 5)
} else {
  cat("Aviso: GA[[5]] não encontrado. Pulando gráfico radial 1.0.\n")
}

cat("\n✓ Todos os gráficos foram gerados!\n")
