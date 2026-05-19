# Script para salvar os dados balanceados após SMOTE e gerar resumo
# Este script conta as classes diretamente dos arquivos CSV:
# - dados.csv (ou dados_original_antes_SMOTE.csv) para dados ANTES do SMOTE
# - dados_balanceado_SMOTE.csv para dados DEPOIS do SMOTE

# Definindo a trilha de dados
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cat("\n========================================\n")
cat("GERANDO RESUMO DO BALANCEAMENTO COM SMOTE\n")
cat("========================================\n\n")

#------------------------------------------------------------
# Função auxiliar para carregar e processar dados.csv
#------------------------------------------------------------
carregar_dados_original <- function(arquivo) {
  cat("Carregando", arquivo, "...\n")
  dados_temp <- read.csv(arquivo, stringsAsFactors = FALSE)
  
  # Verifica se dados.csv foi salvo como lista (res_exact do scrime)
  # Nesse caso, as colunas começam com "data." (ex: "data.SNP1", "data.fenotipo")
  if (any(grepl("^data\\.", names(dados_temp)))) {
    cat("  Detectado formato de lista (res_exact). Extraindo colunas 'data.*'...\n")
    # Extrai apenas as colunas que começam com "data."
    colunas_data <- grep("^data\\.", names(dados_temp), value = TRUE)
    dados_original <- dados_temp[, colunas_data, drop = FALSE]
    
    # Remove o prefixo "data." dos nomes das colunas
    names(dados_original) <- sub("^data\\.", "", names(dados_original))
    
    cat("  Colunas extraídas:", length(colunas_data), "\n")
  } else {
    # Formato normal de data.frame
    dados_original <- dados_temp
  }
  
  # Remove colunas extras se existirem
  if ("linpred" %in% names(dados_original)) {
    dados_original$linpred <- NULL
  }
  if ("prob" %in% names(dados_original)) {
    dados_original$prob <- NULL
  }
  
  # Verifica se fenotipo existe
  if (!"fenotipo" %in% names(dados_original)) {
    stop("Coluna 'fenotipo' não encontrada em ", arquivo, " após processamento.")
  }
  
  # Converte fenotipo para numérico se necessário
  if (is.factor(dados_original$fenotipo)) {
    dados_original$fenotipo <- as.numeric(as.character(dados_original$fenotipo))
  } else if (is.character(dados_original$fenotipo)) {
    dados_original$fenotipo <- as.numeric(dados_original$fenotipo)
  }
  
  cat("  Total de observações:", nrow(dados_original), "\n")
  cat("  Distribuição de classes:\n")
  print(table(dados_original$fenotipo))
  cat("\n")
  
  return(dados_original)
}

#------------------------------------------------------------
# Carregando dados ANTES do SMOTE
#------------------------------------------------------------
cat("1. CARREGANDO DADOS ANTES DO SMOTE\n")
cat("   ----------------------------------------\n")

# Prioridade: dados_original_antes_SMOTE.csv > dados.csv
if (file.exists("dados_original_antes_SMOTE.csv")) {
  dados_original <- carregar_dados_original("dados_original_antes_SMOTE.csv")
} else if (file.exists("dados.csv")) {
  dados_original <- carregar_dados_original("dados.csv")
} else {
  stop("ERRO: Nenhum arquivo de dados originais encontrado.\n",
       "  Procurados: dados_original_antes_SMOTE.csv, dados.csv")
}

# Verifica se há dados válidos
if (nrow(dados_original) == 0) {
  stop("ERRO: dados_original está vazio. Verifique o arquivo de origem.")
}

#------------------------------------------------------------
# Carregando dados DEPOIS do SMOTE
#------------------------------------------------------------
cat("2. CARREGANDO DADOS DEPOIS DO SMOTE\n")
cat("   ----------------------------------------\n")

if (!file.exists("dados_balanceado_SMOTE.csv")) {
  stop("ERRO: Arquivo 'dados_balanceado_SMOTE.csv' não encontrado.\n",
       "  Execute o script principal primeiro para gerar este arquivo.")
}

cat("Carregando dados_balanceado_SMOTE.csv...\n")
dados_balanceados <- read.csv("dados_balanceado_SMOTE.csv", stringsAsFactors = FALSE)

# Verifica se fenotipo existe
if (!"fenotipo" %in% names(dados_balanceados)) {
  stop("Coluna 'fenotipo' não encontrada em dados_balanceado_SMOTE.csv")
}

# Converte fenotipo para numérico se necessário
if (is.factor(dados_balanceados$fenotipo)) {
  dados_balanceados$fenotipo <- as.numeric(as.character(dados_balanceados$fenotipo))
} else if (is.character(dados_balanceados$fenotipo)) {
  dados_balanceados$fenotipo <- as.numeric(dados_balanceados$fenotipo)
}

cat("  Total de observações:", nrow(dados_balanceados), "\n")
cat("  Distribuição de classes:\n")
print(table(dados_balanceados$fenotipo))
cat("\n")

#------------------------------------------------------------
# Gerando resumo das contagens
#------------------------------------------------------------
cat("3. GERANDO RESUMO\n")
cat("   ----------------------------------------\n")

# Cria tabelas de distribuição
tab_antes <- table(dados_original$fenotipo)
tab_depois <- table(dados_balanceados$fenotipo)

cat("Tabela ANTES (dados originais):\n")
print(tab_antes)
cat("\nTabela DEPOIS (dados balanceados):\n")
print(tab_depois)
cat("\n")

# Função auxiliar para obter contagem de classe (trata NA)
get_count <- function(tab, classe) {
  count <- as.numeric(tab[as.character(classe)])
  if (is.na(count)) return(0)
  return(count)
}

# Verifica se as tabelas têm dados válidos
if (length(tab_antes) == 0) {
  stop("ERRO: tab_antes está vazia. Verifique se dados_original foi carregado corretamente.")
}
if (length(tab_depois) == 0) {
  stop("ERRO: tab_depois está vazia. Verifique se dados_balanceados foi carregado corretamente.")
}

# Obtém os valores das classes presentes nas tabelas
classes_antes <- as.numeric(names(tab_antes))
classes_depois <- as.numeric(names(tab_depois))
todas_classes <- sort(unique(c(classes_antes, classes_depois)))

# Garante que temos pelo menos as classes 0 e 1
if (!0 %in% todas_classes) todas_classes <- c(0, todas_classes)
if (!1 %in% todas_classes) todas_classes <- c(todas_classes, 1)
todas_classes <- sort(unique(todas_classes))

# Cria vetores para o resumo
classe_vec <- c(as.character(todas_classes), "Total")
antes_vec <- c(sapply(todas_classes, function(c) get_count(tab_antes, c)), sum(tab_antes))
depois_vec <- c(sapply(todas_classes, function(c) get_count(tab_depois, c)), sum(tab_depois))
diff_vec <- depois_vec - antes_vec

resumo_smote <- data.frame(
  Classe = classe_vec,
  Antes_SMOTE = antes_vec,
  Depois_SMOTE = depois_vec,
  Diferenca = diff_vec
)

# Salva o resumo em CSV
write.csv(resumo_smote, file = "resumo_SMOTE.csv", row.names = FALSE)
cat("Arquivo 'resumo_SMOTE.csv' salvo com sucesso.\n\n")

#------------------------------------------------------------
# Exibindo resumo no console
#------------------------------------------------------------
cat("========================================\n")
cat("RESUMO DO BALANCEAMENTO COM SMOTE\n")
cat("========================================\n")
cat("\nDistribuição ANTES do SMOTE:\n")
cat("  Classe 0:", get_count(tab_antes, 0), "observações\n")
cat("  Classe 1:", get_count(tab_antes, 1), "observações\n")
cat("  Total   :", sum(tab_antes), "observações\n")
cat("\nDistribuição DEPOIS do SMOTE:\n")
cat("  Classe 0:", get_count(tab_depois, 0), "observações\n")
cat("  Classe 1:", get_count(tab_depois, 1), "observações\n")
cat("  Total   :", sum(tab_depois), "observações\n")
cat("\nAmostras sintéticas criadas:\n")
cat("  Classe 0:", get_count(tab_depois, 0) - get_count(tab_antes, 0), "observações\n")
cat("  Classe 1:", get_count(tab_depois, 1) - get_count(tab_antes, 1), "observações\n")
cat("  Total   :", sum(tab_depois) - sum(tab_antes), "observações\n")
cat("========================================\n\n")

cat("Processo concluído com sucesso!\n")
