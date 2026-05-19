# ==============================================================================
#  Busca contagem REAL de publicações do NHGRI-EBI GWAS Catalog por ano
#  Fonte: API v2 — https://www.ebi.ac.uk/gwas/rest/api/v2/publications
# ==============================================================================

library(httr)
library(jsonlite)

BASE_URL  <- "https://www.ebi.ac.uk/gwas/rest/api/v2/publications"
PAGE_SIZE <- 500
SSL_OPTS  <- config(ssl_verifypeer = FALSE)

fetch_json <- function(url) {
  resp <- GET(url, SSL_OPTS)
  if (status_code(resp) != 200) stop("HTTP ", status_code(resp), " em: ", url)
  fromJSON(rawToChar(resp$content), simplifyVector = FALSE)
}

cat("Consultando GWAS Catalog API v2...\n")

# ---- 1. Descobrir total de publicações ----------------------------------------
first <- fetch_json(paste0(BASE_URL, "?size=1&page=0"))
total       <- first$page$totalElements
total_pages <- ceiling(total / PAGE_SIZE)
cat(sprintf("Total de publicações: %d | Páginas a buscar: %d\n",
            total, total_pages))

# ---- 2. Baixar todas as páginas -----------------------------------------------
all_dates <- character(0)

for (pg in 0:(total_pages - 1)) {
  url <- paste0(BASE_URL, "?size=", PAGE_SIZE, "&page=", pg)
  cat(sprintf("  Página %02d / %d ...", pg + 1, total_pages))

  resp <- tryCatch(
    fetch_json(url),
    error = function(e) {
      cat(" ERRO:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(resp)) next

  pubs  <- resp[["_embedded"]][["publications"]]
  dates <- sapply(pubs, function(p) p[["publication_date"]])
  all_dates <- c(all_dates, dates)
  cat(sprintf(" %d registros (total acum.: %d)\n", length(dates), length(all_dates)))

  Sys.sleep(0.4)
}

cat(sprintf("\nDatas coletadas: %d de %d\n", length(all_dates), total))

# ---- 3. Contar publicações por ano --------------------------------------------
anos     <- as.integer(substr(all_dates, 1, 4))
contagem <- table(anos)
df <- data.frame(
  ano    = as.integer(names(contagem)),
  anuais = as.integer(contagem),
  stringsAsFactors = FALSE
)
df <- df[order(df$ano), ]
df$acumulado <- cumsum(df$anuais)

cat("\nPublicações por ano (GWAS Catalog — dados reais da API):\n")
print(df, row.names = FALSE)

# ---- 4. Salvar CSV -----------------------------------------------------------
OUT_CSV <- "gwas_publicacoes_por_ano.csv"
write.csv(df, OUT_CSV, row.names = FALSE)
cat(sprintf("\nDados salvos em: %s\n", OUT_CSV))
