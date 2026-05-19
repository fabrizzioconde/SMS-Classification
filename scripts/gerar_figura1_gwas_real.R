# ==============================================================================
#  Figura 1 — Publicações GWAS Catalog (2006–2025)
#  Barras: publicações ACUMULADAS (eixo esquerdo)
#  Linha:  publicações ANUAIS     (eixo direito)
# ==============================================================================

library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Dados reais da API
# ------------------------------------------------------------------------------
dados <- data.frame(
  ano    = 2006:2025,
  anuais = c(8, 89, 146, 233, 331, 398, 393, 392,
             348, 369, 356, 444, 554, 566, 536, 691,
             622, 554, 369, 220)
)
dados$acumulado <- cumsum(dados$anuais) + 10   # +10 dos 2 artigos de 2005

# Fator de escala: o eixo primário (barras) usa a escala do acumulado;
# a linha (anuais) é escalonada para caber nesse mesmo eixo.
escala <- max(dados$acumulado) / max(dados$anuais)

# ------------------------------------------------------------------------------
# 2. Gráfico combinado: barras (acumulado) + linha (anuais)
# ------------------------------------------------------------------------------
COR_BARRA <- "#1B5E20"   # verde-escuro  — acumulado (barras)
COR_LINHA <- "#B71C1C"   # vermelho-escuro — anuais (linha)

p <- ggplot(dados, aes(x = factor(ano))) +

  # Barras — publicações ACUMULADAS (eixo esquerdo)
  geom_col(aes(y = acumulado), fill = COR_BARRA, alpha = 0.75,
           width = 0.70, colour = NA) +

  # Linha + pontos — publicações ANUAIS (eixo direito, escalonado)
  geom_line(aes(y = anuais * escala, group = 1),
            colour = COR_LINHA, linewidth = 1.0) +
  geom_point(aes(y = anuais * escala),
             colour = COR_LINHA, size = 2.0, shape = 19) +

  # Eixo Y esquerdo: acumulado
  scale_y_continuous(
    name   = "Publicações acumuladas",
    limits = c(0, max(dados$acumulado) * 1.08),
    breaks = seq(0, 8000, 1000),
    labels = function(x) format(x, big.mark = ".", decimal.mark = ",",
                                scientific = FALSE),
    expand = c(0, 0),
    # Eixo Y direito: anuais (inverso da transformação)
    sec.axis = sec_axis(
      trans  = ~ . / escala,
      name   = "Publicações anuais",
      breaks = seq(0, 800, 100)
    )
  ) +

  scale_x_discrete(guide = guide_axis(angle = 45)) +

  labs(
    x       = "Ano de publicação",
    caption = paste0(
      "Fonte: NHGRI-EBI GWAS Catalog REST API v2 (consulta: maio/2026). ",
      "Barras verdes: publicações acumuladas (eixo esquerdo).\n",
      "Linha vermelha: publicações anuais (eixo direito). ",
      "Referência original: Welter et al. (2014)."
    )
  ) +

  theme_classic(base_size = 12) +
  theme(
    axis.title.x       = element_text(margin = margin(t = 8)),
    axis.title.y.left  = element_text(colour = COR_BARRA, margin = margin(r = 8)),
    axis.title.y.right = element_text(colour = COR_LINHA, margin = margin(l = 8)),
    axis.text.y.left   = element_text(colour = COR_BARRA, size = 9),
    axis.text.y.right  = element_text(colour = COR_LINHA, size = 9),
    axis.line.y.left   = element_line(colour = COR_BARRA, linewidth = 0.4),
    axis.line.y.right  = element_line(colour = COR_LINHA, linewidth = 0.4),
    axis.ticks.y.left  = element_line(colour = COR_BARRA, linewidth = 0.4),
    axis.ticks.y.right = element_line(colour = COR_LINHA, linewidth = 0.4),
    axis.line.x        = element_line(colour = "black",   linewidth = 0.4),
    axis.ticks.x       = element_line(colour = "black",   linewidth = 0.4),
    axis.text.x        = element_text(size = 9),
    plot.caption       = element_text(size = 7, colour = "grey40",
                                      hjust = 0, margin = margin(t = 8)),
    plot.margin        = margin(t = 10, r = 20, b = 6, l = 6)
  )

print(p)

# ------------------------------------------------------------------------------
# 3. Exportar PDF
# ------------------------------------------------------------------------------
OUTPUT_FILE <- "d:/Dropbox/SMS - Classificação/Códigos em R/Códigos usados no TCC/TCC - Izabela/figura1_gwas_real.pdf"

ggsave(filename = OUTPUT_FILE, plot = p, device = "pdf",
       width = 17, height = 10, units = "cm")

message("Figura exportada: ", OUTPUT_FILE)
