# SMS – Simulação 6 de Oliveira (2015) com melhores técnicas do TCC

## Objetivo
Reproduzir a **Simulação 6** (classificação binária desbalanceada) da tese de
Oliveira (2015) e aplicar as melhores técnicas identificadas nas simulações
deste TCC, para comparar os resultados com o SMS original.

## Modelo simulado (Seção 7.1.6 de Oliveira 2015)

| Grupo | SNPs | Tipo de efeito | β | Codificação scrime |
|-------|------|----------------|---|-------------------|
| 1 | SNP1 | Indicador I[x=1] (AA) | 2.0 | list.ia = 1 |
| 2 | SNP2 | Indicador I[x=2] (Aa) | 1.3 | list.ia = 2 |
| 3 | SNP3 | Indicador I[x=3] (aa) | 0.9 | list.ia = 3 |
| 4 | SNP4, SNP5 | Interação ord.2: I[x4≠1]·I[x5=3] | 2.0 | list.ia = c(-1,3) |
| 5 | SNP6, SNP7, SNP8 | Interação ord.3: I[x6=1]·I[x7=2]·I[x8=3] | 3.0 | list.ia = c(1,2,3) |

- n = 1000, 100 SNPs, MAF ~ U[0.1, 0.4]
- Distribuição de classes (Oliveira 2015): **862 casos × 138 controles**
- β₀ calibrado via `uniroot` para atingir essa proporção exata

## Inovações em relação ao SMS original

1. **SMOTE-NC** aplicado exclusivamente no treino de cada fold (sem data leakage)
2. **k-fold estratificado** (k=10), mantendo a proporção 862/138 em cada fold
3. **F1-Score** como métrica de corte e fitness do GA (vs. AUC-ROC em Oliveira 2015)
4. **SVM binário explícito** (`C-classification`) vs. SVR adaptado
5. **Avaliação dual**: métricas ao nível do indivíduo (F1 obs.) **e** da variável
   (precisão, sensibilidade e F1 dos SNPs selecionados vs. causais verdadeiros)

## Kernels avaliados

| # | Kernel | γ |
|---|--------|---|
| 1 | Linear | — |
| 2 | Radial (RBF) | 0.001 |
| 3 | Radial (RBF) | 0.01 |
| 4 | Radial (RBF) | 0.1 |
| 5 | Radial (RBF) | 1.0 |

## Arquivos gerados

| Arquivo | Conteúdo |
|---------|----------|
| `dados_sim6_oliveira2015.csv` | Base simulada (1000 × 101: 100 SNPs + fenotipo) |
| `Grafico_F1_SVM_<Kernel>.pdf` | Curva F1 na trilha RF (etapa de corte) |
| `Grafico_GA_F1_<Kernel>.pdf` | Evolução do fitness (F1) no GA |
| `Resultados_SMS_Oliveira2015_Sim6.csv` | Tabela resumo por kernel |
| `SNPs_uniao_GA.csv` | União dos SNPs selecionados por todos os kernels |
| `SNPs_intersecao_GA.csv` | Interseção dos SNPs selecionados |

## Como executar

Abra o arquivo `.R` no RStudio e execute (`Source` ou `Ctrl+Shift+Enter`).
O script detecta automaticamente a pasta de trabalho.

## Referência

OLIVEIRA, Fabrízzio Condé de. *Um método para seleção de atributos em dados
genômicos*. Tese (Doutorado em Modelagem Computacional) – Universidade Federal
de Juiz de Fora, 2015. Seções 7.1.6 e 8.8.
