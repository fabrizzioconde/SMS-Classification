# -*- coding: utf-8 -*-
"""
Insere citações bibliográficas nos pontos relevantes do texto.
"""

with open('Relatorio_Completo_SMS_TCC.tex', encoding='utf-8') as f:
    content = f.read()

# Lista de substituições: (trecho único, novo trecho com citação)
subs = [

    # ---- Geração das bases: pacote scrime --------------------------------
    (
        r'função \texttt{simulateSNPglm} do pacote \texttt{scrime} do R.',
        r'função \texttt{simulateSNPglm} do pacote \texttt{scrime} do R \cite{Schwender2013,scrime2013}.'
    ),

    # ---- Pipeline SMS: definição do método -------------------------------
    (
        r'\novo{O método SMS (\textit{Sequential Multi-Step for SNP selection}) é executado em três etapas sequenciais para cada kernel SVM:}',
        r'\novo{O método SMS (\textit{Sequential Multi-Step for SNP selection}) \cite{de2014snps,Souza2022,Baptista2024} é executado em três etapas sequenciais para cada kernel SVM:}'
    ),

    # ---- RF: Ranking por Floresta Aleatória ------------------------------
    (
        r'\novo{\textbf{Ranking por Floresta Aleatória (RF):} treina-se uma RF com \texttt{ntree} = 4000 árvores sobre o conjunto completo de SNPs e ordena-se os marcadores pela importância média (\textit{Mean Decrease in Gini}), do mais ao menos relevante.}',
        r'\novo{\textbf{Ranking por Floresta Aleatória (RF) \cite{breiman2001random,RF}:} treina-se uma RF com \texttt{ntree} = 4000 árvores sobre o conjunto completo de SNPs e ordena-se os marcadores pela importância média (\textit{Mean Decrease in Gini}), do mais ao menos relevante.}'
    ),

    # ---- SVM: Cutting Step -----------------------------------------------
    (
        r'\novo{\textbf{Etapa de Corte do SVM (\textit{Cutting Step}):} o SVM é treinado e avaliado iterativamente: a cada iteração, o SNP de menor rank é removido e o desempenho (acurácia ou F1 Score, conforme o experimento) é registrado. O subconjunto de SNPs que maximiza essa métrica é retido para a próxima etapa.}',
        r'\novo{\textbf{Etapa de Corte do SVM (\textit{Cutting Step}) \cite{vapnik95,chang2011libsvm,e1071}:} o SVM é treinado e avaliado iterativamente: a cada iteração, o SNP de menor rank é removido e o desempenho (acurácia ou F1 Score, conforme o experimento) é registrado. O subconjunto de SNPs que maximiza essa métrica é retido para a próxima etapa.}'
    ),

    # ---- GA: Refinamento -------------------------------------------------
    (
        r'\novo{\textbf{Refinamento por Algoritmo Genético (GA):} sobre o subconjunto retido, o GA busca o subconjunto ótimo de SNPs que maximiza o fitness do SVM (acurácia ou F1 Score). A aptidão é avaliada por validação (hold-out ou $k$-fold estratificado, conforme o experimento). O resultado final de cada kernel é o subconjunto identificado pelo GA.}',
        r'\novo{\textbf{Refinamento por Algoritmo Genético (GA) \cite{Goldberg1989,scrucca2012ga}:} sobre o subconjunto retido, o GA busca o subconjunto ótimo de SNPs que maximiza o fitness do SVM (acurácia ou F1 Score). A aptidão é avaliada por validação (hold-out ou $k$-fold estratificado, conforme o experimento). O resultado final de cada kernel é o subconjunto identificado pelo GA.}'
    ),

    # ---- Kernel RBF: definição formal ------------------------------------
    (
        r'\novo{No \textit{Support Vector Machine} com kernel de Função de Base Radial (RBF), a função de kernel é definida como:}',
        r'\novo{No \textit{Support Vector Machine} (SVM) com kernel de Função de Base Radial (RBF) \cite{vapnik95,burges1998tutorial,cristianini2000introduction}, a função de kernel é definida como:}'
    ),

    # ---- gamma como regularizador implícito -----------------------------
    (
        r'\novo{Embora a regularização explícita no SVM seja controlada pelo parâmetro $C$ (custo de violação da margem), $\gamma$ atua como um \textbf{regulador implícito da complexidade do modelo} no espaço de características.',
        r'\novo{Embora a regularização explícita no SVM seja controlada pelo parâmetro $C$ (custo de violação da margem) \cite{smola2004tutorial}, $\gamma$ atua como um \textbf{regulador implícito da complexidade do modelo} no espaço de características.'
    ),

    # ---- Métricas: F1 Score é mais adequado para dados desbalanceados ---
    (
        r'\novo{Neste relatório são apresentadas duas métricas complementares para avaliar os métodos de seleção de SNPs. Cada uma actua em um nível de análise distinto, e sua interpretação conjunta é fundamental para uma avaliação adequada dos resultados.}',
        r'\novo{Neste relatório são apresentadas duas métricas complementares para avaliar os métodos de seleção de SNPs \cite{fawcett2006introduction,Adeline2015}. Cada uma actua em um nível de análise distinto, e sua interpretação conjunta é fundamental para uma avaliação adequada dos resultados.}'
    ),

    # ---- k-fold estratificado: validação cruzada ------------------------
    (
        r'A aptidão é avaliada por validação (hold-out ou $k$-fold estratificado, conforme o experimento).',
        r'A aptidão é avaliada por validação (hold-out ou $k$-fold estratificado, conforme o experimento) \cite{Kohavi1995}.'
    ),

    # ---- SNPs em estudos de associação (linha de introdução ao texto) ---
    (
        r'As tabelas evidenciam os dados obtidos a partir da primeira simulação, realizada com apenas oito efeitos aditivos, sendo SNP1, SNP2, SNP3, SNP4, SNP5, SNP6, SNP7 e SNP8 utilizados nessa simulação como SNPs informativos para a geração do fenótipo.',
        r'As tabelas evidenciam os dados obtidos a partir da primeira simulação \cite{de2014snps,Souza2022}, realizada com apenas oito efeitos aditivos, sendo SNP1, SNP2, SNP3, SNP4, SNP5, SNP6, SNP7 e SNP8 utilizados nessa simulação como SNPs informativos para a geração do fenótipo.'
    ),

    # ---- Feature selection: seleção de variáveis em bioinformática ------
    (
        r'\novo{As métricas de avaliação das variáveis (SNPs) são \textbf{matematicamente idênticas} às métricas clássicas de avaliação de classificadores, diferindo apenas no \emph{universo} e na definição de ``positivo''.',
        r'\novo{As métricas de avaliação das variáveis (SNPs) \cite{saeys2007review,guyon2003introduction} são \textbf{matematicamente idênticas} às métricas clássicas de avaliação de classificadores, diferindo apenas no \emph{universo} e na definição de ``positivo''.'
    ),

    # ---- RF para estudos genéticos --------------------------------------
    (
        r'\novo{\textbf{Ranking por Floresta Aleatória (RF) \cite{breiman2001random,RF}:} treina-se uma RF com \texttt{ntree} = 4000 árvores sobre o conjunto completo de SNPs e ordena-se os marcadores pela importância média (\textit{Mean Decrease in Gini}), do mais ao menos relevante.}',
        r'\novo{\textbf{Ranking por Floresta Aleatória (RF) \cite{breiman2001random,RF,goldstein2011random}:} treina-se uma RF com \texttt{ntree} = 4000 árvores sobre o conjunto completo de SNPs e ordena-se os marcadores pela importância média (\textit{Mean Decrease in Gini}), do mais ao menos relevante.}'
    ),

    # ---- Softare R -------------------------------------------------------
    (
        r'do pacote \texttt{scrime} do R \cite{Schwender2013,scrime2013}.',
        r'do pacote \texttt{scrime} do R \cite{Team2013,Schwender2013,scrime2013}.'
    ),
]

ok = 0
fail = 0
for old, new in subs:
    if old in content:
        content = content.replace(old, new, 1)
        ok += 1
    else:
        print(f"AVISO: não encontrado: {repr(old[:70])}")
        fail += 1

print(f"\nCitações inseridas: {ok}  |  Não encontradas: {fail}")

with open('Relatorio_Completo_SMS_TCC.tex', 'w', encoding='utf-8') as f:
    f.write(content)
print("Arquivo salvo.")
