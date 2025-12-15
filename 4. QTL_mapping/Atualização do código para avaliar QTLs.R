##########################################
# Atualização do código para avaliar QTLs
# Bruna Marques Moreno 10/10/2025
##########################################

# Configuração inicial
library(qtl)
library(dplyr)
library(knitr)

# Verificar se os arquivos existem

# Função para análise QTL
analisar_qtl <- function(trait, window_size, step_size, n_perm = 1000, alpha = 0.05, covar) {
  
  cat("Analisando trait:", trait, "| Window:", window_size, "| Step:", step_size, "\n")
  
  # Carregar e preparar dados
  mimulus.qtl <- read.cross("csv", file = "~/Mestrado/Artigo/artigo-jamapa-puebla/Linkage_map/qtl_4_1000_5.csv", 
                            genotypes = c("A", "H", "B"))
  
  mimulus.qtl <- jittermap(mimulus.qtl)
  
  # Inputação de genótipo
  set.seed(22334455)
  mimulus.qtl <- fill.geno(mimulus.qtl)
  
  # Calcular probabilidades genotípicas
  mimulus.qtl.prob <- calc.genoprob(mimulus.qtl, step = step_size, 
                                    map.function = "kosambi", stepwidth = "fixed")
  
  resultados <- list()
  
  # 1. IM (Interval Mapping)
  cat("  Executando IM...\n")
  scan_result_hk <- scanone(mimulus.qtl.prob, method = "hk", pheno.col = trait)
  operm.hk <- scanone(mimulus.qtl.prob, method = "hk", n.perm = n_perm, pheno.col = trait)
  threshold_im <- as.numeric(summary(operm.hk, alpha = alpha))
  
  QTLs_im <- summary(scan_result_hk, threshold = threshold_im)
  
  if(nrow(QTLs_im) > 0) {
    qtlsim <- makeqtl(mimulus.qtl.prob, chr = QTLs_im$chr, pos = QTLs_im$pos, what = "prob")
    model1 <- fitqtl(mimulus.qtl.prob, pheno.col = trait, qtl = qtlsim, get.ests = TRUE)
    resumo <- summary(model1)
    
    IM_QTLs <- data.frame(
      "marcador" = rownames(QTLs_im),
      "cromossomo" = QTLs_im$chr, 
      "posição" = QTLs_im$pos, 
      "LOD" = QTLs_im$lod,
      "efeito_aditivo" = if(!is.null(resumo$ests)) resumo$ests[2,1] else NA,
      "efeito_dominância" = if(!is.null(resumo$ests)) resumo$ests[3,1] else NA,
      "R2" = if(!is.null(resumo$result.full)) resumo$result.full[1,5] else NA
    )
  } else {
    IM_QTLs <- data.frame()
  }
  
  # 2. CIM (Composite Interval Mapping)
  cat("  Executando CIM...\n")
  cim_result <- cim(mimulus.qtl.prob, pheno.col = trait, n.marcovar = covar, 
                    method = "hk", window = window_size, map.function = "kosambi")
  
  threshold.cim <- cim(mimulus.qtl.prob, pheno.col = trait, n.marcovar = covar, 
                       method = "hk", window = window_size, n.perm = n_perm, 
                       map.function = "kosambi")
  
  LOD_cim <- as.numeric(summary(threshold.cim, alpha = alpha))
  QTLs_cim <- summary(cim_result, threshold = LOD_cim)
  
  if(nrow(QTLs_cim) > 0) {
    qtl_data <- makeqtl(mimulus.qtl.prob, chr = QTLs_cim$chr, pos = QTLs_cim$pos)
    model_cim <- fitqtl(mimulus.qtl.prob, pheno.col = trait, qtl = qtl_data, get.ests = TRUE)
    resumo_cim <- summary(model_cim)
    
    CIM_QTLs <- data.frame(
      "marcador" = rownames(QTLs_cim),
      "cromossomo" = QTLs_cim$chr,
      "posição" = QTLs_cim$pos,
      "LOD" = QTLs_cim$lod,
      "efeito_aditivo" = if(!is.null(resumo_cim$ests)) resumo_cim$ests[2,1] else NA,
      "efeito_dominância" = if(!is.null(resumo_cim$ests)) resumo_cim$ests[3,1] else NA,
      "R2" = if(!is.null(resumo_cim$result.full)) resumo_cim$result.full[1,5] else NA
    )
  } else {
    CIM_QTLs <- data.frame()
  }
  
  # 3. MIM (Multiple Interval Mapping)
  cat("  Executando MIM...\n")
  mimulus.sim <- sim.geno(mimulus.qtl.prob, n.draws = 200, error.prob = 0.01, map.function = "kosambi", stepwidth = "fixed")
  scanbi <- scantwo(mimulus.sim, method = "hk", pheno.col = trait, 
                    use = "complete.obs")
  
  # Para MIM, usar menos permutações para maior velocidade
  threshold_mim <- scantwo(mimulus.sim, method = "hk", pheno.col = trait, 
                           n.perm = 1000)  # Reduzido para performance
  
  QTLs_mim <- summary(scanbi, perms = threshold_mim, alpha = alpha, pvalues = TRUE)
  
  # Compilar resultados
  resultados <- list(
    trait = trait,
    window = window_size,
    step = step_size,
    IM = list(
      scan_result = scan_result_hk,
      threshold = threshold_im,
      QTLs = IM_QTLs
    ),
    CIM = list(
      scan_result = cim_result,
      threshold = LOD_cim,
      QTLs = CIM_QTLs
    ),
    MIM = list(
      scan_result = scanbi,
      threshold = threshold_mim,
      QTLs = QTLs_mim
    )
  )
  
  return(resultados)
}

# Função para salvar gráficos
salvar_graficos <- function(resultados, prefixo) {
  trait_name <- resultados$trait
  window_size <- resultados$window
  step_size <- resultados$step
  
  # Gráfico IM
  png(paste0(prefixo, "_IM_", trait_name, "_w", window_size, "_s", step_size, ".png"),
      width = 800, height = 600)
  plot(resultados$IM$scan_result, main = paste("IM -", trait_name))
  abline(h = resultados$IM$threshold, col = "darkred")
  dev.off()
  
  # Gráfico CIM
  png(paste0(prefixo, "_CIM_", trait_name, "_w", window_size, "_s", step_size, ".png"),
      width = 800, height = 600)
  plot(resultados$CIM$scan_result, main = paste("CIM -", trait_name))
  abline(h = resultados$CIM$threshold, col = "darkred")
  dev.off()
  
  # Gráfico MIM
  png(paste0(prefixo, "_MIM_", trait_name, "_w", window_size, "_s", step_size, ".png"),
      width = 800, height = 600)
  plot(resultados$MIM$scan_result, main = paste("MIM -", trait_name))
  dev.off()
}

# Loop principal para análise
analisar_multiplos_params <- function(traits, windows, steps, covars, n_perm = 1000, alpha = 0.05) {
  todos_resultados <- list()
  contador <- 1
  
  for(trait in traits) {
    for(window_size in windows) {
      for(step_size in steps) {
        for(covar in covars) {  # NOVO LOOP para covar
          cat("=== Configuração", contador, "===\n")
          cat("Trait:", trait, "| Window:", window_size, "| Step:", step_size, "| Covar:", covar, "\n")
          
          resultado <- analisar_qtl(trait, window_size, step_size, n_perm, alpha, covar)
          resultado$covar <- covar  # Adicionar covar aos resultados
          todos_resultados[[contador]] <- resultado
          
          # Salvar gráficos com info do covar
          salvar_graficos(resultado, "analise")
          
          # Salvar sumários em arquivos com info do covar
          if(nrow(resultado$IM$QTLs) > 0) {
            write.csv(resultado$IM$QTLs, 
                      paste0("IM_", trait, "_w", window_size, "_s", step_size, "_c", covar, ".csv"))
          }
          if(nrow(resultado$CIM$QTLs) > 0) {
            write.csv(resultado$CIM$QTLs, 
                      paste0("CIM_", trait, "_w", window_size, "_s", step_size, "_c", covar, ".csv"))
          }
          
          contador <- contador + 1
          cat("--- Concluído ---\n\n")
        }
      }
    }
  }
  
  return(todos_resultados)
}

# Atualizar função de resumo para incluir covar
resumir_resultados <- function(resultados) {
  resumo <- data.frame()
  
  for(res in resultados) {
    linha <- data.frame(
      trait = res$trait,
      window = res$window,
      step = res$step,
      covar = res$covar,  # NOVA COLUNA
      QTLs_IM = nrow(res$IM$QTLs),
      QTLs_CIM = nrow(res$CIM$QTLs),
      LOD_threshold_IM = res$IM$threshold,
      LOD_threshold_CIM = res$CIM$threshold
    )
    resumo <- rbind(resumo, linha)
  }
  
  return(resumo)
}

# Atualizar função salvar_graficos para incluir covar
salvar_graficos <- function(resultados, prefixo) {
  trait_name <- resultados$trait
  window_size <- resultados$window
  step_size <- resultados$step
  covar <- resultados$covar  # NOVO
  
  # Gráfico IM
  png(paste0(prefixo, "_IM_", trait_name, "_w", window_size, "_s", step_size, "_c", covar, ".png"),
      width = 800, height = 600)
  plot(resultados$IM$scan_result, main = paste("IM -", trait_name, "| Covar:", covar))
  abline(h = resultados$IM$threshold, col = "darkred")
  dev.off()
  
  # Gráfico CIM
  png(paste0(prefixo, "_CIM_", trait_name, "_w", window_size, "_s", step_size, "_c", covar, ".png"),
      width = 800, height = 600)
  plot(resultados$CIM$scan_result, main = paste("CIM -", trait_name, "| Covar:", covar))
  abline(h = resultados$CIM$threshold, col = "darkred")
  dev.off()
  
  # Gráfico MIM
  png(paste0(prefixo, "_MIM_", trait_name, "_w", window_size, "_s", step_size, "_c", covar, ".png"),
      width = 800, height = 600)
  plot(resultados$MIM$scan_result, main = paste("MIM -", trait_name, "| Covar:", covar))
  dev.off()
}

# Definir parâmetros incluindo diferentes valores de covar
traits <- c("NMO", "IG")
windows <- c(10, 25)
steps <- c(0.5, 1)
covars <- c(1,2,3,4)  # Diferentes valores de n.marcovar para testar
n_perm <- 1000

# Executar análise com todas as combinações
resultados_completos <- analisar_multiplos_params(
  traits = traits,
  windows = windows,
  steps = steps,
  covars = covars,
  n_perm = n_perm
)

# Gerar resumo comparativo
resumo_geral <- resumir_resultados(resultados_completos)
print(resumo_geral)

# Salvar resumo completo
write.csv(resumo_geral, "resumo_comparativo_covar.csv", row.names = FALSE)