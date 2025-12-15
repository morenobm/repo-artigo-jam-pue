library(qtl)

dados <- read.cross(format = "csv", 
                    file ="C:\\Users\\bruna\\Documents\\Mestrado\\Artigo\\artigo-jamapa-puebla\\3. Linkage_map\\qtl_4_1000_5.csv", 
                    sep=",", 
                    dec=".", 
                    crosstype = "f2", 
                    genotypes = c("A","H", "B") )
summary(dados)

dados <- calc.genoprob(dados, error.prob = 0.05)
head(dados$pheno)

## Mapeando sem cofatores

dupmar <- findDupMarkers(dados,adjacent.only=TRUE, exact.only = TRUE)
fill_cross <- drop.markers(dados, unlist(dupmar))
summary(fill_cross)


fill_cross <- fill.geno(fill_cross,
                        method = "argmax", 
                        error.prob = 0.05, 
                        map.function = "kosambi")

head(fill_cross$pheno)

result_nulo <- mqmscan(fill_cross, 
                       pheno.col = 5, 
                       model = "dominance",
                       window.size = 5, 
                       step.size = 1) ##MO
summary(result_nulo)


## Primeiro cofator (Pv01)
marker_do_qtl <- find.marker(fill_cross, chr=1, pos=3)
cofactor_index <- find.markerindex(fill_cross, marker_do_qtl)
cofactor_list <- mqmsetcofactors(fill_cross, cofactors = cofactor_index)
resultado_1cofactor <- mqmscan(fill_cross, 
                               cofactors = cofactor_list, 
                               pheno.col = 5,
                               window.size = 10,
                               model = "dominance") #MO

summary(resultado_1cofactor)
plot(resultado_1cofactor)

perm_test <- mqmpermutation(fill_cross,
                            scanfunction = mqmscan,
                            pheno.col = 5,
                            model = "dominance",
                            n.perm = 1000,
                            cofactors = cofactor_list,
                            n.cluster = 5)
limiar <- mqmprocesspermutation(perm_test)
summary(limiar)

plot(resultado_1cofactor);abline(h = 3.19, col = "darkred")


qtl_model <- makeqtl(dados, 
                     chr = c("1"), 
                     pos = c(5),
                     what = "prob")

fitted_model <- fitqtl(dados, 
                       qtl = qtl_model, 
                       pheno.col = 5, 
                       get.ests = TRUE,
                       formula = y ~ Q1) 
summary(fitted_model)


## Segundo cofator (Pv07)

marker_chr7 <- find.marker(fill_cross, chr=7, pos=50) # Encontrar o marcador no chr7
cofactor_index_chr7 <- find.markerindex(fill_cross, marker_chr7)

# Criar uma NOVA lista de cofactors que inclui o ANTIGO e o NOVO
new_cofactor_list <- mqmsetcofactors(fill_cross, 
                                     cofactors = c(cofactor_index, cofactor_index_chr7))

# Executar o MQMScan com DOIS cofactors
resultado_com_2_cofactors <- mqmscan(fill_cross, 
                                     cofactors = new_cofactor_list,
                                     pheno.col = 5, 
                                     model = "dominance",
                                     window.size = 5, 
                                     step.size = 1)

summary(resultado_com_2_cofactors)
plot(resultado_com_2_cofactors)

perm_2cofactor <- mqmpermutation(fill_cross,
                                 scanfunction = mqmscan,
                                 cofactors = cofactor_list,
                                 pheno.col = 5,
                                 model = "dominance",
                                 n.perm = 1000, 
                                 n.cluster = 10)
limiar_2cof <- mqmprocesspermutation(perm_2cofactor)
summary(limiar_2cof)
plot(resultado_1cofactor);abline(h = 3.18, col = "darkred")

qtl_model <- makeqtl(dados, 
                     chr = c("1","7"), 
                     pos = c(5, 50),
                     what = "prob")

fitted_model <- fitqtl(dados, 
                       qtl = qtl_model, 
                       pheno.col = 5, 
                       get.ests = TRUE,
                       formula = y ~ Q1 + Q2) 
summary(fitted_model)


## Mapeando 3 cofatores (Pv11)

marker_chr11 <- find.marker(fill_cross, chr=11, pos=56) # Encontrar o marcador no chr11
cofactor_index_chr11 <- find.markerindex(fill_cross, marker_chr11)

# Criar uma NOVA lista de cofactors que inclui o ANTIGO e o NOVO
new_cofactor_list <- mqmsetcofactors(fill_cross, cofactors = c(cofactor_index, cofactor_index_chr7,
                                                               cofactor_index_chr11))

# Executar o MQMScan com TRÊS cofactors
resultado_com_3_cofactors <- mqmscan(fill_cross, 
                                     cofactors = new_cofactor_list,
                                     pheno.col = 5, 
                                     model = "dominance",
                                     window.size = 5, 
                                     step.size = 1)

summary(resultado_com_3_cofactors)

perm_3cofactor <- mqmpermutation(fill_cross,
                                 scanfunction = mqmscan,
                                 cofactors = new_cofactor_list,
                                 pheno.col = 5,
                                 model = "dominance",
                                 n.perm = 1000, 
                                 n.cluster = 5)
limiar_3cof <- mqmprocesspermutation(perm_3cofactor)
summary(limiar_3cof)
plot(resultado_com_3_cofactors)
plot(resultado_com_3_cofactors);abline(h = 3.25 , col = "darkred")


qtl_model <- makeqtl(dados, 
                     chr = c("1","7","11"), 
                     pos = c(3,52,56),
                     what = "prob")

fitted_model <- fitqtl(dados, 
                       qtl = qtl_model, 
                       pheno.col = 5, 
                       get.ests = TRUE,
                       formula = y ~ Q1 + Q2 +Q3) 
summary(fitted_model)

## Mapeando 4 cofatores (Pv02)

marker_chr2 <- find.marker(fill_cross, chr=2, pos=16) # Encontrar o marcador no chr2
# Encontrar seu índice
cofactor_index_chr2 <- find.markerindex(fill_cross, marker_chr2)

# Criar uma NOVA lista de cofactors que inclui o ANTIGO e o NOVO
new_cofactor_list <- mqmsetcofactors(fill_cross, cofactors = c(cofactor_index,
                                                               cofactor_index_chr11,
                                                               cofactor_index_chr7,
                                                               cofactor_index_chr2))

# Executar o MQMScan com TRÊS cofactors
resultado_com_4_cofactors <- mqmscan(fill_cross, 
                                     cofactors = new_cofactor_list,
                                     pheno.col = 5, 
                                     model = "dominance",
                                     window.size = 5, 
                                     step.size = 1)

summary(resultado_com_4_cofactors)

perm_4cofactor <- mqmpermutation(fill_cross,
                                 scanfunction = mqmscan,
                                 cofactors = new_cofactor_list,
                                 pheno.col = 5,
                                 model = "dominance",
                                 n.perm = 1000, 
                                 n.cluster = 10)
limiar_4cof <- mqmprocesspermutation(perm_4cofactor)
summary(limiar_4cof)
plot(resultado_com_4_cofactors)
plot(resultado_com_4_cofactors);abline(h = 3.25 , col = "darkred")


qtl_model <- makeqtl(dados, 
                     chr = c("1","7","2","11"), 
                     pos = c(4,50,15,56),
                     what = "prob")

fitted_model <- fitqtl(dados, 
                       qtl = qtl_model, 
                       pheno.col = 5, 
                       get.ests = TRUE,
                       formula = y ~ Q1 + Q2 +Q3 +Q4) 
summary(fitted_model)

### QTLs válidos Pv01 e Pv07
marker_chr1 <- find.flanking(fill_cross, 1, 3)
marker_chr2 <- find.flanking(fill_cross, chr=2, pos=15) # Encontrar o marcador no chr2
marker_chr7 <- find.flanking(fill_cross, chr=7, pos=53) # Encontrar o marcador no chr7
marker_chr11 <- find.flanking(fill_cross, chr=11, pos=55.2) # Encontrar o marcador no chr11

effectplot(fill_cross, pheno.col=5, mname1 = "Chr11_3920499")

