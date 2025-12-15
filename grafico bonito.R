# VERSÃO COM SELEÇÃO POR NÚMERO DA COLUNA

library(ggplot2)
library(dplyr)

# 1. Ler dados
dados <- read.csv("C:\\Users\\bruna\\Documents\\Mestrado\\Artigo\\artigo-jamapa-puebla\\3. Linkage_map\\qtl_4_1000_5.csv", header = TRUE, stringsAsFactors = FALSE)

# Boxplot dos marcadores

nmo <- dados[3:153, 5]
ig <- dados[3:153, 1]
geno <- dados[3:153,10]

data <- data.frame(nmo = nmo, geno = geno) # Correto: data.frame
dados[154:157,] <- NA
dados[154:157,5] <- c(134.77,75.58,43.29,51.24) #bra, jam, pue, tyb - NMO
dados[154:157,1] <- c(2.65,2.48,1.40, 1.74) #bra, jam, pue, tyb - IG


## Definir a ordem desejada
ordem_desejada <- c("A", "H", "B")  

### Boxplot dos marcadores importantes

## Converter geno para factor com ordem específica
data$geno <- factor(data$geno, levels = ordem_desejada)
data <- data[!is.na(data$nmo) & data$geno != "-" & data$geno != "" & !is.na(data$geno), ]

ggplot(data, aes(x = geno, y = nmo, fill = geno)) +
  geom_boxplot(width = 0.6, alpha = 0.7) +
  scale_fill_manual(values = c("A" = "gray", "H" = "gray", "B" = "gray")) +
  scale_x_discrete(labels = c("A" = "TT",  
                              "H" = "TG",   
                              "B" = "GG")) +
  labs(x = "Chr01_2340316",
       y = "Egg mass") +
  theme_minimal() +
  theme(legend.position = "none")

data$geno <- factor(data$geno, levels = ordem_desejada)
data <- data[!is.na(data$nmo) & data$geno != "-" & data$geno != "" & !is.na(data$geno), ]

### Pv11

geno <- dados[3:153,635]
ggplot(data, aes(x = geno, y = nmo, fill = geno)) +
  geom_boxplot(width = 0.6, alpha = 0.7) +
  scale_fill_manual(values = c("A" = "gray", "H" = "gray", "B" = "gray")) +
  scale_x_discrete(labels = c("A" = "GG",  
                              "H" = "AG",   
                              "B" = "AA")) +
  labs(x = "Chr11_3920499",
       y = "Egg mass") +
  theme_minimal() +
  theme(legend.position = "none")


#### Observação dos dados fenotípicos dos 151 indivíduos

mean(nmo)  #61.97554
# Cálculo da correlação de Pearson
correlacao_pearson <- cor(ig, nmo, method = "pearson", use = "complete.obs")
print(correlacao_pearson)
cor.test(nmo, ig)

#	Pearson's product-moment correlation

#data:  nmo and ig
#t = 6.5513, df = 149, p-value = 8.718e-10
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.3387576 0.5881968
#sample estimates:
#  cor 
#0.472898 

#correlação de 0,4670 com p-valor de 1.697e-09 e alfa de 5% de confiança para o intervalo 0.3315813 0.5836061
# Para PNG (menor qualidade)
ggsave("grafico.png", dpi = 300, width = 8.5, height = 6)

# Para EPS (algumas revistas)
ggsave("grafico.eps", device = cairo_ps, width = 8.5, height = 6)


