
#'Projeto
#'Como a qualidade da água influencia o uso e as práticas adaptativas
#'nas em comunidades do semiárido?

#'Hipóteses
#'H1 - A qualidade da água influencia o desenvolvimento de práticas adaptativas
#'H2 - O acesso a água potável influencia o desenvolvimento de práticas adaptativas
#'H3 - A qualidade da água influencia seu uso
#'H4 - O acesso a água potável influencia o consumo humano

### ANÁLISES EXPLORATÓRIAS ###

#Definir diretório
setwd(choose.dir())
dir()

#Pacotes
library(tidyverse) #Exploração e Transformação dos dados
library(FactoMineR) #Analise exploratória multivariada
library(factoextra) #Vizualização dos resultados da PCA
library(psych) # Análise da PCA
library(corrplot) #Verificação da PCA
library(vegan) #Padronização e Transformação da Análise multidimensional
library(ggplot2) #Visualização dos dados
library(DHARMa) # Diagnose dos Modelos Lineares Generalizados
library(performance) # Análise da sobredispersão GLM
library(MASS) # Ajuste da distribuição do Modelo
library(piecewiseSEM) # Análise de distribuição não normal
library(gridExtra) # Visualização dos modelos


## PCA ##

#Padronizar WQI para inserir no modelo

#Carregar os dados
dados <- read.csv("Data_H2O.csv",header=TRUE, check.names = FALSE)

# Criar subconjuntos ativos (Linhas) e variáveis ativas (Colunas)
Use_active<- dados [1:31, 22:26]

#Visualizar dados
View(Use_active)

#Padronização dos dados
dados_pad <- decostand(Use_active, "total")
View(dados_pad)

#Gerar PCA
env.pca <- PCA(dados_pad, scale.unit = TRUE, graph = FALSE)
var_env <- get_pca_var(env.pca)

# Contribuição (%) das variáveis para cada eixo
var_env$contrib

# Loadings - correlação das variáveis com os eixos
var_env$cor

# Qualidade da representação da variável. Esse valor é obtido multiplicado var_env$coord por var_env$coord
var_env$cos2

# Escores (posição) das localidades ("site scores") em cada eixo 
ind_env <- get_pca_ind(env.pca)

env.pca$eig

# Passo 1: obter os primeiros eixos
pred.env <- ind_env$coord[, 1:3]

# Combinar os dois valores em um único data.frame
dat <- data.frame(pred.env, dados)
#dat <- data.frame(Dim.1 = pred.env$Dim.1, dados)
View(dat)


## VISUALIZAÇÃO DA PCA ##

# Criar subconjuntos ativos (Linhas) e variáveis ativas (Colunas)
Use_active<- dados [1:31, 22:26]

#Visualizar dados
View(Use_active)

#Padronização dos dados
dados_pad <- decostand(Use_active, "total")
View(dados_pad)

#Gerar PCA
Use_PCA <- PCA(dados_pad, graph = F)

#Extrair variancia dos valores
#'porcentagem de explicação de cada eixo na variação total
Use_val <- get_eigenvalue(Use_PCA)
Use_val

windows(title="PCA biplots", 14, 7)
par(mfrow = c(1, 1))

#Obter os primerios eixos da PCA
pred.env <- Use_PCA$coord[, 1:3] 

## combinar os primeiros eixos em um único data.frame
dat <- data.frame(pred.env, dados)

#Plotar variancia
fviz_eig(Use_PCA, addlabels = T, ylim = c(0,50))

#Extrair resultados das variaveis
var <- get_pca_var(Use_PCA)
ind <- get_pca_ind(Use_PCA)

#Plotar grafico PCA
fviz_pca_var(Use_PCA, col.var = "blue")

#Criar grupo para cluster
grupo <- as.factor(dados[,14])
grupo

#Plotar Grafico biplot
fviz_pca_biplot(Use_PCA)

fviz_pca_biplot(Use_PCA, habillage = grupo, title = "Uso da água")

fviz_pca_biplot(Use_PCA,
                geom.ind = "point",
                pointshape = 21,
                pointsize = 5,
                fill.ind = dados$Area)+
  theme_classic()+
  labs(title= "Gráfico PCA qualidade da água",
       fill = "Local de coleta",
       x= "PCA1 (41%)",
       y ="PCA2 (25%)")

#Verificar qualidade da representação
var$cos2

corrplot(var$cos2, is.corr = F)


### MODELOS ###


#'Escolha das variaveis para o modelo da Regressão múltipla

mod1 <- lm(Ha_Irrigate ~ Dim.1 + Dim.2 + Dim.3, data = dat)
par(mfrow = c(2, 2))
plot(mod1) # verificar pressupostos dos modelos lineares
summary(mod1) # resultados do  teste
dimdesc(env.pca)$Dim.1 #Verificar variaveis mais importantes

anova(mod1)

mod2 <- lm(Ha_Irrigate ~ Dim.1 + Dim.3, data = dat)
par(mfrow=c(2,2))
plot(mod2)

anova(mod1,mod2)
anova(mod2)

mod3 <- lm(Ha_Irrigate  ~ Dim.1, data = dat)
par(mfrow=c(2,2))
plot(mod3)

anova(mod2,mod3)
anova(mod3)

#' Hipótese I e II
#' A qualidade da água influencia o desenvolvimento de práticas adaptativas.
#' O acesso a água potável influencia o desenvolvimento de práticas adaptativas.
mod_pois1 <- glm(Ha_Irrigate ~ Dim.1 + Temp_C + Precip_mm + Pop_Size + GDP + Hdi
                 + Gini + Pop_Acc_H2O, family = poisson(link = "log"), data = dat)
par(mfrow = c(2, 2))
plot(mod_pois1)

#Diagnose
simulationOutput <- simulateResiduals(fittedModel = mod_pois1, plot = TRUE)

par(mfrow = c(1, 1))

#Verificar overdispersion
testDispersion(mod_pois1) # linha à direita da distribuição,indica overdispersion

# Verificar overdispersion com distribuição chi-quadrado
check_overdispersion(mod_pois1) #P<0.05 indica overdispersion

summary(mod_pois1)

# Dispersion parameter próximo de 1,indica que não há overdispersion.
deviance(mod_pois1) / df.residual(mod_pois1) #559263.6

### Ajuste do modelo ###

names(dat) #Selecionar variáveis

mod_nb <- glm.nb(Ha_Irrigate ~ Dim.1 + Temp_C + Precip_mm + Pop_Size + Hdi + Pop_Acc_H2O, data = dat)

# Diagnose dos resíduos
par(mfrow = c(2, 2))
plot(mod_nb)
par(mfrow = c(1, 1))
# Dispersion parameter próximo de 1,indica que não há overdispersion.
(chat <- deviance(mod_nb) / df.residual(mod_nb)) #1.560679

# Diagnose avançada
simulationOutput <- simulateResiduals(fittedModel = mod_nb, plot = TRUE)
rsquared(mod_nb)

summary(mod_nb)

R2 <- rsq(mod_nb)

# Plot do modelo predito

#' Hipótese I
#' A qualidade da água influencia o desenvolvimento de práticas adaptativas.
ggplot(dat, aes(Hdi, Ha_Irrigate)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(method = "glm.nb", formula = y~x, se = TRUE) +
  labs(x = "IDH",
       y = "Área Irrigada") +
  theme_bw() -> plot_hdi
plot_hdi

ggplot(dat, aes(x = Precip_mm, y = Ha_Irrigate)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(method = "glm.nb", formula = y~x, se = TRUE) +
  labs(x = "Precipitação",
       y = "Área Irrigada") +
  theme_bw() -> plot_precip
plot_precip

ggplot(dat, aes(x = Pop_Acc_H2O, y = Ha_Irrigate)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(method = "glm.nb", formula = y~x, se = TRUE) +
  labs(x = "População acesso a água potável (%)",
       y = "Área Irrigada") +
  theme_bw() -> plot_acess

plot_acess

#Visualização
windows(title="H1 e H2", 14, 7)
par(mfrow = c(2, 2))

grid.arrange(plot_hdi, plot_precip, plot_acess, ncol=2)

### retirar o  outlier

dat %>% dplyr::select(Study, Ha_Irrigate) %>% arrange(-Ha_Irrigate)

dat2 <- dat %>% dplyr::filter(Study != "Reddy_2019")
mod_nb2 <- glm.nb(Ha_Irrigate ~ Dim.1 + Temp_C + Precip_mm + Pop_Size + Hdi + Pop_Acc_H2O, data = dat2)

# Diagnose dos resíduos
par(mfrow = c(2, 2))
plot(mod_nb2)
par(mfrow = c(1, 1))
# Dispersion parameter próximo de 1,indica que não há overdispersion.
(chat2 <- deviance(mod_nb2) / df.residual(mod_nb2)) #1.560679

# Diagnose avançada
simulationOutput2 <- simulateResiduals(fittedModel = mod_nb2, plot = TRUE)
rsquared(mod_nb2)

summary(mod_nb2)

### MODELOS III E IV ###
View(dat)

## Gráfico
ggplot(dat, aes(Dim.1, Agriculture)) +
  geom_point(aes(shape = Area, colour = Area), size = 4, alpha = 0.4) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(y = "Uso na Agricultura", x = "Qualidade da água (Dim.1 PCA Padronizada)", shape = "Area", colour = "Area") +
  theme_classic() -> Plot_H3
Plot_H3

ggplot(dat, aes(Pop_Acc_H2O, Agriculture)) +
  geom_point(aes(shape = Area, colour = Area), size = 4, alpha = 0.4) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(y = "Uso na Agricultura", x = "População com acesso a água potável (%)", shape = "Area", colour = "Area") +
  theme_classic() -> Plot_H4
Plot_H4

grid.arrange(Plot_H3, Plot_H4, ncol=2)
