#### landscape and color - Camila Cason & Filipe C. Serrano ####

library(ggplot2)
library(dplyr)
library(tidyr)
library(terra)
library(sf)
library(mapview)
library(raster)
library(geobr)


#### load data #####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## carregar tabela
tabela_landuse = read.csv("tabela_completa.csv", sep = ";") %>% 
  dplyr::mutate(
    latitude = as.numeric(gsub(",", ".",latitude)),
    longitude = as.numeric(gsub(",", ".",longitude))) %>% 
  dplyr::distinct(latitude, longitude, .keep_all = T) 

table(tabela_landuse$morfotipo)


# filtrar só floresta e transformar a tabela noutro formato
tabela_landuse_forest = tabela_landuse %>% 
  dplyr::select(ID, morfotipo, starts_with("X3.Forest")) %>% # selecionar só colunas com "Forest"
  pivot_longer(cols = starts_with("X3.Forest"),
               names_to = "Year",
               values_to = "Forest_Cover",
               names_pattern = "Forest_(\\d+)") %>% 
  dplyr::filter(morfotipo != "wine") %>% # removendo "wine", so 1 registro
  dplyr::left_join(tabela_landuse %>% 
                     dplyr::select(ID, latitude, longitude))  %>% 
  dplyr::mutate(Year = as.numeric(gsub(".*_", "", Year))) %>% 
  dplyr::mutate(latitude = case_when(ID == 99 ~ -9.58445,
                                     ID == 100 ~ -7.09264,
                                     TRUE ~ latitude),
                longitude = case_when(ID == 99 ~ -55.92094,
                                      ID == 100 ~ -48.20129,
                                      TRUE ~ longitude))# %>% dplyr::filter(morfotipo != "white-spotted")

# definindo cores
frog_colors = c(
  "orange" = "orange",
  "red" = "red3",
  "yellow" = "yellow2",
  "light sky blue" = "skyblue",
  "black" = "black",
  "yellow with black vermiculations" = "gold3",
  "black with a varied cream pattern" = "grey40",
  "black with blue to green spots" = "darkcyan",
  "white-spotted" = "grey80"
)


unique(tabela_landuse_forest$morfotipo)

table(tabela_landuse_forest$morfotipo)

# plotar a mudança de forest cover ao longo dos anos
ggplot(tabela_landuse_forest, aes(x = Year, y = Forest_Cover,color = as.factor(morfotipo), group = morfotipo)) +
  stat_smooth (geom="line", alpha=0.8, size=1.3, method = "glm") +
  # geom_smooth(se = F, method = "glm", size = 2.2, aes(alpha = .4)) +
  scale_color_manual(values = frog_colors, name = "Morphotype") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)
  ) +
  labs(x = "Year",
       y = "Forest cover (nr of pixels)")


library(janitor)
library(nlme)

# corrigir alguns pontos para as distâncias entre eles não serem 0
tabela_landuse_forest_jittered = tabela_landuse_forest %>%
  mutate(
    latitude = latitude + runif(n(), -1e-7, 1e-7),
    longitude = longitude + runif(n(), -1e-7, 1e-7)
  )



# corremos um modelo linear misto 
# (misto porque estamos a incluir o ID - a avisar que nem todos os pontos são iguais porque os pontos de cada ID devem ser mais semelhantes entre si)
# incluindo essa "correlation" para avisar que alguns pontos são próximos logo também devem ser mais semelhantes do que pontos distantes


lme_spatial0 = lme(Forest_Cover ~ 1, 
                   random = ~ 1 | ID, 
                   correlation = corExp(form = ~ longitude + latitude, nugget = T),
                   data = tabela_landuse_forest_jittered)

lme_spatial1 = lme(Forest_Cover ~ Year, 
                   random = ~ 1 | ID, 
                   correlation = corExp(form = ~ longitude + latitude, nugget = T),
                   data = tabela_landuse_forest_jittered)

lme_spatial2 = lme(Forest_Cover ~ -1 +morfotipo, 
                   random = ~ 1 | ID, 
                   correlation = corExp(form = ~ longitude + latitude, nugget = T),
                   data = tabela_landuse_forest_jittered)

lme_spatial3 = lme(Forest_Cover ~ -1 + Year + morfotipo, 
                   random = ~ 1 | ID, 
                   correlation = corExp(form = ~ longitude + latitude, nugget = T),
                   data = tabela_landuse_forest_jittered)

lme_spatial_full = lme(Forest_Cover ~ -1 + Year * morfotipo, 
                       random = ~ 1 | ID, 
                       correlation = corExp(form = ~ longitude + latitude, nugget = T),
                       data = tabela_landuse_forest_jittered)



# aqui é o resultado do modelo

AIC(lme_spatial0, lme_spatial1, lme_spatial2, lme_spatial3,lme_spatial_full )

summary(lme_spatial_full)

library(emmeans)
frog_slopes = emtrends(lme_spatial_full, 
                       specs = ~ morfotipo, 
                       var = "Year")

# nesta parte comparamos as cores, usando os resultados do modelo

df_slopes_results = as.data.frame(frog_slopes)

ggplot(df_slopes_results, aes(x = reorder(morfotipo, Year.trend), y = Year.trend, 
                              color = morfotipo)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, size = 1.5) +
  geom_point(size = 3) +
  scale_color_manual(values = frog_colors) +
  coord_flip() +
  theme_classic() +
  labs(title = "",
       subtitle = "",
       x = "Morphotype",
       y = "Annual Change in Forest Pixels (Slope)") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)
  ) 


tabela_landuse_forest_sf = tabela_landuse_forest %>% 
  distinct(ID, .keep_all = T) %>% 
  sf::st_as_sf( coords = c("longitude", "latitude"), crs = 4326)

mapView(tabela_landuse_forest_sf, z = "morfotipo")


br_states = st_read("estados_BR/BR_UF_2024.shp") %>% 
  dplyr::filter(NM_REGIA == "Norte" | SIGLA_UF %in% c("MA", "MT"))
  

ggplot() +
  geom_sf(data = br_states, size = 1, fill = "NA") +
  geom_sf(data = tabela_landuse_forest_sf, color = "black", shape = 20, size = 3.4) +
  geom_sf(data = tabela_landuse_forest_sf, aes(color = morfotipo), shape = 20, size = 3) +
  scale_color_manual(values = frog_colors, name = "Morphotype") +
  theme_classic()



