# ══════════════════════════════════════════════════════════════════════════════
# PIPELINE COMPLETO — Estatística Espacial para CP e AP
# Ordem de execução: bloco a bloco, de cima para baixo
# ══════════════════════════════════════════════════════════════════════════════

# ── Pacotes ───────────────────────────────────────────────────────────────────
library(geobr)
library(sf)
library(spdep)
library(tidyverse)
library(sidrar)
library(ipeadatar)
library(basedosdados)
library(spatialreg)
library(GWmodel)

# ══════════════════════════════════════════════════════════════════════════════
# BLOCO 1 — Malha municipal
# ══════════════════════════════════════════════════════════════════════════════

sf::sf_use_s2(FALSE)

municipios <- geobr::read_municipality(year = 2010, simplified = TRUE) |>
  dplyr::mutate(code_muni = as.character(code_muni)) |>
  sf::st_make_valid()

# ══════════════════════════════════════════════════════════════════════════════
# BLOCO 2 — Coleta e integração dos dados
# ══════════════════════════════════════════════════════════════════════════════

# ── Analfabetismo ─────────────────────────────────────────────────────────────
analfabetismo_raw <- sidrar::get_sidra(
  x          = 1383,
  variable   = 1646,
  geo        = "City",
  geo.filter = NULL,
  period     = "2010"
)

analfabetismo <- analfabetismo_raw |>
  filter(`Sexo (Código)` == 6794) |>
  transmute(
    code_muni        = as.character(`Município (Código)`),
    tx_analfabetismo = 100 - as.numeric(Valor)
    # variável 1646 = taxa de ALFABETIZAÇÃO; invertendo para analfabetismo
  )

# ── IDH-M renda ───────────────────────────────────────────────────────────────
atlas <- ipeadatar::ipeadata(code = "ADH_IDHM_R") |>
  filter(
    date == as.Date("2010-01-01"),
    is.na(uname)   # uname NA = municípios; Brazil = 1 obs; States = 27 obs
  ) |>
  transmute(
    code_muni  = as.character(tcode),
    idhm_renda = value
  )

# ── Cobertura ESF ─────────────────────────────────────────────────────────────
basedosdados::set_billing_id("nomadic-asset-440300-k2")  # substituir pelo seu ID

esf <- basedosdados::read_sql(
  "SELECT
     id_municipio AS code_muni,
     AVG(proporcao_cobertura_estrategia_saude_familia) AS cobertura_esf
   FROM `basedosdados.br_ms_atencao_basica.municipio`
   WHERE ano = 2010
   GROUP BY id_municipio"
) |>
  dplyr::mutate(
    code_muni     = as.character(code_muni),
    cobertura_esf = as.numeric(cobertura_esf)
  )

# ── População total ───────────────────────────────────────────────────────────
pop <- ipeadatar::ipeadata(code = "POPTOT") |>
  filter(
    date == as.Date("2010-01-01"),
    nchar(as.character(tcode)) == 7  # 7 dígitos = municípios
  ) |>
  transmute(
    code_muni = as.character(tcode),
    pop_total = as.numeric(value)
  )

# ── Join completo ─────────────────────────────────────────────────────────────
dados_municipios <- municipios |>
  left_join(analfabetismo, by = "code_muni") |>
  left_join(atlas,         by = "code_muni") |>
  left_join(esf,           by = "code_muni") |>
  left_join(pop,           by = "code_muni") |>
  filter(!is.na(idhm_renda)) |>  # remove Nazária (PI), criada após 2010
  sf::st_make_valid()

# Diagnóstico
dados_municipios |>
  sf::st_drop_geometry() |>
  summarise(
    n_total          = n(),
    na_analfabetismo = sum(is.na(tx_analfabetismo)),
    na_idhm          = sum(is.na(idhm_renda)),
    na_esf           = sum(is.na(cobertura_esf)),
    na_pop           = sum(is.na(pop_total))
  ) |>
  print()
# esperado: n_total = 5564, todos os NAs = 0

# ══════════════════════════════════════════════════════════════════════════════
# BLOCO 3 — Matriz de pesos W (construída APÓS o filtro de dados)
# ══════════════════════════════════════════════════════════════════════════════

# W deve ter o mesmo N de dados_municipios — sempre reconstruir após filtros
vizinhos <- spdep::poly2nb(dados_municipios, queen = TRUE)
W        <- spdep::nb2listw(vizinhos, style = "W", zero.policy = TRUE)

summary(vizinhos)

# ══════════════════════════════════════════════════════════════════════════════
# BLOCO 4 — Índice de Moran global e diagrama de dispersão
# ══════════════════════════════════════════════════════════════════════════════

set.seed(42)

moran_analf <- spdep::moran.mc(
  x           = dados_municipios$tx_analfabetismo,
  listw       = W,
  nsim        = 999,
  zero.policy = TRUE
)

print(moran_analf)
# resultado real: I = 0.876, p = 0.001
# atualizar Tabela 1 e texto da Seção 4.2 no artigo

spdep::moran.plot(
  x      = c(scale(dados_municipios$tx_analfabetismo)),
  listw  = W,
  labels = FALSE,
  pch    = 20,
  col    = "steelblue",
  xlab   = "Analfabetismo (padronizado)",
  ylab   = "Lag espacial (padronizado)",
  main   = "Diagrama de Moran – Analfabetismo Municipal (2010)"
)

# ══════════════════════════════════════════════════════════════════════════════
# BLOCO 5 — LISA e mapeamento
# ══════════════════════════════════════════════════════════════════════════════

lisa <- spdep::localmoran_perm(
  x           = dados_municipios$tx_analfabetismo,
  listw       = W,
  nsim        = 999,
  zero.policy = TRUE
)

col_pval <- grep("Pr\\(z", colnames(lisa), value = TRUE)[1]

dados_municipios$lisa_Ii   <- lisa[, "Ii"]
dados_municipios$lisa_pval <- lisa[, col_pval]

z_analf  <- c(scale(dados_municipios$tx_analfabetismo))
wz_analf <- spdep::lag.listw(W, z_analf, zero.policy = TRUE)

dados_municipios$quadrante <- dplyr::case_when(
  z_analf >  0 & wz_analf >  0 & dados_municipios$lisa_pval < 0.05 ~ "Alto-Alto",
  z_analf <  0 & wz_analf <  0 & dados_municipios$lisa_pval < 0.05 ~ "Baixo-Baixo",
  z_analf >  0 & wz_analf <  0 & dados_municipios$lisa_pval < 0.05 ~ "Alto-Baixo",
  z_analf <  0 & wz_analf >  0 & dados_municipios$lisa_pval < 0.05 ~ "Baixo-Alto",
  TRUE ~ "Não significativo"
)

table(dados_municipios$quadrante)  # diagnóstico antes de mapear

cores <- c(
  "Alto-Alto"         = "#d73027",
  "Baixo-Baixo"       = "#4575b4",
  "Alto-Baixo"        = "#fdae61",
  "Baixo-Alto"        = "#abd9e9",
  "Não significativo" = "#f0f0f0"
)

estados <- geobr::read_state(year = 2010, simplified = TRUE)

dados_municipios |>
  dplyr::select(code_muni, quadrante) |>
  sf::st_drop_geometry() |>
  dplyr::right_join(municipios |> filter(code_muni != "2206720"), by = "code_muni") |>
  sf::st_as_sf() |>
  ggplot() +
  geom_sf(aes(fill = quadrante), color = NA) +
  geom_sf(data = estados, fill = NA, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = cores, name = "Cluster LISA") +
  theme_void()

# ══════════════════════════════════════════════════════════════════════════════
# BLOCO 6 — OLS, testes LM, SLM, SEM e GWR
# ══════════════════════════════════════════════════════════════════════════════

# ── OLS e testes LM ───────────────────────────────────────────────────────────
modelo_ols <- lm(
  tx_analfabetismo ~ idhm_renda + cobertura_esf + log(pop_total),
  data = dados_municipios
)

spdep::lm.morantest(modelo_ols, W, zero.policy = TRUE)

testes_lm <- spdep::lm.LMtests(
  model = modelo_ols,
  listw = W,
  test  = c("LMlag", "LMerr", "RLMlag", "RLMerr")
)
summary(testes_lm)

# ── SLM ───────────────────────────────────────────────────────────────────────
modelo_slm <- spatialreg::lagsarlm(
  tx_analfabetismo ~ idhm_renda + cobertura_esf + log(pop_total),
  data        = dados_municipios,
  listw       = W,
  zero.policy = TRUE
)

spatialreg::impacts(modelo_slm, listw = W, R = 500) |> summary(zstats = TRUE)

# ── SEM ───────────────────────────────────────────────────────────────────────
modelo_sem <- spatialreg::errorsarlm(
  tx_analfabetismo ~ idhm_renda + cobertura_esf + log(pop_total),
  data        = dados_municipios,
  listw       = W,
  zero.policy = TRUE
)

# nobs não tem método para spatialreg — usar length() no lugar
cat("N OLS:", nobs(modelo_ols), "\n")
cat("N SLM:", length(modelo_slm$residuals), "\n")
cat("N SEM:", length(modelo_sem$residuals), "\n")

# Resumo dos modelos
summary(modelo_slm)
summary(modelo_sem)

# Moran nos resíduos — verificar se a autocorrelação foi absorvida
spdep::moran.mc(
  x           = residuals(modelo_ols),
  listw       = W,
  nsim        = 999,
  zero.policy = TRUE
) |> print()

spdep::moran.mc(
  x           = residuals(modelo_slm),
  listw       = W,
  nsim        = 999,
  zero.policy = TRUE
) |> print()

spdep::moran.mc(
  x           = residuals(modelo_sem),
  listw       = W,
  nsim        = 999,
  zero.policy = TRUE
) |> print()

AIC(modelo_ols, modelo_slm, modelo_sem)


# Extrair coeficientes para montar a tabela
library(broom)

ols_tidy <- broom::tidy(modelo_ols) |>
  mutate(modelo = "OLS")

slm_tidy <- broom::tidy(modelo_slm) |>
  mutate(modelo = "SLM")

sem_tidy <- broom::tidy(modelo_sem) |>
  mutate(modelo = "SEM")

# Tabela 2 corrigida — pivot por termo e modelo sem duplicação
tabela2 <- bind_rows(ols_tidy, slm_tidy, sem_tidy) |>
  filter(term %in% c("(Intercept)", "idhm_renda", 
                     "cobertura_esf", "log(pop_total)")) |>
  mutate(
    coef = paste0(
      round(estimate, 4),
      case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ ""
      ),
      "\n(", round(std.error, 4), ")"
    )
  ) |>
  select(modelo, term, coef) |>
  tidyr::pivot_wider(names_from = modelo, values_from = coef)

print(tabela2)

# Linha de parâmetros espaciais e ajuste — montar manualmente
ajuste <- tibble::tibble(
  term  = c("ρ (rho – lag)", "λ (lambda – erro)",
            "I Moran resíduos", "AIC", "Log-verossimilhança",
            "N"),
  OLS   = c("—", "—", "0,569***", "31.893", "−15.942", "5.564"),
  SLM   = c("0,6824***", "—", "0,054***", "27.199", "−13.594", "5.564"),
  SEM   = c("—", "0,8994***", "−0,119 n.s.", "27.112", "−13.550", "5.564")
)

tabela_final <- bind_rows(tabela2, ajuste)
print(tabela_final, n = Inf)


library(knitr)

tabela_final |>
  dplyr::rename(
    "Variável"             = term,
    "OLS"                  = OLS,
    "SLM (Lag)"            = SLM,
    "SEM (Erro)"           = SEM
  ) |>
  knitr::kable(
    format  = "pipe",
    caption = "Tabela 2 – Resultados dos modelos de regressão (VD: taxa de analfabetismo, N = 5.564)"
  ) |>
  print()


tabela_word <- tabela_final |>
  mutate(across(everything(), ~stringr::str_replace_all(., "\n", " "))) |>
  dplyr::rename(
    "Variável"   = term,
    "SLM (Lag)"  = SLM,
    "SEM (Erro)" = SEM
  )

knitr::kable(
  tabela_word,
  format  = "pipe",
  caption = "Tabela 2 – Resultados dos modelos de regressão (VD: taxa de analfabetismo, N = 5.564)"
)

# ── GWR ───────────────────────────────────────────────────────────────────────
install.packages("GWmodel")
library(GWmodel)

dados_sp <- dados_municipios |>
  sf::st_transform(crs = 5880) |>
  as("Spatial")

bw <- GWmodel::bw.gwr(
  tx_analfabetismo ~ idhm_renda + cobertura_esf + log(pop_total),
  data     = dados_sp,
  approach = "CV",
  kernel   = "gaussian",
  adaptive = FALSE
)

modelo_gwr <- GWmodel::gwr.basic(
  tx_analfabetismo ~ idhm_renda + cobertura_esf + log(pop_total),
  data     = dados_sp,
  bw       = bw,
  kernel   = "gaussian",
  adaptive = FALSE
)

coefs_gwr <- as.data.frame(modelo_gwr$SDF)
print(head(coefs_gwr))




# Resumo dos coeficientes locais da GWR
summary(coefs_gwr[, c("idhm_renda", "cobertura_esf", "log(pop_total)", "Local_R2")])

# Intervalos por variável para o texto do artigo
cat("\n── idhm_renda (coeficientes locais) ─────────\n")
cat("Mínimo:", round(min(coefs_gwr$idhm_renda), 1), "\n")
cat("Mediana:", round(median(coefs_gwr$idhm_renda), 1), "\n")
cat("Máximo:", round(max(coefs_gwr$idhm_renda), 1), "\n")
cat("Q1:", round(quantile(coefs_gwr$idhm_renda, 0.25), 1), "\n")
cat("Q3:", round(quantile(coefs_gwr$idhm_renda, 0.75), 1), "\n")

cat("\n── cobertura_esf (coeficientes locais) ──────\n")
cat("Mínimo:", round(min(coefs_gwr$cobertura_esf), 4), "\n")
cat("Mediana:", round(median(coefs_gwr$cobertura_esf), 4), "\n")
cat("Máximo:", round(max(coefs_gwr$cobertura_esf), 4), "\n")

cat("\n── Local_R2 ──────────────────────────────────\n")
cat("Mínimo:", round(min(coefs_gwr$Local_R2), 3), "\n")
cat("Mediana:", round(median(coefs_gwr$Local_R2), 3), "\n")
cat("Máximo:", round(max(coefs_gwr$Local_R2), 3), "\n")

# Mapas dos coeficientes locais para o Apêndice A
library(ggplot2)

# Juntar coeficientes com geometria
dados_municipios$gwr_idhm  <- coefs_gwr$idhm_renda
dados_municipios$gwr_esf   <- coefs_gwr$cobertura_esf
dados_municipios$gwr_r2    <- coefs_gwr$Local_R2


# Figura A1 — Distribuição municipal da taxa de analfabetismo (2010)
ggplot(dados_municipios) +
  geom_sf(aes(fill = tx_analfabetismo), color = NA) +
  geom_sf(data = estados, fill = NA, color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low      = "#ffffb2",
    high     = "#bd0026",
    name     = "Analfabetismo (%)",
    breaks   = c(0, 10, 20, 30, 40),
    labels   = c("0", "10", "20", "30", "40+")
  ) +
  labs(title = "Figura A1 – Taxa de analfabetismo municipal (2010)") +
  theme_void() +
  theme(legend.position = "bottom")



# Figura A2 — coeficientes locais IDH-M renda
ggplot(dados_municipios) +
  geom_sf(aes(fill = gwr_idhm), color = NA) +
  geom_sf(data = estados, fill = NA, color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = "#d73027",
    mid      = "white",
    high     = "#4575b4",
    midpoint = median(coefs_gwr$idhm_renda),
    name     = "β IDH-M renda"
  ) +
  labs(title = "Figura A2 – Coeficientes locais GWR: efeito do IDH-M renda") +
  theme_void()

# Figura A3 — coeficientes locais cobertura ESF
ggplot(dados_municipios) +
  geom_sf(aes(fill = gwr_esf), color = NA) +
  geom_sf(data = estados, fill = NA, color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = "#4575b4",
    mid      = "white",
    high     = "#d73027",
    midpoint = 0,
    name     = "β Cobertura ESF"
  ) +
  labs(title = "Figura A3 – Coeficientes locais GWR: efeito da cobertura ESF") +
  theme_void()


# Salvar figuras em alta resolução (300 dpi)
ggsave("figura_A1_analfabetismo.png",
       plot   = last_plot(),
       width  = 18, height = 20,
       units  = "cm", dpi = 300)

# Regerar e salvar A2 e A3
fig_A2 <- ggplot(dados_municipios) +
  geom_sf(aes(fill = gwr_idhm), color = NA) +
  geom_sf(data = estados, fill = NA, color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#d73027", mid = "white", high = "#4575b4",
    midpoint = median(coefs_gwr$idhm_renda),
    name = "β IDH-M renda"
  ) +
  labs(title = "Figura A2 – Coeficientes locais GWR: efeito do IDH-M renda") +
  theme_void()

fig_A3 <- ggplot(dados_municipios) +
  geom_sf(aes(fill = gwr_esf), color = NA) +
  geom_sf(data = estados, fill = NA, color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#4575b4", mid = "white", high = "#d73027",
    midpoint = 0, name = "β Cobertura ESF"
  ) +
  labs(title = "Figura A3 – Coeficientes locais GWR: efeito da cobertura ESF") +
  theme_void()

ggsave("figura_A2_gwr_idhm.png",  plot = fig_A2,
       width = 18, height = 20, units = "cm", dpi = 300)
ggsave("figura_A3_gwr_esf.png",   plot = fig_A3,
       width = 18, height = 20, units = "cm", dpi = 300)
