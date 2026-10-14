#=======================================================================
# META-ANALYSIS OF SPIRITS DEMAND ELASTICITIES — R SCRIPT (revised version)
# Faye & Cardebat, "Spirits Demand Elasticities over Time: A Methodological Illusion for Price, a Structural Decline for Income"
#
# Input: Supplementary_data_revised.xlsx (sheets "Epdata", "Eidata", "Studies Listing")
# Output: tables and figures written to ./output
#
# Step 0  Data and corpus description
# Step 1  Full meta-analysis: model selection, publication bias
#         (Egger, PET-PEESE), CR2 cluster-robust moderators (Table 2)
# Step 2  Temporal homogeneity (Table 1, Tables S3-S4)
#         A) Pettitt (endogenous change point) + Von Neumann ratio
#         B) Distributional comparisons at 1975/1980/1985/1990
#         C) Bai-Perron (endogenous breaks) + Chow, aggregated series
# Step 3  Net period effect after compositional controls (Table 3)
#         A) Cutoff 1985     B) Cutoff grid 1970-2000
#         C) Continuous trend     D) Alternative control sets
# Step 4  Compositional change of the corpus pre/post 1985 (Table S6)
# Step 5  Robustness: estimates significant at 5% only (Table S7)
#
# Analytical decisions
# - All estimates reporting a SE, t-statistic or p-value are retained,
#   whatever their significance (Havranek et al. 2020)
# - Full models: EP random effects (REML); EI fixed effects (the ML
#   estimate of tau^2 is 0 in the full EI model, so FE = RE)
# - Reference models of Step 3: the moderator groups with at least one level
#   significant at 10% in either full model (Table 2), applied to both series
#   (country, data type, estimation method); weighting class (FE/RE)
#   re-selected by ML-consistent AIC on the reference model itself
# - CR2 cluster-robust inference throughout, clustered by study,
#   Satterthwaite small-sample degrees of freedom
# - Single-estimate methods ARC (EP) and GLM (EI) merged into OTHEREST
#=======================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(metafor)
  library(clubSandwich)
  library(strucchange)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(2026)
if (!exists("path_xlsx")) path_xlsx <- "Supplementary_data_revised.xlsx"
out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE)

#=======================================================================
# STEP 0 — DATA
#=======================================================================

prep <- function(d, merge_to_other) {
  d$Method_grp <- as.character(d$Method)
  d$Method_grp[d$Method_grp == merge_to_other] <- "OTHEREST"
  d$Method_grp <- relevel(factor(d$Method_grp), ref = "OLS")
  d$Country    <- relevel(factor(d$Country),   ref = "US")
  d$Datatype   <- relevel(factor(d$Datatype),  ref = "Pooled")
  d$Functform  <- relevel(factor(d$Functform), ref = "LIN")
  d$IDstudy    <- factor(d$IDstudy)
  d$Midpoint_c  <- d$Midpoint - mean(d$Midpoint)
  d$Midpoint_c2 <- d$Midpoint_c^2
  d$post1985    <- as.integer(d$Midpoint >= 1985)
  d[order(d$Midpoint), ]
}
Epdata <- prep(as.data.frame(read_excel(path_xlsx, sheet = "Epdata")), "ARC")
Eidata <- prep(as.data.frame(read_excel(path_xlsx, sheet = "Eidata")), "GLM")
Epdata$vi <- Epdata$SE_EP^2
Eidata$vi <- Eidata$SE_EI^2

cat("\n== Corpus ==\n")
cat(sprintf("EP: %d estimates, %d studies | EI: %d estimates, %d studies\n",
            nrow(Epdata), nlevels(Epdata$IDstudy), nrow(Eidata), nlevels(Eidata$IDstudy)))
cat(sprintf("Estimates with |t| < 1.96: EP %d | EI %d\n",
            sum(abs(Epdata$EP / Epdata$SE_EP) < 1.96), sum(abs(Eidata$EI / Eidata$SE_EI) < 1.96)))

studies <- as.data.frame(read_excel(path_xlsx, sheet = "Studies Listing"))
used    <- studies[studies$`In final dataset` == "Yes", ]
cat(sprintf("Publications screened: %d | contributing estimates: %d (%s-%s)\n",
            nrow(studies), nrow(used), min(used$Year, na.rm = TRUE), max(used$Year, na.rm = TRUE)))
for (y in c(2006, 2009, 2012))
  cat(sprintf("  contributing publications dated after %d: %d\n", y, sum(used$Year > y, na.rm = TRUE)))

formulaEP_full <- ~ Income + Advertizing + Socialdata + Habit +
  Crosselasticities + Longrun + Compensated + Conditional +
  Country + Midpoint_c + Duration + Datatype + Functform + Method_grp + Spirits
formulaEI_full <- ~ Price + Advertizing + Socialdata + Habit +
  Crosselasticities + Longrun + Compensated + Conditional +
  Country + Midpoint_c + Midpoint_c2 + Duration + Datatype + Functform +
  Method_grp + Spirits

#=======================================================================
# STEP 1 — FULL META-ANALYSIS
#=======================================================================

# 1.1 Model selection
fitEP <- list(
  FE    = rma(yi = EP, vi = vi, mods = formulaEP_full, method = "FE",   data = Epdata),
  REML  = rma(yi = EP, vi = vi, mods = formulaEP_full, method = "REML", data = Epdata),
  Multi = rma.mv(yi = EP, V = vi, mods = formulaEP_full, random = ~ 1 | IDstudy, data = Epdata))
fitEI <- list(
  FE    = rma(yi = EI, vi = vi, mods = formulaEI_full, method = "FE",   data = Eidata),
  REML  = rma(yi = EI, vi = vi, mods = formulaEI_full, method = "REML", data = Eidata),
  Multi = rma.mv(yi = EI, V = vi, mods = formulaEI_full, random = ~ 1 | IDstudy, data = Eidata,
                 control = list(optimizer = "optim", optmethod = "BFGS")))
cat("\n== 1.1 Model selection ==\n")
msel <- data.frame(Model = c("Fixed effects", "Random effects (REML)", "Multilevel (study)"),
                   AIC_EP = sapply(fitEP, AIC), logLik_EP = sapply(fitEP, logLik),
                   AIC_EI = sapply(fitEI, AIC), logLik_EI = sapply(fitEI, logLik))
print(msel, row.names = FALSE, digits = 5)
write.csv(msel, file.path(out_dir, "S2_model_selection.csv"), row.names = FALSE)
# FE vs RE on a common ML likelihood (REML and FE likelihoods are not comparable)
for (s in list(list("EP", Epdata, formulaEP_full), list("EI", Eidata, formulaEI_full))) {
  fe <- rma(yi = s[[2]][[s[[1]]]], vi = s[[2]]$vi, mods = s[[3]], method = "FE", data = s[[2]])
  ml <- rma(yi = s[[2]][[s[[1]]]], vi = s[[2]]$vi, mods = s[[3]], method = "ML", data = s[[2]])
  cat(sprintf("%s full model: AIC FE = %.1f | AIC RE (ML) = %.1f | tau2 (ML) = %.4f\n",
              s[[1]], AIC(fe), AIC(ml), ml$tau2))
}
cat(sprintf("EP tau2 (REML) = %.4f\n", fitEP$REML$tau2))
modEP <- fitEP$REML
modEI <- fitEI$FE

# 1.2 Publication bias: Egger test and PET-PEESE
cat("\n== 1.2 Egger tests ==\n")
egEP <- regtest(modEP, predictor = "sei"); egEI <- regtest(modEI, predictor = "sei")
print(egEP); print(egEI)
write.csv(data.frame(Series = c("EP", "EI"), Model = c("Random effects (REML)", "Fixed effects"),
                     z = c(egEP$zval, egEI$zval), p = c(egEP$pval, egEI$pval)),
          file.path(out_dir, "S2_egger.csv"), row.names = FALSE)
png(file.path(out_dir, "funnel.png"), width = 1600, height = 800, res = 150)
par(mfrow = c(1, 2)); funnel(modEP, main = "EP"); funnel(modEI, main = "EI")
dev.off()

# PET regresses the estimate on its SE, PEESE on its variance; the
# intercept is the estimate purged of small-sample (publication) bias.
# PEESE is retained when the PET intercept is significant.
petpeese <- function(d, y) {
  res <- NULL
  for (m in c("FE", "REML", "Multi")) for (type in c("PET", "PEESE")) {
    mods <- if (type == "PET") ~ sqrt(vi) else ~ vi
    fit <- if (m == "Multi") rma.mv(yi = d[[y]], V = d$vi, mods = mods, random = ~ 1 | IDstudy, data = d)
           else rma(yi = d[[y]], vi = d$vi, mods = mods, method = m, data = d)
    cr <- coef_test(fit, vcov = "CR2", cluster = d$IDstudy)
    res <- rbind(res, data.frame(Series = y, Weights = m, Model = type,
                                 Estimate = fit$beta[1], SE_CR2 = cr$SE[1],
                                 CI_low = fit$beta[1] - 1.96 * cr$SE[1],
                                 CI_high = fit$beta[1] + 1.96 * cr$SE[1], p_CR2 = cr$p_Satt[1]))
  }
  res
}
pp <- rbind(petpeese(Epdata, "EP"), petpeese(Eidata, "EI"))
cat("\n== 1.2 PET-PEESE (retained: EP REML, EI FE) ==\n")
print(pp, row.names = FALSE, digits = 3)
write.csv(pp, file.path(out_dir, "S2_petpeese.csv"), row.names = FALSE)

# 1.3 CR2 cluster-robust moderators (Table 2; full tables in Table S1)
robust_EP <- coef_test(modEP, vcov = "CR2", cluster = Epdata$IDstudy)
robust_EI <- coef_test(modEI, vcov = "CR2", cluster = Eidata$IDstudy)
cat("\n== 1.3 Moderators significant at 10% (Table 2) ==\n")
print(robust_EP[robust_EP$p_Satt <= 0.10, ])
print(robust_EI[robust_EI$p_Satt <= 0.10, ])
write.csv(as.data.frame(robust_EP), file.path(out_dir, "S1_moderators_EP.csv"), row.names = FALSE)
write.csv(as.data.frame(robust_EI), file.path(out_dir, "S1_moderators_EI.csv"), row.names = FALSE)

cat("\n== 1.3 Joint CR2 Wald tests ==\n")
groups <- list(Concept = c("Longrun", "Compensated", "Conditional"),
               Datatype = "^Datatype", Time = c("Midpoint_c", "Duration"))
for (g in names(groups)) for (s in list(list("EP", modEP, Epdata), list("EI", modEI, Eidata))) {
  idx <- grep(paste(groups[[g]], collapse = "|"), names(coef(s[[2]])))
  w <- Wald_test(s[[2]], constraints = constrain_zero(idx),
                 vcov = "CR2", cluster = s[[3]]$IDstudy, test = "HTZ")
  cat(sprintf("  %-8s %s : F = %.2f, df = (%d, %.1f), p = %.3f\n",
              g, s[[1]], w$Fstat, w$df_num, w$df_denom, w$p_val))
}

# 1.4 Quadratic trend in the full models (LRT, retained model class)
cat("\n== 1.4 LRT on Midpoint^2 ==\n")
print(anova(rma(yi = EP, vi = vi, mods = update(formulaEP_full, ~ . + Midpoint_c2), method = "ML", data = Epdata),
            rma(yi = EP, vi = vi, mods = formulaEP_full, method = "ML", data = Epdata)))
print(anova(rma(yi = EI, vi = vi, mods = formulaEI_full, method = "FE", data = Eidata),
            rma(yi = EI, vi = vi, mods = update(formulaEI_full, ~ . - Midpoint_c2), method = "FE", data = Eidata)))

#=======================================================================
# STEP 2 — TEMPORAL HOMOGENEITY
#=======================================================================

# Pettitt: rank-based test for a single change in level, which also
# locates the most likely change point. Von Neumann: ratio of the mean
# squared successive difference to the variance (values below 2 signal
# a trend or persistence). Both accept irregularly spaced observations.
pettitt <- function(x) {
  n <- length(x)
  s <- sign(outer(x, x, "-"))
  U <- vapply(1:(n - 1), function(t) sum(s[1:t, (t + 1):n]), numeric(1))
  k <- which.max(abs(U)); K <- abs(U[k])
  list(U = K, k = k, p = min(1, 2 * exp(-6 * K^2 / (n^3 + n^2))), Uvec = abs(U))
}
vonneumann <- function(x) {
  n <- length(x)
  eta <- (sum(diff(x)^2) / (n - 1)) / (sum((x - mean(x))^2) / n)
  z <- (eta - 2 * n / (n - 1)) / sqrt(4 * n^2 * (n - 2) / ((n^2 - 1) * (n - 1)^2))
  list(eta = eta, p = 2 * pnorm(-abs(z)))
}

homog <- NULL
for (s in list(list("EP", Epdata), list("EI", Eidata))) {
  d <- s[[2]]
  subs <- list(full = d, `pre-1985` = d[d$Midpoint < 1985, ], `post-1985` = d[d$Midpoint >= 1985, ])
  for (nm in names(subs)) {
    x <- subs[[nm]][[s[[1]]]]; pt <- pettitt(x); vn <- vonneumann(x)
    homog <- rbind(homog, data.frame(
      Series = s[[1]], Period = nm, n = length(x),
      Pettitt_U = pt$U, Pettitt_p = pt$p, Change_point = subs[[nm]]$Midpoint[pt$k],
      VN_eta = vn$eta, VN_p = vn$p,
      Decision = ifelse(pt$p <= .05 & vn$p <= .05, "Rejected",
                 ifelse(pt$p <= .05 | vn$p <= .05, "Ambiguous", "Homogeneous"))))
  }
}
cat("\n== 2A Homogeneity tests (Table 1) ==\n")
print(homog, row.names = FALSE, digits = 3)
write.csv(homog, file.path(out_dir, "T1_homogeneity.csv"), row.names = FALSE)

# Pettitt profile. Several estimates share the same midpoint, and a split
# inside a group of ties has no temporal meaning: |U| is reported at each
# split between two distinct midpoints. A flat profile means the test dates
# the shift to a period rather than to a precise year.
pt_profile <- NULL
for (s in list(list("EP", Epdata), list("EI", Eidata))) {
  d <- s[[2]]; U <- pettitt(d[[s[[1]]]])$Uvec
  ends <- which(diff(d$Midpoint) > 0)
  pt_profile <- rbind(pt_profile, data.frame(
    Series = s[[1]], Split_after = d$Midpoint[ends], Next_midpoint = d$Midpoint[ends + 1],
    n_before = ends, n_after = nrow(d) - ends, Pettitt_absU = U[ends]))
}
cat("\n== 2A' Pettitt profile, EI splits 1975-1995 ==\n")
print(subset(pt_profile, Series == "EI" & Split_after >= 1975 & Split_after <= 1995), row.names = FALSE)
write.csv(pt_profile, file.path(out_dir, "S4_pettitt_profile.csv"), row.names = FALSE)

dist <- NULL
for (cut in c(1975, 1980, 1985, 1990)) for (s in list(list("EP", Epdata), list("EI", Eidata))) {
  a <- s[[2]][[s[[1]]]][s[[2]]$Midpoint < cut]; b <- s[[2]][[s[[1]]]][s[[2]]$Midpoint >= cut]
  dist <- rbind(dist, data.frame(
    Cutoff = cut, Series = s[[1]], n_before = length(a), n_after = length(b),
    Mean_before = mean(a), Mean_after = mean(b),
    Welch_p = t.test(a, b)$p.value,
    Wilcoxon_p = wilcox.test(a, b, exact = FALSE)$p.value,
    KS_p = suppressWarnings(ks.test(a, b)$p.value)))
}
cat("\n== 2B Distributional comparisons (Table S3) ==\n")
print(dist, row.names = FALSE, digits = 3)
write.csv(dist, file.path(out_dir, "S3_cutoff_comparisons.csv"), row.names = FALSE)

agg <- function(d, y) d %>% group_by(Midpoint) %>%
  summarise(y = mean(.data[[y]]), n = n(), .groups = "drop") %>% arrange(Midpoint)
cat("\n== 2C Bai-Perron and Chow (Table S4) ==\n")
bp_out <- NULL
for (s in list(list("EP", Epdata), list("EI", Eidata))) for (win in c("full", "1975+")) {
  a <- agg(if (win == "full") s[[2]] else s[[2]][s[[2]]$Midpoint >= 1975, ], s[[1]])
  bp <- breakpoints(y ~ 1, data = a, h = max(3, floor(0.15 * nrow(a))))
  m  <- which.min(BIC(bp)) - 1
  dates <- if (m > 0) paste(a$Midpoint[breakpoints(bp, breaks = m)$breakpoints], collapse = ", ") else "none"
  supF <- sctest(Fstats(y ~ 1, data = a, from = 0.15))
  chow <- sapply(c(1980, 1985, 1990), function(yr) {
    pt <- which(a$Midpoint >= yr)[1] - 1
    if (pt >= 3 && nrow(a) - pt >= 3) sctest(y ~ 1, type = "Chow", point = pt, data = a)$p.value else NA })
  bp_out <- rbind(bp_out, data.frame(Series = s[[1]], Window = win, n_dates = nrow(a),
                                     BIC_breaks = m, Break_dates = dates,
                                     supF = unname(supF$statistic), supF_p = supF$p.value,
                                     Chow1980_p = chow[1], Chow1985_p = chow[2], Chow1990_p = chow[3]))
}
print(bp_out, row.names = FALSE, digits = 3)
write.csv(bp_out, file.path(out_dir, "S4_baiperron_chow.csv"), row.names = FALSE)

fig <- function(d, y, col, ttl, ylab) ggplot(d, aes(Midpoint, .data[[y]])) +
  geom_point(alpha = .45, size = 1.8, colour = col) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), colour = "black", fill = "grey80") +
  geom_vline(xintercept = 1985, linetype = "dotted") +
  labs(title = ttl, x = "Median year of the estimation sample", y = ylab) +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "Figure1.png"),
       fig(Epdata, "EP", "steelblue", "(a) Own-price elasticities", "Own-price elasticity") /
       fig(Eidata, "EI", "coral", "(b) Income elasticities", "Income elasticity"),
       width = 7, height = 7, dpi = 300)

#=======================================================================
# STEP 3 — NET PERIOD EFFECT AFTER COMPOSITIONAL CONTROLS
#=======================================================================

controls <- list(
  EP = list(main      = ~ Country + Datatype + Method_grp,
            submitted = ~ Country + Functform + Method_grp,
            saturated = ~ Country + Datatype + Functform + Method_grp),
  EI = list(main      = ~ Country + Datatype + Method_grp,
            submitted = ~ Country + Functform,
            saturated = ~ Country + Datatype + Functform + Method_grp))
dat <- list(EP = Epdata, EI = Eidata)

fit_ref <- function(d, y, mods, m) {
  if (m == "FE") rma(yi = d[[y]], vi = d$vi, mods = mods, method = "FE", data = d)
  else rma(yi = d[[y]], vi = d$vi, mods = mods, method = "REML", test = "knha", data = d)
}
coef_row <- function(fit, d, term) {
  cr <- coef_test(fit, vcov = "CR2", cluster = d$IDstudy)
  i  <- which(rownames(fit$beta) == term)
  data.frame(Coef = cr$beta[i], SE = cr$SE[i], df_Satt = cr$df_Satt[i], p_Satt = cr$p_Satt[i])
}
period_test <- function(d, y, ref, m, cut) {
  d$post <- as.integer(d$Midpoint >= cut)
  cbind(Cutoff = cut, coef_row(fit_ref(d, y, update(ref, ~ . + post), m), d, "post"))
}
trend_test <- function(d, y, ref, m) {
  d$trend <- d$Midpoint_c / 10
  coef_row(fit_ref(d, y, update(ref, ~ . + trend), m), d, "trend")
}

cat("\n== 3.0 Weighting class of the reference models (ML-consistent AIC) ==\n")
classes <- list()
for (y in c("EP", "EI")) {
  d <- dat[[y]]; f <- controls[[y]]$main
  fe <- rma(yi = d[[y]], vi = d$vi, mods = f, method = "FE", data = d)
  ml <- rma(yi = d[[y]], vi = d$vi, mods = f, method = "ML", data = d)
  classes[[y]] <- if (AIC(ml) < AIC(fe)) "RE" else "FE"
  cat(sprintf("  %s: AIC FE = %.1f | AIC RE = %.1f | tau2 (ML) = %.4f -> %s\n",
              y, AIC(fe), AIC(ml), ml$tau2, classes[[y]]))
}

# 3A-3C Main control set: retained class, and the other class as check
t3 <- NULL; grid <- NULL
for (y in c("EP", "EI")) for (m in c("RE", "FE")) {
  d <- dat[[y]]; f <- controls[[y]]$main
  p85 <- period_test(d, y, f, m, 1985)
  ml  <- if (m == "RE") "ML" else "FE"
  lrt <- anova(rma(yi = d[[y]], vi = d$vi, mods = f, method = ml, data = d),
               rma(yi = d[[y]], vi = d$vi, mods = update(f, ~ . + post1985), method = ml, data = d))
  g   <- do.call(rbind, lapply(1970:2000, function(cut) period_test(d, y, f, m, cut)))
  tr  <- trend_test(d, y, f, m)
  grid <- rbind(grid, cbind(Series = y, Weights = m, g))
  t3 <- rbind(t3, data.frame(
    Series = y, Controls = "main", Weights = m, Retained = classes[[y]] == m,
    Post1985_coef = p85$Coef, Post1985_SE = p85$SE, Post1985_df = p85$df_Satt, Post1985_p = p85$p_Satt,
    LRT_chi2 = lrt$LRT, LRT_p = lrt$pval,
    Grid_sig05_of31 = sum(g$p_Satt < .05), Grid_min_p = min(g$p_Satt),
    Grid_median_coef = median(g$Coef),
    Trend_decade = tr$Coef, Trend_SE = tr$SE, Trend_p = tr$p_Satt))
}

# 3D Alternative control sets (retained weighting class)
for (y in c("EP", "EI")) for (cs in c("submitted", "saturated")) {
  d <- dat[[y]]; f <- controls[[y]][[cs]]; m <- classes[[y]]
  p85 <- period_test(d, y, f, m, 1985)
  g   <- do.call(rbind, lapply(1970:2000, function(cut) period_test(d, y, f, m, cut)))
  tr  <- trend_test(d, y, f, m)
  t3 <- rbind(t3, data.frame(
    Series = y, Controls = cs, Weights = m, Retained = FALSE,
    Post1985_coef = p85$Coef, Post1985_SE = p85$SE, Post1985_df = p85$df_Satt, Post1985_p = p85$p_Satt,
    LRT_chi2 = NA, LRT_p = NA,
    Grid_sig05_of31 = sum(g$p_Satt < .05), Grid_min_p = min(g$p_Satt),
    Grid_median_coef = median(g$Coef),
    Trend_decade = tr$Coef, Trend_SE = tr$SE, Trend_p = tr$p_Satt))
}
cat("\n== 3 Net period effect (Table 3; Table S5) ==\n")
print(t3, row.names = FALSE, digits = 3)
write.csv(t3, file.path(out_dir, "T3_net_period_effect.csv"), row.names = FALSE)
write.csv(grid, file.path(out_dir, "S5_cutoff_grid.csv"), row.names = FALSE)

g_ret <- do.call(rbind, lapply(c("EP", "EI"), function(y) grid[grid$Series == y & grid$Weights == classes[[y]], ]))
g_ret$half <- qt(0.975, g_ret$df_Satt) * g_ret$SE          # Satterthwaite-consistent 95% CI
g_ret$Significant <- ifelse(g_ret$p_Satt < 0.05, "p < 0.05", "p >= 0.05")
g_ret$Panel <- factor(ifelse(g_ret$Series == "EP", "(a) Own-price elasticity", "(b) Income elasticity"),
                      levels = c("(a) Own-price elasticity", "(b) Income elasticity"))
ggsave(file.path(out_dir, "FigureS1_cutoff_grid.png"),
       ggplot(g_ret, aes(Cutoff, Coef)) +
         geom_hline(yintercept = 0, colour = "grey40") +
         geom_vline(xintercept = 1985, linetype = "dotted") +
         geom_ribbon(aes(ymin = Coef - half, ymax = Coef + half), fill = "grey85") +
         geom_line(colour = "black") +
         geom_point(aes(shape = Significant), size = 1.8, fill = "black") +
         scale_shape_manual(values = c("p < 0.05" = 21, "p >= 0.05" = 1), name = NULL) +
         facet_wrap(~ Panel, scales = "free_y") +
         labs(x = "Cutoff year (estimates with sample midpoint >= cutoff = post period)",
              y = "Net post-cutoff effect (95% CI)") +
         theme_minimal(base_size = 11) + theme(legend.position = "bottom"),
       width = 8, height = 4, dpi = 300)

#=======================================================================
# STEP 4 — COMPOSITIONAL CHANGE PRE/POST 1985
#=======================================================================
cat("\n== 4 Composition pre/post 1985 (Table S6) ==\n")
comp <- NULL; shares <- NULL
for (y in c("EP", "EI")) {
  d <- dat[[y]]; per <- factor(ifelse(d$post1985 == 1, "post-1985", "pre-1985"),
                               levels = c("pre-1985", "post-1985"))
  for (v in c("Datatype", "Functform", "Method_grp", "Country", "Compensated", "Conditional", "Longrun")) {
    ft <- fisher.test(table(per, d[[v]]), simulate.p.value = TRUE, B = 1e4)
    comp <- rbind(comp, data.frame(Series = y, Variable = v, Fisher_p = ft$p.value))
  }
  for (v in c("Datatype", "Functform")) {
    tb <- as.data.frame.matrix(round(100 * prop.table(table(per, d[[v]]), 1), 1))
    shares <- rbind(shares, cbind(Series = y, Variable = v, Period = rownames(tb),
                                  Level = rep(colnames(tb), each = 2), Share = unlist(tb)))
  }
}
print(comp, row.names = FALSE, digits = 3)
write.csv(comp, file.path(out_dir, "S6_composition_tests.csv"), row.names = FALSE)
write.csv(shares, file.path(out_dir, "S6_composition_shares.csv"), row.names = FALSE)

#=======================================================================
# STEP 5 — ROBUSTNESS: ESTIMATES SIGNIFICANT AT 5% ONLY
#=======================================================================
cat("\n== 5 Significant-only subsample (Table S7) ==\n")
rob <- NULL
for (y in c("EP", "EI")) {
  se <- if (y == "EP") "SE_EP" else "SE_EI"
  d <- droplevels(dat[[y]][abs(dat[[y]][[y]] / dat[[y]][[se]]) >= 1.96, ])
  w <- if (y == "EP") "REML" else "FE"
  pe <- rma(yi = d[[y]], vi = d$vi, mods = ~ vi, method = w, data = d)
  pc <- coef_test(pe, vcov = "CR2", cluster = d$IDstudy)
  f  <- controls[[y]]$main; m <- classes[[y]]
  p85 <- period_test(d, y, f, m, 1985)
  g   <- do.call(rbind, lapply(1970:2000, function(cut) period_test(d, y, f, m, cut)))
  tr  <- trend_test(d, y, f, m)
  rob <- rbind(rob, data.frame(Series = y, n = nrow(d), studies = nlevels(d$IDstudy),
                               PEESE = pe$beta[1], PEESE_SE_CR2 = pc$SE[1],
                               Post1985_coef = p85$Coef, Post1985_p = p85$p_Satt,
                               Grid_sig05_of31 = sum(g$p_Satt < .05),
                               Trend_decade = tr$Coef, Trend_p = tr$p_Satt))
}
print(rob, row.names = FALSE, digits = 3)
write.csv(rob, file.path(out_dir, "S7_significant_only.csv"), row.names = FALSE)

cat("\nDone. Outputs in ./", out_dir, "\n", sep = "")
# END
