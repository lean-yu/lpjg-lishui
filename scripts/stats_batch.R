# stats_batch.R - statistics and figures from batch/results.csv (run after analyze_batch.py)
suppressMessages(library(ggplot2))
X <- file.path(Sys.getenv("HOME"), "lpjg_lishui", "batch")
d <- read.csv(file.path(X, "results.csv"))
d$kind <- factor(d$kind, levels = c("dry50", "dry50x2", "hot3", "comp"))
d$density_f <- factor(d$density)
d$age_f <- factor(d$age)
d$beta_f <- factor(d$beta)

# 1. Which factors matter? Share of total variance in event-year NPP change
fit <- lm(dNPP_event_pct ~ kind * density_f + density_f * mgmt + kind * age_f + beta_f, data = d)
a <- anova(fit)
cat("\n== Share of variance in event-year NPP change explained by each term (%) ==\n")
print(data.frame(term = rownames(a), share_pct = round(100 * a[["Sum Sq"]] / sum(a[["Sum Sq"]]), 1)),
      row.names = FALSE)

# 2. Density effect: change in event-year NPP per +1000 stems/ha, for every scenario x management x age x root depth
keys <- unique(d[, c("kind", "mgmt", "age", "beta")])
sl <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  k <- keys[i, ]
  s <- d[d$kind == k$kind & d$mgmt == k$mgmt & d$age == k$age & abs(d$beta - k$beta) < 1e-9, ]
  data.frame(kind = k$kind, mgmt = k$mgmt, age = k$age, beta = k$beta,
             slope = unname(coef(lm(dNPP_event_pct ~ I(density / 1000), data = s))[2]))
}))
w <- reshape(sl, idvar = c("kind", "mgmt", "age"), timevar = "beta", direction = "wide")
w <- w[order(w$kind, w$mgmt, w$age), ]
cat("\n== Density slope: event-year NPP change (percentage points) per +1000 stems/ha ==\n")
cat("   negative = denser stands lose more; one column per root_beta\n")
print(w, row.names = FALSE, digits = 3)

# 3. Carbon effects, root_beta 0.96
c96 <- d[abs(d$beta - 0.96) < 1e-9, ]
cat("\n== Mean cumulative NPP change over event + 5 years (kg C m-2), root_beta 0.96 ==\n")
print(round(tapply(c96$cum_dNPP, list(density = c96$density, mgmt = c96$mgmt), mean), 3))
cat("\n== Mean NPP change in the 3 years after the event (%), root_beta 0.96 ==\n")
print(round(tapply(c96$dNPP_post_pct, list(density = c96$density, mgmt = c96$mgmt), mean), 2))

# 4. Figure: centre line = root_beta 0.96, ribbon = range over 0.95-0.97
rng <- aggregate(dNPP_event_pct ~ kind + mgmt + age + density, data = d,
                 FUN = function(x) c(lo = min(x), hi = max(x)))
rng <- do.call(data.frame, rng)
p <- ggplot() +
  geom_ribbon(data = rng, aes(x = density, ymin = dNPP_event_pct.lo, ymax = dNPP_event_pct.hi,
                              fill = factor(age)), alpha = 0.15) +
  geom_line(data = c96, aes(x = density, y = dNPP_event_pct, colour = factor(age))) +
  geom_point(data = c96, aes(x = density, y = dNPP_event_pct, colour = factor(age)), size = 1.5) +
  geom_hline(yintercept = 0, colour = "grey50") +
  facet_grid(mgmt ~ kind, labeller = labeller(mgmt = c(A = "no thinning", B = "standard thinning"))) +
  labs(x = "Planting density (stems/ha)", y = "Event-year NPP change vs control (%)",
       colour = "Stand age", fill = "Stand age") +
  theme_bw()
ggsave(file.path(X, "fig_resistance_R.png"), p, width = 11, height = 5.5, dpi = 150)
cat("\nFigure:", file.path(X, "fig_resistance_R.png"), "\n")
