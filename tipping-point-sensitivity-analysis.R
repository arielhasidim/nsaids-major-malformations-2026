tipping_point_reclassify <- function(data,
                                     exposure,  # bare name, e.g., first_trimester_exposure_Ibuprofen
                                     outcome,  # bare name, e.g., Congenital_Malformations_Major
                                     covariates,  # character vector of adjustment terms
                                     condition,  # STRING condition to define the "pool" to reclassify from
                                     ratio    = 1,  # 0-1: יחס מתוך הקבוצה המועברת שיעמוד ב-condition
                                     steps    = seq(0, 0.03, by = 0.01),
                                     reps     = 10,
                                     family   = poisson(link = "log"),
                                     seed     = NULL) {
  # --- deps ---
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("marginaleffects", quietly = TRUE)
  
  # --- capture names ---
  expo <- deparse1(substitute(exposure))
  outc <- deparse1(substitute(outcome))
  
  # --- optional seed for reproducibility ---
  if (!is.null(seed))
    set.seed(seed)
  
  # --- build formula ---
  rhs <- c(expo, covariates)
  frm <- stats::as.formula(sprintf("%s ~ %s", outc, paste(rhs, collapse = " + ")))
  match_frm <- stats::as.formula(sprintf("%s ~ %s", expo, paste(covariates, collapse = " + ")))
  
  # --- validate ratio ---
  if (!is.numeric(ratio) ||
      length(ratio) != 1 || is.na(ratio) || ratio < 0 || ratio > 1) {
    stop("`ratio` must be a single numeric in [0,1].")
  }
  
  # --- evaluate condition (string) over data to define pools among UNEXPOSED ---
  cond_expr <- gsub("&&", "&", condition)
  cond_true_str  <- sprintf("(%s == 0) & (%s)", expo, cond_expr)
  cond_false_str <- sprintf("(%s == 0) & !(%s)", expo, cond_expr)
  
  # print(cond_expr)
  # print(cond_true_str)
  # print(cond_false_str)
  
  cond_true_idx  <- try(eval(str2lang(cond_true_str), envir = data), silent = TRUE)
  cond_false_idx <- try(eval(str2lang(cond_false_str), envir = data), silent = TRUE)
  
  # print(cond_true_idx)
  # print(cond_false_idx)
  # 
  # print(inherits(cond_true_idx, "try-error"))
  # print(!is.logical(cond_true_idx))
  # print(length(cond_true_idx)  != nrow(data))
  # print(inherits(cond_false_idx, "try-error"))
  # print(!is.logical(cond_false_idx))
  # print(length(cond_false_idx) != nrow(data)) 
  
  if (inherits(cond_true_idx, "try-error")  ||
      !is.logical(cond_true_idx)  ||
      length(cond_true_idx)  != nrow(data) ||
      inherits(cond_false_idx, "try-error") ||
      !is.logical(cond_false_idx) ||
      length(cond_false_idx) != nrow(data)) {
    stop(
      "`condition` must be a valid logical expression over `data` columns, provided as a STRING."
    )
  }
  base_pool_idx        <- which(cond_true_idx)         # נשמר לשם תאימות לאחור (pool=condition TRUE)
  base_pool_false_idx  <- which(cond_false_idx)        # unexposed & !condition
  
  N_pool_true  <- length(base_pool_idx)
  N_pool_false <- length(base_pool_false_idx)
  N_pool_all   <- N_pool_true + N_pool_false
  N_total <- nrow(data)
  
  # --- helper: set exposure to '1' preserving type as much as practical ---
  set_exposure_one <- function(vec, idx) {
    if (length(idx) == 0)
      return(vec)
    if (is.factor(vec)) {
      # try to set to level "1"; if not present, use the second level if exists
      lvl <- levels(vec)
      target <- if ("1" %in% lvl)
        "1"
      else if (length(lvl) >= 2)
        lvl[2]
      else
        lvl[1]
      vec[idx] <- target
      return(droplevels(vec))
    } else if (is.logical(vec)) {
      vec[idx] <- TRUE
      return(vec)
    } else {
      # numeric/integer/character -> coerce to 1 (character "1" if character)
      if (is.character(vec)) {
        vec[idx] <- "1"
      } else {
        vec[idx] <- 1
      }
      return(vec)
    }
  }
  
  # --- one run: reclassify k from pool, fit, extract metrics ---
  run_once <- function(p, r) {
    # כמה מעבירים בסה"כ (מתוך כל הלא-חשופים), וכמה מהם עם condition לפי ratio
    k_total <- floor(p * N_pool_all)
    if (k_total <= 0) {
      switch_idx <- integer(0)
    } else {
      k_true_desired  <- floor(ratio * k_total)
      k_true  <- min(k_true_desired, N_pool_true)
      k_false <- min(k_total - k_true, N_pool_false)
      # אם חסר כדי להשלים k_total – נמלא מהקבוצה השנייה לפי היתרה
      if (k_true + k_false < k_total) {
        extra <- k_total - (k_true + k_false)
        extra_from_true <- min(extra, N_pool_true - k_true)
        k_true  <- k_true + extra_from_true
        k_false <- k_total - k_true
      }
      idx_true  <- if (k_true  > 0)
        sample(base_pool_idx, size = k_true, replace = FALSE)
      else
        integer(0)
      idx_false <- if (k_false > 0)
        sample(base_pool_false_idx,
               size = k_false,
               replace = FALSE)
      else
        integer(0)
      switch_idx <- c(idx_true, idx_false)
    }
    
    df2 <- data
    df2[[expo]] <- set_exposure_one(df2[[expo]], switch_idx)
    
    # print(match_frm)
    
    m.out <- matchit(
      match_frm,
      data = df2,
      method = "quick",
      estimand = "ATE"
    )
    
    m.df2 <- match.data(m.out)
    
    # fit model
    fit <- stats::glm(frm, data = m.df2, family = family, weights = weights)
    # print(fit)
    
    # vcov_mat <- sandwich::vcovCL(fit, cluster = m.df2[["subclass"]])
    vcov_mat <- sandwich::vcovCL(fit, cluster = m.df2[, c("subclass", "MPatientIdentity")], type = "HC3")
    
    # avg_comparisons to get RR, robust SE (HC3)
    ac <- try(marginaleffects::avg_comparisons(
      fit,
      variables  = expo,
      # newdata = m.df2,
      type       = "response",
      comparison = "lnratioavg",
      transform = "exp",
      vcov       = vcov_mat
      # vcov       = ~MPatientIdentity + HC3
    ),
    silent = TRUE)
    # print(ac)
    
    if (inherits(ac, "try-error")) {
      return(
        dplyr::tibble(
          step  = p,
          rep = r,
          est   = NA_real_,
          p = NA_real_,
          `ci2.5` = NA_real_,
          `ci97.5` = NA_real_
        )
      )
    }
    
    # --- robust column extraction across marginaleffects versions ---
    get_col <- function(df, candidates) {
      hit <- candidates[candidates %in% names(df)]
      if (length(hit) == 0)
        return(rep(NA_real_, nrow(df)))
      return(df[[hit[1]]])
    }
    
    est_vec <- get_col(ac, c("Estimate", "estimate"))
    p_vec   <- get_col(ac, c("p", "p.value", "Pr(>|z|)"))
    lo_vec  <- get_col(ac, c("2.5 %", "conf.low"))
    hi_vec  <- get_col(ac, c("97.5 %", "conf.high"))
    
    # קח את הרשומה הראשונה (אם יש כמה קונטרסטים)
    dplyr::tibble(
      step  = p,
      rep   = r,
      est   = suppressWarnings(as.numeric(est_vec[1])),
      p     = suppressWarnings(as.numeric(p_vec[1])),
      `ci2.5`  = suppressWarnings(as.numeric(lo_vec[1])),
      `ci97.5` = suppressWarnings(as.numeric(hi_vec[1]))
    )
  }
  
  # --- iterate over steps and reps ---
  out <- lapply(steps, function(p) {
    lapply(seq_len(reps), function(r)
      run_once(p, r)) |> dplyr::bind_rows()
  }) |> dplyr::bind_rows()
  
  # --- attach attributes for bookkeeping ---
  attr(out, "N_total") <- N_total
  attr(out, "N_pool")  <- N_pool_all
  attr(out, "N_pool_condition_true")  <- N_pool_true
  attr(out, "N_pool_condition_false") <- N_pool_false
  attr(out, "sampling_ratio") <- ratio
  attr(out, "pool_condition_true_str")  <- cond_true_str
  attr(out, "pool_condition_false_str") <- cond_false_str
  attr(out, "exposure") <- expo
  attr(out, "outcome")  <- outc
  
  out
}

# קובריאטים למודל
covars <- c(
  "BirthTerminationYear",
  "I(BirthTerminationYear ^ 2)",
  "Mom_Age",
  "LOPC",
  "Bedouin",
  "NSAID_indication",
  "Diabetes",
  "PregnancyNum",
  "obesity",
  "folic_acid",
  "Smoker",
  "first_trimester_exposure_other_antipyretics"
)

cond <- "Congenital_Malformations_Major == 1"
res1 <- tipping_point_reclassify(
  data       = clean_df_excluded,
  exposure   = first_trimester_exposure_Ibuprofen,
  outcome    = Congenital_Malformations_Major,
  covariates = covars,
  condition  = cond,
  ratio    = 0.082,
  steps      = seq(0, 0.03, by = 0.0005),## by = 0.0005),
  reps       = 100, #100
  family     = poisson(link = "log")
)

res1_summary <-res1 %>% summarise(
  .by = step,
  rep = n(),
  est = mean(est, na.rm = TRUE),
  p = mean(p),
  ci2.5 = mean(ci2.5),
  ci97.5 = mean(ci97.5)
) 

res1_summary %>%
  ggplot(aes(x = step, y = est)) +
  geom_ribbon(aes(ymin = ci2.5, ymax = ci97.5), alpha = 0.2) +
  geom_line() +
  labs(
    x = "Proportion of re-classified cases",
    y = "Simulated RR",
    title = "Sensetivity analysis for tipping point and RR",
    subtitle = "Reclassification of unxposed cases in ratio of 20% major congenital malformations",
    caption = "proportion steps resolusion of 0.0005, 100 random samples per step"
  ) +
  ggsci::scale_color_jama	() +
  theme_classic()


x_star <- 0.013
y_star <- approx(res1_summary$step, res1_summary$est, xout = x_star, rule = 2)$y
xmin   <- min(res1_summary$step, na.rm = TRUE)
ymin   <- min(res1_summary$ci2.5, na.rm = TRUE)

x_star2 <- 0.02
y_star2 <- 1
xmin2   <- 0
ymin2   <- min(res1_summary$ci2.5, na.rm = TRUE)

res1_summary %>%
  ggplot(aes(x = step, y = est)) +
  geom_ribbon(aes(ymin = ci2.5, ymax = ci97.5), alpha = 0.15) +
  geom_line() +
  # קווים אדומים לצירים
  geom_segment(aes(x = x_star, xend = x_star, y = ymin,  yend = y_star),
               color = "red") +
  geom_segment(aes(x = xmin,   xend = x_star, y = y_star, yend = y_star),
               color = "red") +
  geom_segment(aes(x = x_star2, xend = x_star2, y = ymin2,  yend = y_star2),
               color = "brown",
               linetype="dashed") +
  geom_segment(aes(x = xmin2,   xend = x_star2, y = y_star2, yend = y_star2),
               color = "brown",
               linetype="dashed") +
  # geom_point(aes(x = x_star, y = y_star), color = "red", size = 2) + # נקודת החיתוך עצמה (אופציונלי)
  # טיקים + לייבלים אדומים בנקודות המפגש
  scale_x_continuous(
    breaks = c(0, x_star, 0.02, 0.03),
    labels = function(l) ifelse(abs(l - x_star) < 1e-12,
                                paste0("<span style='color:red;'>", scales::percent(l, accuracy = 0.01), "</span>"),
                                ifelse(abs(l - x_star2) < 1e-12,
                                       paste0("<span style='color:brown;'>", scales::percent(l, accuracy = 0.01), "</span>"),
                                       scales::percent(l, accuracy = 0.01)))
  ) +
  scale_y_continuous(
    breaks = c(0.90, 1, y_star, 1.15),
    labels = function(l) ifelse(abs(l - y_star) < 1e-12,
                                paste0("<span style='color:red;'>", scales::number(l, accuracy = 0.01), "</span>"),
                                ifelse(abs(l - y_star2) < 1e-12,
                                       paste0("<span style='color:brown;'>", scales::number(l, accuracy = 0.01), "</span>"),
                                       scales::number(l, accuracy = 0.01)))
  ) +
  coord_cartesian(xlim = c(0.001, 0.03), ylim = c(0.91, 1.15)) +
  labs(
    x = "Proportion of reassigned cases",
    y = "Simulated RR",
    # title = "Sensitivity analysis for tipping point and RR",
    # subtitle = "Reclassification of 'Congenital_Malformations_Major == 1' in proportions of 20%",
    # caption = "proportion steps resolution of 0.0005, 100 random samples per step"
  ) +
  theme_classic() +
  theme(
    axis.text.x = ggtext::element_markdown(),
    axis.text.y = ggtext::element_markdown()
  )


# save(res1, file = "sensi-4-2-2026.Rda")
