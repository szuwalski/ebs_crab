library(PBSmodelling)
library(dplyr)
library(ggplot2)
library(mgcv)

outs <- readList("models/bbrkc/test/rkc.rep")
all_output <- read.csv("data/all_output.csv")
alt_metrics_calc <- read.csv("data/alt_metrics_calc.csv")
unc_mort <- read.csv("data/uncertainty_mort.csv")

m_dat <- filter(unc_mort, stock == "BBRKC")
alt_d <- filter(alt_metrics_calc, stock == "BBRKC")
dens  <- filter(all_output, species == "BBRKC" & process == "Abundance")
colnames(dens)[2] <- "Abund"

gam_dat <- merge(m_dat, alt_d, by = "Year")
gam_dat <- merge(gam_dat, dens, by = "Year")

gam_dat <- gam_dat %>%
  dplyr::select(Year, est_m, Temperature, Size, Abund, sd) %>%
  filter(
    is.finite(Year),
    is.finite(est_m),
    is.finite(Size),
    is.finite(Abund),
    is.finite(sd),
    sd > 0
  ) %>%
  mutate(m_resp = 1 - exp(-est_m))

eps <- 1e-6
gam_dat$m_resp <- pmin(pmax(gam_dat$m_resp, eps), 1 - eps)

as_num_matrix <- function(x, ncol) {
  x <- as.numeric(x)
  matrix(x, ncol = ncol, byrow = FALSE)
}

styr  <- outs$styr
endyr <- outs$endyr
sizes <- outs$sizes

hist_yrs <- seq(styr, endyr)
size_n   <- length(sizes)

pred_numbers_mat <- as_num_matrix(outs$`pred numbers at size`, ncol = size_n)
nat_m_hist       <- as_num_matrix(outs$`natural mortality`, ncol = size_n)
ret_size_comp    <- as_num_matrix(outs$`pred_retained_size_comp`, ncol = size_n)
ret_size_comp[!is.finite(ret_size_comp)] <- 0

n_hist <- sweep(pred_numbers_mat, 1, outs$numbers_pred, "*")

retain_fish_sel  <- outs$ret_fish_sel
total_fish_sel   <- outs$total_fish_sel
size_trans       <- outs$size_trans
in_prob_molt     <- outs$in_prob_molt
temp_prop_rec    <- outs$temp_prop_rec
recruits         <- outs$recruits
est_fishing_mort <- outs$est_fishing_mort

advance_one_year <- function(n_start,
                             nat_m_y,
                             f_y,
                             discard_survival = 0.85) {
  temp_n <- n_start * exp(-1 * 0.17 * exp(nat_m_y))
  
  temp_catch_n <- temp_n * (1 - exp(-(f_y * total_fish_sel)))
  temp_n <- temp_n * exp(-(f_y * total_fish_sel))
  temp_n <- temp_n + temp_catch_n * (1 - retain_fish_sel * discard_survival)
  
  trans_n <- size_trans %*% (temp_n * in_prob_molt)
  temp_n <- trans_n + (temp_n * (1 - in_prob_molt))
  
  temp_n[1] <- temp_n[1] + median(recruits, na.rm = TRUE) * temp_prop_rec[1]
  temp_n[2] <- temp_n[2] + median(recruits, na.rm = TRUE) * temp_prop_rec[2]
  temp_n[3] <- temp_n[3] + median(recruits, na.rm = TRUE) * temp_prop_rec[3]
  
  n_next <- temp_n * exp(-1 * 0.83 * exp(nat_m_y))
  as.numeric(n_next)
}

hindcast_nat_m <- function(start_fit_year = 2014,
                           end_proj_year = 2025,
                           discard_survival = 0.85,
                           k_smooth = 4,
                           use_realized_f = TRUE,
                           constant_f = NULL,
                           carry_temp_year = 2020) {
  
  all_years <- seq(styr, end_proj_year)
  ny <- length(all_years)
  
  year_to_row <- function(y) match(y, all_years)
  
  n_size_pred <- matrix(NA_real_, nrow = ny, ncol = size_n)
  nat_m_used  <- matrix(NA_real_, nrow = ny, ncol = size_n)
  f_mort      <- rep(NA_real_, ny)
  
  hist_fill_years <- seq(styr, start_fit_year)
  hist_fill_rows  <- year_to_row(hist_fill_years)
  
  n_size_pred[hist_fill_rows, ] <- n_hist[hist_fill_rows, ]
  nat_m_used[hist_fill_rows, ]  <- nat_m_hist[hist_fill_rows, ]
  
  if (use_realized_f) {
    f_mort[hist_fill_rows] <- est_fishing_mort[hist_fill_rows]
  } else {
    if (is.null(constant_f)) {
      constant_f <- median(est_fishing_mort[hist_fill_rows], na.rm = TRUE)
    }
    f_mort[hist_fill_rows] <- constant_f
  }
  
  first_hc_year <- start_fit_year + 1
  row_fit_year  <- year_to_row(start_fit_year)
  row_first_hc  <- year_to_row(first_hc_year)
  
  f_y <- if (use_realized_f) est_fishing_mort[row_fit_year] else constant_f
  
  n_size_pred[row_first_hc, ] <- advance_one_year(
    n_start = n_size_pred[row_fit_year, ],
    nat_m_y = nat_m_hist[row_fit_year, ],
    f_y = f_y,
    discard_survival = discard_survival
  )
  
  pred_years <- seq(first_hc_year, end_proj_year)
  
  annual_pred <- data.frame(
    Year = pred_years,
    fit_through = pred_years - 1,
    Temperature = NA_real_,
    proj_Size = NA_real_,
    proj_Abund = NA_real_,
    pred_m_resp = NA_real_,
    pred_est_m = NA_real_,
    actual_est_m = gam_dat$est_m[match(pred_years, gam_dat$Year)],
    prediction_type = NA_character_
  )
  
  for (yy in pred_years) {
    row_y <- year_to_row(yy)
    
    tmp_temp <- gam_dat$Temperature[match(yy, gam_dat$Year)]
    tmp_abund <- sum(n_size_pred[row_y, ], na.rm = TRUE) / sum(n_size_pred[1, ], na.rm = TRUE)
    tmp_size  <- weighted.mean(sizes, w = n_size_pred[row_y, ], na.rm = TRUE)
    
    annual_pred$Temperature[annual_pred$Year == yy] <- tmp_temp
    annual_pred$proj_Size[annual_pred$Year == yy]   <- tmp_size
    annual_pred$proj_Abund[annual_pred$Year == yy]  <- tmp_abund
    
    # special carry-forward year: use previous year's predicted M
    if (!is.finite(tmp_temp) && yy == carry_temp_year) {
      prev_pred <- annual_pred$pred_est_m[annual_pred$Year == (yy - 1)]
      pred_est_m <- prev_pred
      pred_resp  <- 1 - exp(-pred_est_m)
      
      annual_pred$prediction_type[annual_pred$Year == yy] <- "carry_forward"
    } else {
      train_dat <- gam_dat %>%
        filter(
          Year <= (yy - 1),
          is.finite(Temperature)
        )
      
      m_mod <- gam(
        m_resp ~ s(Temperature, k = k_smooth) +
          s(Size, k = k_smooth) +
          s(Abund, k = k_smooth),
        data = train_dat,
        family = betar(link = "logit"),
        weights = 1 / sd,
        method = "REML"
      )
      
      pred_resp <- predict(
        m_mod,
        newdata = data.frame(
          Temperature = tmp_temp,
          Size = tmp_size,
          Abund = tmp_abund
        ),
        type = "response"
      )
      
      pred_resp  <- pmin(pmax(as.numeric(pred_resp), eps), 1 - eps)
      pred_est_m <- -log(1 - pred_resp)
      
      annual_pred$prediction_type[annual_pred$Year == yy] <- "gam"
    }
    
    nat_m_used[row_y, ] <- pred_est_m
    
    annual_pred$pred_m_resp[annual_pred$Year == yy] <- pred_resp
    annual_pred$pred_est_m[annual_pred$Year == yy]  <- pred_est_m
    
    if (yy < end_proj_year) {
      f_y <- if (use_realized_f) est_fishing_mort[row_y] else constant_f
      f_mort[row_y] <- f_y
      
      n_size_pred[row_y + 1, ] <- advance_one_year(
        n_start = n_size_pred[row_y, ],
        nat_m_y = nat_m_used[row_y, ],
        f_y = f_y,
        discard_survival = discard_survival
      )
    }
  }
  
  list(
    years = all_years,
    n_size_pred = n_size_pred,
    nat_m_used = nat_m_used,
    annual_pred = annual_pred
  )
}

hc <- hindcast_nat_m(
  start_fit_year = 2014,
  end_proj_year = 2025,
  discard_survival = 0.85,
  k_smooth = 4,
  use_realized_f = TRUE,
  carry_temp_year = 2020
)

plot_dat <- gam_dat %>%
  dplyr::select(Year, est_m) %>%
  rename(actual_est_m = est_m) %>%
  left_join(
    hc$annual_pred %>%
      dplyr::select(Year, pred_est_m, prediction_type),
    by = "Year"
  )

ggplot(plot_dat, aes(x = Year)) +
  geom_line(aes(y = actual_est_m), color = "black", linewidth = 1.1) +
  geom_point(aes(y = actual_est_m), color = "black", size = 2) +
  geom_line(aes(y = pred_est_m), color = "firebrick", linewidth = 1.1) +
  geom_point(aes(y = pred_est_m, shape = prediction_type), color = "firebrick", size = 2.4) +
  geom_vline(xintercept = 2014.5, linetype = 2, color = "grey40") +
  theme_bw() +
  labs(
    x = "Year",
    y = "Natural mortality (est_m scale)",
    title = "BBRKC rolling hindcast of mortality",
    subtitle = "Black = assessment estimates; red = hindcast; 2020 carried forward from 2019, then GAM resumes",
    shape = "Prediction"
  )