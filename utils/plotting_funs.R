library(rstan)
library(matrixStats)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(bayesplot)
library(posterior)
library(ftplottools)
library(ggsci)
library(cowplot)
library(ggExtra)
source(here("two_strain_model", "utils", "geom-stepribbon.R"))

plt_cross <- function(a) {
  mcmc_scatter(a$fit, pars = c("R_difference","cross"), size = 3.5, alpha = 0.1) +
    ylab("Cross-immnunity") +
    xlab("Transmissibility increase") +
    ggtitle(deparse(substitute(a))) +
    theme_pubr(base_size = 26) +
    stat_density_2d(color = "black", size = .5) 
}

plt_rr <- function(a) {
  mcmc_scatter(a$fit, pars = c("RR[1]","cross"), size = 3.5, alpha = 0.1) +
    ylab("Cross-immunity") +
    xlab("Relative risk of mortality") +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(deparse(substitute(a))) +
    theme_pubr(base_size = 26) +
    #     xlim(c(0,2)) +
    stat_density_2d(color = "black", size = .5)
}

plt_ifr <- function(a) {
  mcmc_scatter(a$fit, pars = c("ifr1[1]","ifr2[1]"), size = 3.5, alpha = 0.1) +
    xlab("ifr1") +
    ylab("ifr2") +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(deparse(substitute(a))) +
    theme_pubr(base_size = 26) +
    #     xlim(c(0,2)) +
    stat_density_2d(color = "black", size = .5)
}

plot_deaths_rt_immunity <- function(data){
  stan_data <- data$stan_data
  dates <- data$dates
  selected_region <- data$selected_region
  out <- rstan::extract(data$fit)
  ### plotting deaths
  deaths <- stan_data$deaths[,1]
  deaths[which(deaths<0 | deaths > 250)] <- NA
  
  deaths = data.frame(
    "time" = data$dates[[1]],
    "deaths" =deaths
  )
  
  deaths_df95v1 = data.frame(
    "time" = data$dates[[1]],
    "key" = rep("95v1", length(data$dates[[1]])),
    "deaths" = colMeans(out$E_deaths),
    "deaths" = colMeans(out$E_deaths_v1),
    "deaths_l" = colQuantiles(out$E_deaths_v1[,,], probs=.025),
    "deaths_u" = colQuantiles(out$E_deaths_v1[,,], probs=.975)
  )
  
  deaths_df95v2 = data.frame(
    "time" = data$dates[[1]],
    "key" = rep("95v2", length(data$dates[[1]])),
    "deaths" = colMeans(out$E_deaths),
    "deaths" = colMeans(out$E_deaths_v1),
    "deaths_l" = colQuantiles(out$E_deaths_v2[,,], probs=.025),
    "deaths_u" = colQuantiles(out$E_deaths_v2[,,], probs=.975)
  )
  
  deaths_df50v1 = data.frame(
    "time" = data$dates[[1]],
    "key" = rep("50v1", length(data$dates[[1]])),
    "deaths" = colMeans(out$E_deaths),
    "deaths" = colMeans(out$E_deaths_v1),
    "deaths_l" = colQuantiles(out$E_deaths_v1[,,], probs=.25),
    "deaths_u" = colQuantiles(out$E_deaths_v1[,,], probs=.75)
  )
  
  deaths_df50v2 = data.frame(
    "time" = data$dates[[1]],
    "key" = rep("50v2", length(data$dates[[1]])),
    "deaths" = colMeans(out$E_deaths),
    "deaths" = colMeans(out$E_deaths_v1),
    "deaths_l" = colQuantiles(out$E_deaths_v2[,,], probs=.25),
    "deaths_u" = colQuantiles(out$E_deaths_v2[,,], probs=.75)
  )
  
  deaths_df = rbind(deaths_df95v1,deaths_df95v2,deaths_df50v1,deaths_df50v2)
  g1 <- 
    ggplot() + 
    geom_blank() + 
    theme_bw() +
    geom_point(data = deaths, 
               aes(x = time, y = deaths))+
    geom_ribbon(data = deaths_df, 
                aes(x=time, ymax=deaths_l, ymin=deaths_u, fill = key)) +
    scale_fill_manual(name = "", labels = c("50%","50%","95%","95%"),
                      values = c(alpha("Deepskyblue4",0.85),
                                 alpha("chocolate1",0.85),
                                 alpha("Deepskyblue4",0.45),
                                 alpha("chocolate1",0.45))) +
    theme_light()+
    theme(strip.text.x = element_text(
      color = "black"),
      strip.background = element_blank()
    ) +  # guides(color = guide_legend(override.aes = list(size = c(2,1,1)))) +
    scale_y_continuous(name = "Deaths")+
    scale_x_date(name = "Date", date_labels = "%b %y", date_breaks = "1 month") +
    ft_theme() + 
    scale_color_npg() + 
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 24)
    )
  
  ### plotting Rt
  g2 <- 
    ggplot() + 
    geom_blank() + 
    theme_bw() +
    geom_line(aes(x=dates[[selected_region]],
                  y = colMedians(out$Rt_adj_immune_v1[,,1])),
              alpha = 0.3,
              color="seagreen"
    ) +
    geom_stepribbon(aes(x=dates[[selected_region]],
                        ymin = colQuantiles(out$Rt_adj_immune_v1[,,1], probs = 0.025), 
                        ymax = colQuantiles(out$Rt_adj_immune_v1[,,1], probs = 0.975)),
                    alpha = 0.3,
                    fill="seagreen"
    ) +
    geom_line(aes(x=dates[[selected_region]],
                  y = colMedians(out$Rt_adj_immune_v2[,,1])),
              alpha = 0.3,
              color="chocolate1"
    ) +
    geom_stepribbon(aes(x=dates[[selected_region]],
                        ymin = colQuantiles(out$Rt_adj_immune_v2[,,1], probs = 0.025), 
                        ymax = colQuantiles(out$Rt_adj_immune_v2[,,1], probs = 0.975)),
                    alpha = 0.3,
                    fill="chocolate1"
    ) +
    theme_light()+
    theme(strip.text.x = element_text(
      color = "black"),
      strip.background = element_blank()
    ) +  # guides(color = guide_legend(override.aes = list(size = c(2,1,1)))) +
    scale_y_continuous(name = "Rt")+
    scale_x_date(name = "Date", date_labels = "%b %y", date_breaks = "1 month") +
    ft_theme() + 
    scale_color_npg() + 
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 24)
    )
  
  posterior <- as.matrix(data$fit)
  posterior[,'cross'] <- (1 - posterior[,'cross'])
  x <- data.frame(R_difference = out$R_difference, cross =1-out$cross)
  kd <- ks::kde(x, compute.cont=TRUE,H=matrix(c(0.01,0,0,0.0033),2))
  contour_95 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                                 z=estimate, levels=cont["5%"])[[1]]))
  contour_75 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                                 z=estimate, levels=cont["25%"])[[1]]))
  contour_50 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                                 z=estimate, levels=cont["50%"])[[1]]))
  contour_25 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                                 z=estimate, levels=cont["75%"])[[1]]))
  contour_5 <- data.frame(with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                                z=estimate, levels=cont["95%"])[[1]]))
  pA <- mcmc_scatter(posterior, pars = c("R_difference","cross"), size = 3.5, alpha = 0.1) +
    ylab("Immune evasion (%)") +
    xlab("Transmissibility increase") +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw(base_size = 26) +
    ggplot2::geom_path(aes(x, y), data=contour_95) +
    ggplot2::geom_path(aes(x, y), data=contour_75) +
    ggplot2::geom_path(aes(x, y), data=contour_50) +
    ggplot2::geom_path(aes(x, y), data=contour_25) +
    ggplot2::geom_path(aes(x, y), data=contour_5) +
    theme_light()+
    theme(strip.text.x = element_text(
      color = "black"),
      strip.background = element_blank()
    ) +  # guides(color = guide_legend(override.aes = list(size = c(2,1,1)))) +
    ft_theme() + 
    # scale_y_continuous(breaks = c(0, 0.25,0.5,0.75 ,1.1), labels = c("0%", "25%", "50%", "75%", "100%")) +
    # scale_y_continuous(limits = c(0,1.1), labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(expand=c(0.1,0),labels = function(x) paste0(x*100, "%"))  + 
    scale_color_npg() + 
    theme(axis.text.x = element_text(hjust = 1),
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 24)
    ) 
  pA <- 
    ggMarginal(pA, R_difference, cross, type = c("density"),
               margins = 'both',
               size = 4,
               colour = "#589e73",
               fill = '#589e73',
               alpha = 0.25,
               xparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01)),
               yparams = list(colour = "#589e73", size = 1,bw=sqrt(0.01))) 
  p <- plot_grid(plot_grid(g1,g2, ncol = 1, labels=c("A", "B")),pA,ncol = 2, labels= c("", "C"))
  title <- ggdraw() + draw_label(deparse(substitute(data)), fontface='bold', size = "24") + ft_theme()
  # return(   plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) )
  return(p)
  
}

plt_sero_conv <- function(data)
{
  stan_data <- data$stan_data
  dates <- data$dates
  selected_region <- data$selected_region
  out <- rstan::extract(data$fit)
  ### plotting deaths
  deaths <- stan_data$deaths[,1]
  deaths[which(deaths<0 | deaths > 250)] <- NA
  mumbaiPopulation = 2219580
  seroconv_df95v1 = data.frame(
    "time" = dates[[selected_region]],
    "key" = rep("95v1", length(dates[[selected_region]])),
    "seroconv" = colMeans(out$seroconv_v1)/stan_data$pop[[1]],
    "seroconv_l" = colQuantiles(out$seroconv_v1[,,1], probs=.025)/stan_data$pop[[1]],
    "seroconv_u" = colQuantiles(out$seroconv_v1[,,1], probs=.975)/stan_data$pop[[1]]
  )
  seroconv_df95v2 = data.frame(
    "time" = dates[[selected_region]],
    "key" = rep("95v2", length(dates[[selected_region]])),
    "seroconv" = colMeans(out$seroconv_v1)/stan_data$pop[[1]],
    "seroconv_l" = colQuantiles(out$seroconv_v2[,,1], probs=.025)/stan_data$pop[[1]],
    "seroconv_u" = colQuantiles(out$seroconv_v2[,,1], probs=.975)/stan_data$pop[[1]]
  )
  seroconv_df50v1 = data.frame(
    "time" = dates[[selected_region]],
    "key" = rep("50v1", length(dates[[selected_region]])),
    "seroconv" = colMeans(out$seroconv_v1)/stan_data$pop[[1]],
    "seroconv_l" = colQuantiles(out$seroconv_v1[,,1], probs=.25)/stan_data$pop[[1]],
    "seroconv_u" = colQuantiles(out$seroconv_v1[,,1], probs=.75)/stan_data$pop[[1]]
  )
  seroconv_df50v2 = data.frame(
    "time" = dates[[selected_region]],
    "key" = rep("50v2", length(dates[[selected_region]])),
    "seroconv" = colMeans(out$seroconv_v2)/stan_data$pop[[1]],
    "seroconv_l" = colQuantiles(out$seroconv_v2[,,1], probs=.25)/stan_data$pop[[1]],
    "seroconv_u" = colQuantiles(out$seroconv_v2[,,1], probs=.75)/stan_data$pop[[1]]
  )
  seroconv_df = rbind(seroconv_df95v1,seroconv_df95v2,seroconv_df50v1,seroconv_df50v2)
  
  mumbai_sero <- tibble(dates = dates$Mumbai[stan_data$sero_N],
                        sero_prev = stan_data$sero_prev/100)
  p1 <-
    seroconv_df %>% 
    ggplot() +
    geom_ribbon(aes(x = time, ymin = seroconv_l, ymax = seroconv_u, fill =key)) +
    scale_fill_manual(name = "", labels = c("50%", "95%","50%", "95%"),
                      values = c(alpha("Deepskyblue4",0.85),
                                 alpha("tan2",0.85),
                                 alpha("Deepskyblue4",0.55),
                                 alpha("tan2",0.55))) +
    geom_point(data = mumbai_sero, aes(x=dates,y=sero_prev)) +
    xlab("") +
    ylab("Cumulative incidence per capita\n") +
    theme_light()+
    theme(strip.text.x = element_text(
      color = "black"),
      strip.background = element_blank()
    ) +  # guides(color = guide_legend(override.aes = list(size = c(2,1,1)))) +
    scale_y_continuous(name = "Deaths")+
    scale_x_date(name = "Date", date_labels = "%b %y", date_breaks = "1 month") +
    ft_theme() + 
    scale_color_npg() + 
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 24))
  
  p1
}

posterior_interest <- function(data, UR = "30%", VDATE="2020-12-01"){
  out <- rstan::extract(data$fit)
  ######################################################### 
  posteriors_of_interest <- 
    bind_rows(
      tibble(parameter = "Transmissibility Increase(%)",
             `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(out$R_difference - 1,probs = c(0.25)),0), round(100*quantile(out$R_difference - 1,probs = c(0.75)),0)),
             `value(95% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(out$R_difference - 1,probs = c(0.025)),0), round(100*quantile(out$R_difference - 1,probs = c(0.975)),0)),
             median = round(100 * median(out$R_difference - 1))
      ),
      tibble(parameter = "Relative Risk Increase(%)",
             `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(out$RR - 1,probs = c(0.25)),0), round(100*quantile(out$RR - 1,probs = c(0.75)),0)),
             `value(95% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(out$RR - 1,probs = c(0.025)),0), round(100*quantile(out$RR - 1,probs = c(0.975)),0)),
             median = round(100 * median(out$RR - 1))
      ),
      tibble(parameter = "Immunity Evasion(%)",
             `value(50% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(1 - out$cross,probs = c(0.25)),0), round(100*quantile( 1 - out$cross,probs = c(0.75)),0)),
             `value(95% CI)` = sprintf("%1.0f%% - %1.0f%%", round(100*quantile(1 - out$cross,probs = c(0.025)),0), round(100*quantile( 1 - out$cross,probs = c(0.975)),0)),
             median = round(100 * median(1 - out$cross))
      )
    )
  as_hux(posteriors_of_interest) %>%
    set_all_padding(4) %>% 
    set_outer_padding(0) %>% 
    set_number_format(2) %>% 
    set_bold(row = 1, col = everywhere) %>% 
    set_bottom_border(row = 1, col = everywhere) %>% 
    set_width(0.4) %>% 
    set_caption(glue("Epidemiological Properties of B1671.2 (UR = {UR}, VDATE= {VDATE})",))
}
