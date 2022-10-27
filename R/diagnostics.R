#' .test.ph.aft
#'
#' Test the PH and AFT assumptions
#'
#' Uses IPD/pseudo IPD to fit a cox model and visualize the results of the
#' cox.zph test including p-value. Also fits a log cumulative hazard plot
#' to visualize PH.
#'
#' @param  guyot_data pseudo-ipd with the format as output by simsurv function
#'   within this project.
#'
#' @return  The results of cox.zph, a nicer version of the cox.zph plot,
#'   a plot of log cumulative hazard against log time, and a plot of quantiles
#'   of each arm across time.
#'
#' @examples
#'  Put some examples of how it works or just delete this
#'
#' @export
.test.ph.aft <- function(dat){


  ############################################################################## #
  ############################################################################## #
  #
  # 1. Input error catches----
  #
  #     Section Notes
  #
  ############################################################################## #
  ############################################################################## #


    ############################################################################## #
    ############################################################################## #
    #
    # 2. Cox zph----
    #
    #     Section Notes
    ############################################################################## #
    ############################################################################## #

    cox <- survival::coxph(survival::Surv(futime, fustat) ~ trt, data = dat)
    zph <- survival::cox.zph(cox)
    p <- zph$table[,"p"][["trt"]]

    p <- ifelse(p < 0.001, "< 0.001", round(p,3))

    zph.plot <- data.frame(y = zph$y[,1],
                           x = zph$time) %>%
      ggplot2::ggplot(ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(size = 0.5, colour = cols[[1]]) +
      ggplot2::geom_smooth(se = TRUE, colour = cols[[1]], fill = cols[[2]]) +
      ggplot2::theme_minimal(base_size = 18) +
      ggplot2::labs(x = "Time", y = "Beta(t) for Treatment",
                    caption = glue::glue("Test for violation of proportional hazards = {p}"),
                    title = glue::glue("Hazard Ratio over Time"))


    mods <- c( "exp","weibull", "weibullPH", "gompertz",
              "gamma", "llogis", "lnorm")



    joint <- purrr::map(mods, ~ {
      dist <- .
      flexsurv::flexsurvreg(survival::Surv(futime, fustat) ~ trt,dist = ., data = dat)
    }) %>%
      purrr::set_names(glue::glue("{mods} - Joint"))


    sep <- purrr::map(mods[-1], ~ {
      dist <- .

      if(!dist %in% c("lnorm")){
        flexsurv::flexsurvreg(survival::Surv(futime, fustat) ~ trt,
                              anc = list("shape" = ~ trt),
                              dist = dist, data = dat)
      } else{
        flexsurv::flexsurvreg(survival::Surv(futime, fustat) ~ trt,
                              anc = list("sdlog" = ~ trt),
                              dist = dist, data = dat)
      }


    }) %>%
      purrr::set_names(glue::glue("{mods[-1]} - Separate"))

    all.mods <- c(joint, sep)

    aic_comp <- tibble::tibble(Model = names(all.mods),
                               AIC = purrr::map_dbl(all.mods, AIC),
                               BIC = purrr::map_dbl(all.mods, BIC)) %>%
      dplyr::arrange(AIC) %>%
      dplyr::mutate(AICrank = 1:dplyr::n()) %>%
      dplyr::arrange(BIC) %>%
      dplyr::mutate(BICrank = 1:dplyr::n()) %>%
      dplyr::arrange(AIC)

    out <- dplyr::lst(fits = aic_comp,
                      zph.plot,
                      all.mods)

    return(out)
    }




#' Function name
#'
#' What does the function do in one sentence
#'
#' More detailed explanation of assumptions etc
#'
#' @param  First argument (repeat for more)
#'
#' @return  What does the function return?
#'
#' @examples
#'  Put some examples of how it works or just delete this
#'
#' @export
.km.plots <- function(fits,dat,include){


  ############################################################################## #
  ############################################################################## #
  #
  # 3.  Cumulative Hazard and KM Curves ----
  #
  #     Section Notes
  #
  ############################################################################## #
  ############################################################################## #
  cols <- c("#002F6C", "#ED8B00")
  m <- survival::survfit(survival::Surv(futime, fustat) ~ trt, data = dat)

  p <- survminer::ggsurvplot(m, conf.int = TRUE, risk.table = TRUE,
                             ggtheme = ggplot2::theme_minimal(base_size = 18),
                             palette = cols,
                             surv.median.line = "v",
                             pval = TRUE,
                             data = dat)

  p$plot <- p$plot +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_color_manual(values = cols, name = "Treatment",
                                labels = c("Control", "Treatment")) +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(x = ggplot2::element_blank(),
                  title = glue::glue("KM Curve"))

  p$table <- p$table + ggplot2::labs(y = "Treatment") +
    ggplot2::scale_y_discrete(labels = c("Treatment", "Control"))

  tt <- data.frame(times = round(as.numeric(unique(p$plot$data$time))))


  pred <- function(fit, model){

    p <-  predict(fit, times = seq(from = 0.01, to = max(p$plot$data$time), length.out = 30), newdata  = data.frame(trt = c(0, 1)), type = "survival")

    rbind(p[[1]][[1]] %>% dplyr::mutate(trt = 0),
          p[[1]][[2]] %>% dplyr::mutate(trt = 1)) %>%
      dplyr::rename(time = .time,
                    surv = .pred_survival) %>%
      dplyr::mutate(model = model)
  }
    #---------------------------------------------------------------------- -
    # ..3b. Plot parametric fits over survivor function ----
    #---------------------------------------------------------------------- -

  preds <- purrr::map_df(include, ~{
    model <- .
    pred(all.mods[[model]], model)
  })
  pred(all.mods$`exp - Joint`, "exp")
    overlay <- preds %>%
      ggplot2::ggplot(ggplot2::aes(x = time, y = surv, colour = factor(trt))) +
      ggplot2::geom_line(size = 1) +
      ggplot2::geom_step(data = p$plot$data,ggplot2::aes(x = time, y = surv, group = trt, colour = "KM"), inherit.aes = FALSE, size = 0.5) +
      ggplot2::facet_wrap(~ model) +
      ggplot2::theme_minimal(base_size = 18) +
      ggplot2::scale_color_manual(values = c(cols, "black"), labels = c("Control", "Treatment", "KM"),
                                  name = "Treatment") +
      ggplot2::labs(y = "Survival", x = "Time",
                    title = glue::glue("Parametric Fits vs KM"))



    loglog_gomp <- p$plot$data %>%
      dplyr::filter(time != 0) %>%
      dplyr::mutate(surv = log(-log(surv))) %>%
      ggplot2::ggplot(ggplot2::aes(x = time, y = surv, colour = strata)) +
      ggplot2::geom_step(size = 1) +
      ggplot2::scale_color_manual(values = cols, name = "Treatment",
                                  labels = c("Control", "Treatment")) +
      ggplot2::theme_minimal(base_size = 18) +
      ggplot2::labs(y = "Log Cumulative Hazard",
                    x = "Time (linear scale)",
                    title = "Log Cumulative Hazard  vs Time for assessment of Gompertz")


    loglog_weib <- loglog_gomp +
      ggplot2::scale_x_continuous(trans = "log", labels = function(x) round(x,2)) +
      ggplot2::labs(x = "Time (log scale)",
                    y = "Log Cumulative Hazard",
                    title = "Log-log Plot for assessment of Weibull")


    quant_dat <- quantile(m, probs = seq(from = 0, to = 1, by = 0.01), conf.int = FALSE)  %>%
      t() %>%
      as.data.frame() %>%
      tidyr::drop_na() %>%
      purrr::set_names("control", "trt") %>%
      tibble::rownames_to_column("quantile")

    m1 <- lm(trt ~ control, data = quant_dat)

    aft_plot <- quant_dat %>%
      ggplot2::ggplot(ggplot2::aes(x = control, y = trt)) +
      ggplot2::geom_point(size = 2, colour = cols[[1]]) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, colour =cols[[2]]) +
      ggplot2::labs(caption = glue::glue("Estimated accelaration factor (line slope) = {round(coef(m1)[[2]], 2)}"),
                    y = "Percentile for Control",
                    x = "Percentile for Treatment",
                    title = "Q-Q Plot for assessment of constant acceleration factor") +
      ggplot2::theme_minimal(base_size = 16)


    dplyr::lst(km = p,
               loglog_gomp,
               loglog_weib,
               aft_plot,
               overlay)
}


#' Function name
#'
#' What does the function do in one sentence
#'
#' More detailed explanation of assumptions etc
#'
#' @param  First argument (repeat for more)
#'
#' @return  What does the function return?
#'
#' @examples
#'  Put some examples of how it works or just delete this
#'
#' @export
.smooth.haz<- function(dat,fits, nbreaks = 10, include){
  cols <- c("#002F6C", "#ED8B00")
  tlabs <- c("Control", "Treatment")
  m <- survival::survfit(survival::Surv(futime, fustat) ~ trt, data = dat)

  p <- survminer::ggsurvplot(m, conf.int = TRUE, risk.table = TRUE,
                             ggtheme = ggplot2::theme_minimal(base_size = 18),
                             palette = cols,
                             surv.median.line = "v",
                             pval = TRUE,
                             data = dat)

  risk_test <- summary(m)

  mt <- max(risk_test$time)


  breaks <- seq(from = 0.01, to = mt, length.out = nbreaks)

  km_sum <- summary(m, times = breaks, extend = TRUE)


  risk_table <- lapply(c("time", "n.risk", "strata") , function(x) km_sum[x]) %>%
    do.call(data.frame, .) %>%
    dplyr::mutate(strata = ifelse(strata == "trt=0", tlabs[[1]], tlabs[[2]]),
                  strata = factor(strata, levels = rev(tlabs)))






  muhaz_bws <- muhaz_fits(data = dat,
                          time_var = "futime",
                          event_var = "fustat",
                          nsamp = 100)

  muhaz_pws <- .muhaz_pw(data = dat,
                         time_var = "futime",
                         event_var = "fustat")

  muhaz_pw_plot <- .muhaz.pehaz.plot(muhaz = muhaz_bws,
                                     pehaz = muhaz_pws,
                                     out = "Survival",
                                     survfit = m,
                                     risk_tab = TRUE,
                                     breaks = 4,
                                     textsize = 8,
                                     basesize = 18,
                                     xlab = "Time")


  muhaz_ol_p <- .hazard.overlay(
    mod_type = "Weibull",
    fits = fits[include],
    muhaz = muhaz_bws,
    out = "Survival",
    subt = "",
    survfit = m,
    risk_tab = TRUE,
    breaks = 10,
    textsize = 8,
    basesize = 18,
    xlab = "Time"
  )



  dplyr::lst(
    muhaz_ol_p,
    muhaz_pw_plot,

  )

}




