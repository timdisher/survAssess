#' hazplot_dat
#'
#' @param trt muhaz output for treatment
#' @param control muhaz output for control
#'
#' @return A dataframe of smoothed hazards for plotting
hazplot_dat <- function(trt, control, settings, t_lab, c_lab){

  plot <- out %>% ggplot2::ggplot(ggplot2::aes(x = time, y = haz, colour = Treatment)) +
    ggplot2::geom_line() +
    ggplot2::coord_cartesian(xlim = c(0,75),
                             ylim = c(0, 0.05)) +
    ggplot2::scale_colour_es(values = c("#002F6C", "#ED8B00")) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous(breaks = seq(from = 0, to = 75, by = 5)) +
    ggplot2::labs(y = "Hazard",
                  x = "Time (months)",
                  title = glue::glue("Smoothed hazard - Settings: {settings}"))

  list(data = out, plot = plot)
}






#' muhaz_fits
#'
#' Created smoothed hazard plots and save them into a directory of interest.
#'
#' These plots were created to support HTA submissions that wanted to use smoothed
#' hazards as a guide for when to start parametric modeling of models that use
#' a mix of KM estimate and parametric extrapolation.
#'
#' @param data the ipd/pseudo ipd data you would like to use.
#' @param dir The output directory
#' @param name The name of the saved excel file for example OUTCOME-DBL_smoothed-hazard-plots
muhaz_fits <- function(data,
                       event_var = "out_os_event",
                       time_var = "out_os_time",
                       nsamp = 1) {
  cols <- c("#002F6C", "#ED8B00")
  if(! "trt" %in% colnames(data)) stop("Code requires a treatment variable named 'trt'")
  if(any(!data$trt %in% c(0,1))) stop("trt must be a numeric variable with values 0/1 only")
  if(! is.numeric(data[[time_var]])) stop("timevar must be numeric")

  if(!is.numeric(data[[event_var]])){
    data[[event_var]] <- as.numeric(data[[event_var]]) - 1

  }

  if(any(! data[[event_var]] %in% c(0,1))){
    stop("Muhaz requires numeric conversion of event variable to result in 0/1 only. Check data")
  }

  data <- data %>%
    dplyr::rename(event_var := !! event_var,
                  time_var := !! time_var)


  data_trt <- data %>% dplyr::filter(trt == "1")
  data_pbo <- data %>% dplyr::filter(trt == "0")


  ests <- purrr::map(1:nsamp, ~ {

    t_samp <- dplyr::sample_frac(data_trt, size = 1, replace = TRUE)
    c_samp <- dplyr::sample_frac(data_pbo, size = 1, replace = TRUE)

    ts <- survival::survfit(survival::Surv(time_var, event_var) ~ 1, data = t_samp)
    cs <- survival::survfit(survival::Surv(time_var, event_var) ~ 1, data = c_samp)



    haztreat = with(t_samp, muhaz::muhaz(time_var, event_var,min.time = 0, max.time = ts$time[max(cumsum(ts$n.risk >= 10))]))
    hazcontrol = with(c_samp, muhaz::muhaz(time_var, event_var,min.time = 0, max.time = cs$time[max(cumsum(cs$n.risk >= 10))]))

    t_df <-data.frame(haz = haztreat$haz.est,time = haztreat$est.grid, Treatment = 1)
    c_df <- data.frame(haz = hazcontrol$haz.est,time = hazcontrol$est.grid, Treatment = 0)



    rbind(t_df, c_df)


  }) %>%
    dplyr::bind_rows(.id = "iteration") %>%
    dplyr::mutate(iteration = as.numeric(iteration))

  return(ests)
}


#' muhaz_pw
#'
#' Fit peacewise hazards estimate with muhaz
#'
#'
#' @param  First argument (repeat for more)
#'
#' @return  What does the function return?
#'
#' @examples
#'  Put some examples of how it works or just delete this
#'
#' @export
.muhaz_pw <-function(data,
                     event_var = "out_os_event",
                     time_var = "out_os_time",
                     nsamp = 1) {

  if(! "trt" %in% colnames(data)) stop("Code requires a treatment variable named 'trt'")
  if(any(!data$trt %in% c(0,1))) stop("trt must be a numeric variable with values 0/1 only")
  if(! is.numeric(data[[time_var]])) stop("timevar must be numeric")

  if(!is.numeric(data[[event_var]])){
    data[[event_var]] <- as.numeric(data[[event_var]]) - 1

  }

  if(any(! data[[event_var]] %in% c(0,1))){
    stop("Muhaz requires numeric conversion of event variable to result in 0/1 only. Check data")
  }

  data <- data %>%
    dplyr::rename(event_var := !! event_var,
                  time_var := !! time_var)


  data_trt <- data %>% dplyr::filter(trt == "1")
  data_pbo <- data %>% dplyr::filter(trt == "0")

  ts <- survival::survfit(survival::Surv(time_var, event_var) ~ 1, data = data_trt)
  cs <- survival::survfit(survival::Surv(time_var, event_var) ~ 1, data = data_pbo)


  haztreat = with(data_trt, muhaz::pehaz(time_var, event_var, max.time = ts$time[max(cumsum(ts$n.risk >= 10))]))
  hazcontrol = with(data_pbo, muhaz::pehaz(time_var, event_var, max.time = cs$time[max(cumsum(cs$n.risk >= 10))]))

  t_df <-data.frame(haz = c(haztreat$Hazard, haztreat$Hazard[length(haztreat$Hazard)]),time = c(haztreat$Cuts), Treatment = 1)
  c_df <- data.frame(haz = c(hazcontrol$Hazard, hazcontrol$Hazard[length(haztreat$Hazard)]),time = hazcontrol$Cuts, Treatment = 0)



  ests <- rbind(t_df, c_df)

  return(ests)
}


#' .hazard.overlay
#'
#' Overlay hazard from flexsurvreg model and muhaz smoothed non-parametric hazard
#'
#'
#' @param  mod_type A string, the type of flexsurvreg model. Used to label the legend.
#' @param model the flexsurvreg model object
#' @param muhaz output from muhaz_fits stripped down to data output object
#' @param out A string. Outcome being modeled. Used in title.
#' @param subt A string. Subtitle (usually trial name).
#' @param survfit Survfit object used to make risk table
#' @param risk_tab logical indicating whether to add numbers at risk
#' @param breaks number of x axis breaks
#' @param textsize size of risk table values
#' @param basesize base size for ggplot theme

#'
#' @return  A plot of modeled vs non-parametrically smoothed hazards
#'
#' @export
.hazard.overlay <- function(mod_type,
                            fits,
                            muhaz,
                            out,
                            subt,
                            survfit = surv_os,
                            risk_tab = FALSE,
                            breaks = 4,
                            textsize = 8,
                            basesize = 18,
                            xlab = "Time (months)"){


  if(any(! c("Treatment", "haz", "time", "iteration") %in% colnames(muhaz))){

    stop("Column names for muhaz must be 'Treatment', 'haz', 'time', and 'iteration'")

  }


  mt <- ceiling(max(muhaz$time)) # Max time is always rounded up
  t_lab <- c("Control", "Treatment")
  breaks <- round(seq(from = 0, to = mt, length.out = 10))



  muhaz_p <- muhaz %>%
    dplyr::group_by(iteration) %>%
    dplyr::mutate(et = 1:dplyr::n()) %>%
    dplyr::group_by(et, Treatment) %>%
    dplyr::summarize(med = quantile(haz, probs = 0.5),
                     lwr = quantile(haz, probs = 0.025),
                     upr = quantile(haz, probs = 0.975),
                     time = mean(time)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-et) %>%
    dplyr::rename(haz = med) %>%
    dplyr::mutate(Treatment = ifelse(Treatment == 1, t_lab[[2]], t_lab[[1]]),
                  Treatment = factor(Treatment, levels = t_lab)#,
                  #type = "Flexible non-parametric (muhaz)"
                  )


 m.haz <- purrr::map_df(names(fits), ~{
   model <- fits[[.]]
   m.name <- .
   mod.haz <- predict(model, newdata = data.frame(trt = c(0, 1)), times = unique(muhaz$time) + 0.01,
                      type = "hazard") %>%
     dplyr::pull(.pred) %>%
     dplyr::bind_rows(.id = "trt") %>%
     dplyr::mutate(trt = ifelse(trt == 1, t_lab[[1]], t_lab[[2]]),
                   trt = factor(trt, levels = t_lab),
                   type = m.name,
                   lwr = NA,
                   upr = NA) %>%
     dplyr::rename(Treatment = trt, haz = .pred_hazard, time = .time)

   mod.haz
 })


  out_p <- m.haz %>%
    ggplot2::ggplot(ggplot2::aes(x = time, y = haz, colour = Treatment, linetype = type)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::facet_wrap(~type) +
    ggplot2::scale_colour_manual(values =  c("#002F6C", "#ED8B00")) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr,x = time, ymax = upr, y = haz, group = Treatment, fill = Treatment),alpha = 0.2, data = muhaz_p,
                         inherit.aes = FALSE,
                         fill = "black") +
    ggplot2::geom_line(ggplot2::aes(x = time, y = haz, group = Treatment, colour = "Flexible non-parametric hazard"),colour = "black",
                       alpha = 1, data = muhaz_p,
                         inherit.aes = FALSE) +
    ggplot2::theme_minimal(base_size = basesize) +
    ggplot2::labs(x = "Time", y = "Hazard", linetype = "Model Type",
                  title = glue::glue("Modeled vs Non-parametrically smoothed hazards {out}"),
                  subtitle = subt) +
    ggplot2::scale_x_continuous(breaks = breaks) +
    ggplot2::coord_cartesian(xlim = c(0, mt))


  return(out_p + ggplot2::labs(caption = "*Smoothed hazards restricted to timepoints with at least 10 patients at risk"))





}


#' .muhaz.pehaz.plot
#'
#' Assess whether smoothed hazards are reasonable representation of piecewise hazards
#'
#'
#' @param pehaz output from .muhaz_pw
#' @param muhaz output from muhaz_fits stripped down to data output object
#' @param out A string. Outcome being modeled. Used in title.
#' @param subt A string. Subtitle (usually trial name).
#' @param survfit surfvit object used to make risk table
#' @param risk_tab logical indicating whether to add numbers at risk
#' @param breaks number of x axis breaks
#' @param textsize size of risk table values
#' @param basesize base size for ggplot theme
#'
#' @return  A plot of modeled vs non-parametrically smoothed hazards
#'
#' @export
.muhaz.pehaz.plot <- function(muhaz,
                              pehaz,
                              out,
                              subt = "",
                              survfit = surv_os,
                              risk_tab = FALSE,
                              breaks = 4,
                              textsize = 8,
                              basesize = 18,
                              xlab = "Time (months)"){


  if(any(! c("Treatment", "haz", "time") %in% colnames(muhaz))){

    stop("Column names for muhaz must be 'Treatment', 'haz', 'time'")

  }


  mt <- ceiling(max(muhaz$time)) # Max time is always rounded up
  t_lab <- c("Control", "Treatment")
  breaks <- round(seq(from = 0, to = mt, length.out = 10))



  muhaz_p <- muhaz %>%
    dplyr::group_by(iteration) %>%
    dplyr::mutate(et = 1:dplyr::n()) %>%
    dplyr::group_by(et, Treatment) %>%
    dplyr::summarize(med = quantile(haz, probs = 0.5),
                     lwr = quantile(haz, probs = 0.025),
                     upr = quantile(haz, probs = 0.975),
                     time = mean(time)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-et) %>%
    dplyr::rename(haz = med) %>%
    dplyr::mutate(Treatment = ifelse(Treatment == 1, t_lab[[2]], t_lab[[1]]),
                  Treatment = factor(Treatment, levels = t_lab),
                  type = "Flexible non-parametric (muhaz)")

  pehaz_p <- pehaz %>%
    dplyr::mutate(Treatment = ifelse(Treatment == 1, t_lab[[2]], t_lab[[1]]),
                  Treatment = factor(Treatment, levels = t_lab),
                  type = "Piecewise hazard")


  out_p <- pehaz_p %>%
    ggplot2::ggplot(ggplot2::aes(x = time, y = haz, colour = Treatment, linetype = type)) +
    ggplot2::geom_step(size = 1) +
    ggplot2::scale_colour_manual(values = c("#002F6C", "#ED8B00")) +
    ggplot2::geom_line(data = muhaz_p) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr,x = time, ymax = upr, y = haz, group = Treatment, fill = Treatment),alpha = 0.2, data = muhaz_p) +
    ggplot2::scale_fill_manual(values = c("#002F6C", "#ED8B00")) +
    ggplot2::theme_minimal(base_size = basesize) +
    ggplot2::labs(x = "Time", y = "Hazard", linetype = "Model Type",
                  title = glue::glue("Piecewise vs Non-parametrically smoothed hazards {out}"),
                  subtitle = subt) +
    ggplot2::scale_x_continuous(breaks = breaks) +
    ggplot2::coord_cartesian(xlim = c(0, mt))




  if(risk_tab){

    km_sum <- summary(survfit, times = breaks, extend = TRUE)
    risk_table <- lapply(c("time", "n.risk", "strata") , function(x) km_sum[x]) %>%
      do.call(data.frame, .) %>%
      dplyr::mutate(strata = ifelse(strata == "trt=0", t_lab[[1]], t_lab[[2]]),
                    strata = factor(strata, levels = rev(t_lab)))



    out_p <- out_p + ggplot2::labs(x = ggplot2::element_blank()) # Remove xlab from upper plot

    rt_p_strata <- risk_table  %>%
      ggplot2::ggplot(ggplot2::aes(x = time, y = strata, label = n.risk)) +
      ggplot2::geom_text(size = textsize,ggplot2::aes(colour = strata), show.legend = FALSE) +
      ggplot2::scale_color_manual(values = c("#002F6C", "#ED8B00")) +
      ggplot2::theme_minimal(base_size = basesize) +
      ggplot2::coord_cartesian(xlim = c(0, 60)) +
      ggplot2::labs(y = "No. at Risk", x = xlab) +
      ggplot2::scale_x_continuous(breaks = breaks) +
      ggplot2::coord_cartesian(xlim = c(0, mt))


    comb <- cowplot::plot_grid(out_p, rt_p_strata, align = "v", ncol = 1, rel_heights = c(0.75, 0.25),
                               axis = "rl")
    return(cowplot::ggdraw(cowplot::add_sub(comb, "*Smoothed hazards restricted to timepoints with at least 10 patients at risk")) )



  }

  return(out_p + ggplot2::labs(caption = "*Smoothed hazards restricted to timepoints with at least 10 patients at risk"))





}
