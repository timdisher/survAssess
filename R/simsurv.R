#' .simsurv
#'
#' Simulate survival data
#'
#' @param n Number of patients in each trial.
#' @param lambda value of lambda (scale) for PH Weibull. Must be > 0.
#' @param gamma value of gamma (shape) for PH weibull. Must be > 0 and gives exponential if 1 and dist is weibull.
#' @param maxt Time at which patients are administratively censored
#' @param beta a named vector with the log hazard ratio for pertuzumab pre-treatment
#' @param ntrial number of trials to simulate
#' @param dist either weibull or gompertz
#'
#'
#' @return  What does the function return?
#'
#' @examples
#'
#' .sim.surv(n = 1000,
#'           lambda = 0.2,
#'           gamma = 1,
#'           maxt = 5,
#'           beta = c(pert = log(1.5)),
#'           ntrial = 25)
#'
#' @export
.simsurv <- function(n){

############################################################################## #
############################################################################## #
#
# 1. Set up simulated survival based on random criteria ----
#
#     Section Notes
#
############################################################################## #
############################################################################## #


 ph.aft <- sample(c(TRUE), 1)
 dist <- sample(c("exp", "weibullPH", "gompertz"), 1)
 dataset <- sample(c("bc"), 1)

 if(dataset == "bc"){
   dat <- flexsurv::bc %>%
     dplyr::mutate(trt = ifelse(group == "Good", 1, 0)) %>%
     dplyr::rename(
       futime =recyrs,
       fustat = censrec
     )

   maxt = max(dat$futime)
 }

 trt <- rbinom(n, 1, 0.5)
 x <- data.frame(trt)



 if(ph.aft | dist == "exp"){

  m <- flexsurv::flexsurvreg(survival::Surv(futime, fustat) ~ trt, dist = dist, data = dat)

 } else{
  m.c  <- flexsurv::flexsurvreg(survival::Surv(futime, fustat) ~ 1, dist = dist, data = dat %>% dplyr::filter(trt == 0))
  m.t  <- flexsurv::flexsurvreg(survival::Surv(futime, fustat) ~ 1, dist = dist, data = dat %>% dplyr::filter(trt == 1))
 }

 if(dist == "exp"){
 sim.dist <- "weibull"
 gamma <- 1
 lambda <- exp(coef(m)["rate"])
 }

 if(dist == "weibullPH"){
 sim.dist <- "weibull"
  if(ph.aft){

  lambda <- exp(coef(m)[["scale"]])
  gamma <- exp(coef(m)[["shape"]])

  } else{

  lambda.c <- exp(coef(m.c)[["scale"]])
  gamma.c <- exp(coef(m.c)[["shape"]])

  lambda.t <- exp(coef(m.t)[["scale"]])
  gamma.t <- exp(coef(m.t)[["shape"]])


  }

 }




 if(dist == "gompertz"){
   sim.dist <- "gompertz"
   if(ph.aft){

     lambda <- exp(coef(m)[["rate"]])
     gamma <- coef(m)[["shape"]]

   } else{

     lambda.c <- exp(coef(m.c)[["rate"]])
     gamma.c <- coef(m.c)[["shape"]]

     lambda.t <- exp(coef(m.t)[["rate"]])
     gamma.t <- coef(m.t)[["shape"]]


   }

 }


 ############################################################################## #
 ############################################################################## #
 #
 # 2. Simulations done with simsurv ----
 #
 #     Section Notes
 #
 ############################################################################## #
 ############################################################################## #



 if(!dist %in% c("weibull"))

if(ph.aft | dist == "exp"){
  surv <- simsurv::simsurv(
    dist = sim.dist,
    lambdas = lambda,
    gammas = gamma, # Exponential if 1
    x = x,
    maxt = maxt,
    betas = c(trt = coef(m)[["trt"]])
  ) %>%
    dplyr::rename(
      futime = eventtime,
      fustat = status
    )

  dat.out <- cbind(surv, x)

} else {




  cntl <- simsurv::simsurv(
    dist = sim.dist,
    lambdas = lambda.c,
    gammas = gamma.c, # Exponential if 1
    x = x %>% dplyr::filter(trt == 0),
    maxt = maxt
  ) %>%
    dplyr::rename(
      futime = eventtime,
      fustat = status
    )

  trt <- simsurv::simsurv(
    dist = sim.dist,
    lambdas = lambda.t,
    gammas = gamma.t, # Exponential if 1
    x = x %>% dplyr::filter(trt == 1),
    maxt = maxt
  ) %>%
    dplyr::rename(
      futime = eventtime,
      fustat = status
    )

  dat.out <- rbind(
    cbind(cntl, x %>% dplyr::filter(trt == 0)),
    cbind(trt, x %>% dplyr::filter(trt == 1))
  )
}

 ############################################################################## #
 ############################################################################## #
 #
 # 3. Simulations done without simsurv----
 #
 #     Section Notes
 #
 # Required for anything other than PH Weibull, gompertz, and exponential
 #
 ############################################################################## #
 ############################################################################## #



 if(dist %in% c("weibull")){


 if(dist == "weibull"){
   sim.dist <- "weibull"
   if(ph.aft){


     gamma <- exp(coef(m)[["shape"]])
     lambda <- exp(coef(m)[["scale"]])


     dat.out <- tibble::tibble(
       futime = rweibull(n,
                         shape = exp(coef(m)[["shape"]]),
                         scale = exp(coef(m)[["scale"]] + coef(m)[["trt"]]*x$trt)),
       fustat = 1,
       trt = x$trt) %>%
       dplyr::mutate(cen = ifelse(futime > maxt, 1, 0),
                     futime2 = ifelse(cen, maxt, futime),
                     fustat2 = ifelse(cen, 0, 1))


   } else{


     cntl <- tibble::tibble(
       trt = 0,
       futime = rweibull(n,
                         shape = exp(coef(m.c)[["shape"]]),
                         scale = exp(coef(m.c)[["scale"]]))) %>%
       dplyr::mutate(cen = ifelse(futime > maxt, 1, 0),
                     futime = ifelse(cen, maxt, futime),
                     fustat = ifelse(cen, 0, 1))


     trt <- tibble::tibble(
       trt = 1,
       futime = rweibull(n,
                         shape = exp(coef(m.t)[["shape"]]),
                         scale = exp(coef(m.t)[["scale"]]))) %>%
       dplyr::mutate(cen = ifelse(futime > maxt, 1, 0),
                     futime = ifelse(cen, maxt, futime),
                     fustat = ifelse(cen, 0, 1))


     dat.out <- rbind(cntl, trt)

     }

 }

 }


 dplyr::lst(
  dat = dat.out,
  dist,
  ph.aft,
  m, m.c, m.t
  )


  }


