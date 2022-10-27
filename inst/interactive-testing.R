dat <- .simsurv(n = 150, maturity = "max")

tests <- .test.ph.aft(dat$dat)

tests$fits

plots <- .km.plots(fits = tests$all.mods,
                   dat = dat$dat,
                   include = names(tests$all.mods))

plots$loglog_gomp
plots$loglog_weib
plots$aft_plot
plots$overlay
plots$km

sm.haz.p <- .smooth.haz(
  dat = dat$dat,
  nbreaks = 10,
  fits = tests$all.mods,
  include = names(tests$all.mods)
)

sm.haz.p$muhaz_ol_p
sm.haz.p$muhaz_pw_plot

dat$dist
dat$ph.aft
