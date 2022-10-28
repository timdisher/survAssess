dat <- .simsurv(n = 100, maturity = "max")


tests <- .test.ph.aft(dat$dat)


plots <- .km.plots(fits = tests$all.mods,
                   dat = dat$dat,
                   include = names(tests$all.mods))

sm.haz.p <- .smooth.haz(
  dat = dat$dat,
  nbreaks = 10,
  fits = tests$all.mods,
  include = names(tests$all.mods)
)


plots$km
plots$loglog_gomp
plots$loglog_weib
plots$aft_plot
plots$overlay
tests$fits


sm.haz.p$muhaz_pw_plot
sm.haz.p$muhaz_ol_p

dat$dist
dat$ph.aft
