inv_logit <-function (x) {
    exp(x)/(exp(x) + 1)
}
 

my_ilink <- function(x, link) {
  switch(link, identity = x, log = exp(x), logm1 = expp1(x), 
         log1p = expm1(x), inverse = 1/x, sqrt = x^2, `1/mu^2` = 1/sqrt(x), 
         tan_half = 2 * atan(x), logit = inv_logit(x), probit = pnorm(x), 
         cauchit = pcauchy(x), cloglog = inv_cloglog(x), probit_approx = pnorm(x), 
        softplus = log1p_exp(x), stop2("Link '", link, "' not supported."))
  }

inv_logit <- function(x) {
    x[x>500] <- 500
    exp(x)/(exp(x) + 1)
}

posterior_predict.inla <- function(object, newdata = NULL, ndraws=100, include_re = TRUE, new_random_levels = FALSE) {
    draws <- inla.posterior.sample(n=ndraws, result = object) %>%
        suppressWarnings() %>%
        suppressMessages()
    contents <- object$misc$configs$contents

  contrasts <- object$.args$contrasts
  form <- object$.args$formula
  gf <- INLA:::inla.interpret.formula(form)
  ## get the linear predictor
  if (is.null(newdata)) {
      newdata <- object$.args$data %>%
          filter(!is.na(!!sym(gf$response)))
  }
  
  ## wch <- which(contents$tag=='Predictor')
  ## i_data <- contents$start[wch]:(1+contents$length[wch])
  ## mm <- mm[i_data,]
  ## } else {
  ## Fixed effects
  Xmat <- model.matrix(update(gf$fixf, NULL ~ .), newdata, contrasts = contrasts)
  
  ## Random effects
  if (!is.null(gf$randf)) {
      nms <- unlist(lapply(gf$random.spec, `[[`, "label"))
      ## newdata <- newdata %>% mutate(across(!!nms, ~ paste0('NEW',.x) %>% as.factor()))  
      Zmat <- lapply(nms, function(n) model.matrix(as.formula(paste0("~ 0 + ",n)), newdata))
      Zmat <- do.call('cbind', Zmat)
  }
  ## }
  ## Fixed effects
  nms <- colnames(Xmat)
  i_lp <- contents$start[contents$tag %in% nms]
  lp <- t(sapply(draws, function(x) x$latent[i_lp]))
  if (nrow(lp)==1) lp <- t(lp)
  b <- tcrossprod(Xmat, lp)
  link <- object$misc$linkfunctions$names
  out <- my_ilink(b, link) 

  ## family theta
  if ((ntheta <- object$misc$configs$ntheta) > 0) {
    theta <- vector('list', ntheta)
    for (t in 1:ntheta) {
      theta[[t]] <- sapply(draws, function(x) 1/sqrt(x$hyperpar[[t]]))
    }
  }

  ## Random effects
  if (gf$n.random>0 & include_re) {
      nms <- unlist(lapply(gf$random.spec, `[[`, "label"))
      i_lp <- contents$start[contents$tag %in% nms]
      lp <- ii_lp <- vector('list', length(nms))
      for (i in 1:length(nms)) {
          ii_lp[[i]] <- i_lp[[i]]:(i_lp[[i]] + (contents$length[contents$tag %in% nms] - 1)[i])
          lp[[i]] <- t(sapply(draws, function(x) x$latent[ii_lp[[i]]]))
          if (new_random_levels) lp[[i]] <- apply(lp[[i]], 2,
                                                  function(x) rnorm(n = length(x), mean = 0, sd = theta[[i]]))
      }
      lp <- do.call(`cbind`, lp)
      ## i_lp <- unlist(lapply(1:length(i_lp),
      ##                       function(i)
      ##                           i_lp[i]:(i_lp[i] + (contents$length[contents$tag %in% nms] - 1)[i])))
      
      ## i_lp <- i_lp:(i_lp + contents$length[contents$tag %in% nms] - 1)
      ## lp <- t(sapply(draws, function(x) x$latent[i_lp]))
      ## if (new_random_levels) lp <- apply(lp, 2, function(x) rnorm(n = length(x), mean = 0, sd = theta[[2]]))
      r <- tcrossprod(Zmat, lp)
      link <- object$misc$linkfunctions$names
      out <- my_ilink(b+r, link) 
  }

  fam <- object$.args$family
  N <- nrow(out)
  if (fam == 'gaussian') {
    rdist <- function(N, x, sigma) rnorm(N, x, sigma)
    out <- apply(out, 2, rdist, N=N, sigma=theta[[1]])
  }
  if (fam=='poisson') {
    rdist <- function(N, x) rpois(N, x)
    out <- apply(out, 2, rdist, N=N)
  }
  if (fam=="nbinomial") {
      wch <- grep(
          "size for the nbinomial observations",
          names(draws[[1]]$hyperpar)
      )
      size <- sapply(draws, function(x) x$hyperpar[[wch]])
      ## size <- inla.hyperpar.sample(n = ndraws, result = mod.inla)[, wch]
      ## rdist <- function(N, x, size) MASS::rnegbin(N, mu = x, size) #rnbinom(N, mu = x, prob)
      rdist <- function(N, x, size) rnbinom(N, mu = x, size) #rnbinom(N, mu = x, prob)
      out <- apply(out, 2, rdist, N=N, size=size)
  }
  if (fam=="binomial") {
    ## size <- ifelse(is.null(object$.args$Ntrials), 1, object$.args$Ntrials)
      if (is.null(object$.args$Ntrials)) {
          size <- rep(1, N)
      } else {
          size <- object$.args$Ntrials[1:N]
      }
    rdist <- function(N, x, size) rbinom(N, prob = x, size)
    out <- apply(out, 2, rdist, N=N, size=size)
  }
  if (fam=="beta") {
      wch <- grep(
          "precision parameter for the beta observations",
          names(draws[[1]]$hyperpar)
      )
      phi <- sapply(draws, function(x) x$hyperpar[[wch]])
      rdist <- function(N, x, phi) rbeta(N, shape1 = x*phi, shape2 = (1 - x)*phi)
      out <- apply(out, 2, rdist, N=N, phi=phi)
  }
  if (fam=="betabinomial") {
      wch <- grep(
          "overdispersion for the betabinomial observations",
          names(draws[[1]]$hyperpar)
      )
      phi <-1/sapply(draws, function(x) x$hyperpar[[wch]])
      if (is.null(object$.args$Ntrials)) {
          size <- rep(1, N)
      } else {
          size <- object$.args$Ntrials[1:N]
      }
      rdist <- function(N, x, phi, size) {
          rbinom(N,
                 prob = rbeta(N, shape1 = x*phi, shape2 = (1 - x)*phi),
                 size = size)
          }
      out <- sapply(seq_len(ncol(out)), function(i)
          rdist(N = N, x = out[,i], phi = phi[i], size = size)) %>%
          matrix(ncol = ncol(out), byrow = TRUE)
      ## a <- mapply(FUN = rdist, out, phi, MoreArgs = list(N = N, size = size))
  }
  if (fam=="zeroinflatedpoisson1") {
      wch <- grep(
          "zero-probability parameter for zero-inflated poisson_1",
          names(draws[[1]]$hyperpar)
      )
      phi <- sapply(draws, function(x) x$hyperpar[[wch]])
      rdist <- function(N, x, phi) rbinom(N, size = 1, prob = 1 - phi) * rpois(N, lambda = x)
      out <- apply(out, 2, rdist, N=N, phi=phi)
  }
  
  
  t(out)
}





pit_hist_plot <- function(mod, i.mod) {
    pit <- mod$cpo$pit[i.mod]
    g <- ggplot(data = NULL, aes(x = pit)) +
        ## geom_histogram() +
        geom_density() +
        scale_y_continuous("Density") +
        scale_x_continuous("PIT") +
        theme_bw()
    return(g)
}
pit_resid_plot <- function(mod, i.mod) {
    pit <- mod$cpo$pit[i.mod]
    g <- ggplot(data = NULL, aes(y = pit, x = mod$summary.fitted.values$mean[i.mod])) +
        geom_point() +
        scale_y_continuous("PIT values") +
        scale_x_continuous("Posterior fitted values") +
        theme_bw()
    return(g)
}
pit_plot <- function(mod, i.mod) {
    pit <- na.omit(mod$cpo$pit[i.mod])
    data <- na.omit(mod$.args$data[i.mod,])
    n <- nrow(data)
    g <- ggplot(data = NULL, aes(y = pit, x = 1:n)) +
        geom_point() +
        scale_y_continuous("PIT values") +
        scale_x_continuous("1:n") +
        theme_bw()
    return(g)
}
pit_qq_plot <- function(mod, i.mod, logit_scale = TRUE) {
    logit <- function(x) log(x/(1-x))
    pit <- na.omit(mod$cpo$pit[i.mod])
    data <- na.omit(mod$.args$data[i.mod,])
    n <- nrow(data)
    uniquant <- (1:n)/(n+1)
    g <- ggplot(data = NULL, aes(y = sort(pit), x = uniquant))
    if (logit_scale) 
        g <- ggplot(data = NULL, aes(y = logit(sort(pit)), x = logit(uniquant)))
    g <- g + geom_abline(slope = 1) +
        geom_point() +
        scale_y_continuous("Sorted PIT values") +
        scale_x_continuous("Uniform quantiles") +
        theme_bw()
    return(g)
}
pp_check.foo <- function(object, type = c("multiple", "overlaid"), ...) {
  type <- match.arg(type)
  y <- object[["y"]]
  yrep <- object[["yrep"]]

  if (type == "overlaid") {
    ppc_dens_overlay(y, yrep, ...) 
  } else {
    ppc_hist(y, yrep, ...)
  }
}
