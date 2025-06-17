
CI__clean_inla_formula <- function(form, data) {
    term.labels <- attr(terms(form), 'term.labels')
    vars <- all.vars(form)[-1]
    vars <- vars[-grep('spatial|spde', vars)]
    for (v in vars) {
        if ((data %>%
             pull(!!sym(v)) %>%
             levels() %>%
             length() < 2) &
            (data %>%
             pull(!!sym(v)) %>% class() != "numeric")) {
            ## fixed effects
            if (v %in% term.labels) {
                wch <- which(v == term.labels)
                if (length(wch) ==0) next
                term.labels <- term.labels[-wch]
                form <- reformulate(term.labels, response = all.vars(form)[1], intercept = FALSE)
            }
            ## random effects
            wch <- grep(paste0('^f.',v), term.labels, perl = TRUE)
            if (length(wch)==0) next
            term.labels <- term.labels[-wch]
            form <- reformulate(term.labels, response = all.vars(form)[1], intercept = FALSE)
        }
    }
    form<-update(form, .~.+Intercept)
    form
}

fitModel <- function(form, stack.est, spde) {
  environment(form)<-environment()  
  inla(form,
         data = inla.stack.data(stack.est),
         family= 'binomial',
         Ntrials=Total,
         control.predictor = list(compute = TRUE,
                                  link = 1,
                                  A = inla.stack.A(stack.est)
                                  ),
         control.compute = list(config = TRUE),
         verbose = FALSE)
}
