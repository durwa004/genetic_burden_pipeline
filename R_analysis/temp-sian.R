library(dplyr)
library(forcats)
set.seed(50)
d <- expand.grid(Season=c("Fall", "Winter", "Spring", "Summer"),
                 Sex=c("M","F"), rep=1:3) %>% 
  select(-rep) %>%
  mutate(Season = fct_inorder(Season), Sex=fct_inorder(Sex)) %>%
  mutate(Age=round(runif(nrow(.), 3,12),1)) %>%
  mutate(y = 3*as.numeric(Season) + 2*as.numeric(Sex) + 1*Age + rnorm(nrow(.)))

library(emmeans)
library(dplyr)
library(forcats)
library(ggplot2)
modelemms <- function(model, vs, dat) {
  stopifnot(sapply(dat[vs], "class") %in% c("numeric", "factor"))
  lapply(vs, function(v) {
    if(class(dat[[v]])=="factor") { atq <- levels(dat[[v]])
    } else { atq <- quantile(dat[[v]], c(0.25, 0.75)) }
    atq <- atq %>% list() %>% setNames(v)
    model %>% emmeans(v, at=atq, data=dat) %>%
      summary() %>% mutate(variable=v) %>% rename(value=!!v) %>% 
      mutate(value=paste(value))
  }) %>% bind_rows() %>% mutate(variable=factor(variable, levels=vs))
}

vars <- c("Season", "Sex", "Age")
f <- paste("y", paste(vars, collapse="+"), sep="~")
modelx <- do.call("lm", list(formula=as.formula(f), data=as.name("d")))
emms <- modelemms(modelx, vars, d) %>% 
  mutate(valueX = paste0(variable, ": ", value) %>% fct_inorder,
         valueX = factor(valueX, levels=rev(levels(valueX))))
ggplot(emms) + aes(valueX, emmean, ymin=lower.CL, ymax=upper.CL) + 
  geom_pointrange() + facet_grid(variable~., scales="free", space="free") + 
  coord_flip() + xlab("") + ylab("Estimated Marginal Mean of Y")
