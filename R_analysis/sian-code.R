library(emmeans)
library(dplyr)
library(forcats)
library(ggplot2)

modelemms2 <- function(model, dat, main=TRUE, type="response", ...) {
  rr <- list(...)
  vs <- terms(model) %>% attr("factors") %>% colnames()
  vs <- gsub(":", "|", vs, fixed=TRUE)
  xs <- vs[!grepl("|", vs, fixed=TRUE)]
  vsi <- vs[grepl("|", vs, fixed=TRUE)] %>% strsplit("|", fixed=TRUE) %>% 
    unlist() %>% unique()
  if(!main) { vs <- vs[!vs %in% vsi] }
  atq <- lapply(xs, function(v) {
    if(class(dat[[v]])=="factor") { levels(dat[[v]])
    } else { quantile(dat[[v]], c(0.25, 0.75),na.rm=TRUE) } %>% unname()
  }) %>% setNames(xs)
  out <- lapply(vs, function(v) {
    model %>% emmeans(as.formula(paste0("~", v)), at=atq, data=dat) %>% 
      summary(type=type) %>% as.data.frame() %>% mutate(term=v)
  }) %>% bind_rows()
  names(out)[names(out)=="response"] <- "emmean"
  out <- out %>% mutate(term=gsub("|", ":", term, fixed=TRUE) %>% fct_inorder())
  xx <- c("emmean", "SE", "df", "lower.CL", "upper.CL")
  terms <- names(out)[!names(out) %in% c("term", xx)]
  out <- out[c("term", terms, xx)]
  out
}

## round to specified digits, maintaining NA values
format_digits <- function(x, digits=2) {
  k <- !is.na(x)
  x[k] <- formatC(x[k], digits=digits, format="f")
  factor(x, levels=unique(x[k]))
}

## create a summary variable across all terms
## the pattern is how to create the text for each variable; 
##      the first %s will be replaced with the variable, the second %s with the value
## the sep is how to separate two variables
add_value <- function(x, pattern="%s: %s", sep=", ") {
  xx <- c("term", "emmean", "SE", "df", "lower.CL", "upper.CL")
  terms <- colnames(x)[!colnames(x) %in% xx]
  foo <- lapply(terms, function(y) ifelse(is.na(x[[y]]), NA, sprintf(pattern, y, x[[y]])))
  foo <- do.call(cbind, foo)  
  foo <- apply(foo, 1, function(y) paste(y[!is.na(y)], collapse=sep))
  factor(foo, levels=rev(foo))
}

plotemms <- function(emms, color=FALSE, shape=TRUE, digits=2, pattern="%s (%s)", sep=" * ", ...) {
  rr <- list(...)
  terms <- names(emms)[!names(emms) %in% c("term", "emmean", "SE", "df", "lower.CL", "upper.CL")]
  terms.numeric <- terms[sapply(emms[terms], class)=="numeric"]
  emms$X <- NA
  if(identical(FALSE, color)) { color <- "X" }
  else if(isTRUE(color)) if (length(terms.numeric)==1) { color <- terms.numeric } else { color <- NA }
  if(identical(FALSE, shape)) { shape <- NULL }
  else if(isTRUE(shape)) if (length(terms.numeric)==1) { shape <- terms.numeric } else { shape <- NULL }
  for(v in terms.numeric) {
    k <- match(v, names(rr))
    if(!is.na(k)) { emms[[v]] <- format_digits(emms[[v]], digits=rr[[k]]) }
    else { emms[[v]] <- format_digits(emms[[v]], digits=digits) }
  }
  emms$value <- add_value(emms, pattern=pattern, sep=sep)
  colors <- scales::hue_pal(h=c(0,360) + 15, c=100, l=65, h.start=0, direction=1)(2)
  ggplot(emms) + 
    aes(x = value, y = emmean, ymin=lower.CL, ymax=upper.CL)  + aes_string(color=color, shape=shape, fill=color) + 
    geom_pointrange() + facet_grid(term~., scales="free", space="free") + 
    theme(strip.text.y = element_text(size=8), axis.text.y = element_text(size=2), 
          axis.text.x = element_text(size=8), axis.title.x = element_text(size=10)) +
    coord_flip() + xlab("") + 
    scale_shape_manual(values=c(2,24), na.value=21) +
    discrete_scale("colour", "manual", function(n) colors, na.value="black") +
    discrete_scale("fill", "manual", function(n) colors, na.value="black") +
    theme(legend.position="none")
}

getsig <- function(emms, model) {
  sigs <- model %>% car::Anova() %>% broom::tidy() %>% filter(p.value<0.05) %>% select(term) %>% unlist()
  emms %>% filter(term %in% sigs)
}