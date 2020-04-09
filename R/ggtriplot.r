# 
#  ggtriplot.r
#  
#  Copyright 2011 Vincent Q. Vu.
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# 
#'  Triplot for RDA using ggplot2
#'
#' @param pcobj               an object returned by prcomp() or princomp() or rda() (from the vegan package)
#' @param choices             which PCs to plot
#' @param scale               covariance biplot (scale = 1), form biplot (scale = 0). When scale = 1, the inner product between the variables approximates the covariance and the distance between the points approximates the Mahalanobis distance.
#' @param pc.biplot           for compatibility with biplot.princomp()
#' 
#' @param obs.scale           scale factor to apply to observations
#' @param var.scale           scale factor to apply to variables
#' 
#' @param groups              optional factor variable indicating the groups that the observations belong to. If provided the points will be colored according to groups
#' @param ellipse             draw a stat_ellipse() for each group?
#' @param ellipse.level       the confidence level at which to draw an ellipse (default is 0.95), or, if type="euclid", the radius of the circle to be drawn.
#' 
#' @param ellipse.type        what type of stat_ellipse to draw : 't', 'norm' or 'euclid'
#' @param ellipse.size        width of the path of the ellipse drawn (size aes of stat_ellipse)
#' @param ellipse.size.range  if size set to 'groups' or 'gr', will use this as size range for the ellipse line size
#' @param ellipse.show.legend logical. Should this layer be included in the legends? NA, the default, includes if any aesthetics are mapped. FALSE never includes, and TRUE always includes.
#' @param ellipse.alpha       alpha transparency value for the fill aes of the ellipse
#' 
#' @param labels              optional vector of labels for the observations
#' @param labels.size         size of the text used for the labels
#' @param alpha               alpha transparency value for the points (0 = transparent, 1 = opaque)
#' 
#' @param circle              draw a correlation circle? (only applies when prcomp was called with scale = TRUE and when var.scale = 1)
#' @param axes.lang           in which language should the axis text be displayed : 'FR' for french, everything else is english for now
#' @param name.adjust         adjustment factor the placement of the names, >= 1 means farther from the arrow
#' 
#' @param var.arrows          logical; draw arrows for the variables?
#' @param varname.size        size of the text for variable names
#' @param varname.abbrev      whether or not to abbreviate the variable names
#' 
#' @param var.arrow.scaling   factor to which the length of the variables arrows will be scaled
#' @param var.arrow.color     color aes of the variables arrows
#' @param var.arrow.alpha     alpha transparency value for the variables arrows
#' 
#' @param bi.arrows           logical; draw arrows for the biplot?
#' @param biname.size         size of the text for biplot names
#' @param biname.abbrev       whether or not to abbreviate the biplot names
#' 
#' @param bi.arrow.scaling    factor to which the length of the biplot arrows will be scaled
#' @param bi.arrow.color      color aes of the biplot arrows
#' @param bi.arrow.alpha      alpha transparency value for the biplot arrows
#' 
#' @param obs.size            size aes of observation
#' param obs.size.binned     weither to use scale_size_binned() or not
#' @param obs.color           color aes of observation; default 'black' help distinguishing the points imo
#'
#' @return                    a ggplot2 plot
#' @export
#' @examples
#'   data(wine)
#'   wine.pca <- prcomp(wine, scale. = TRUE)
#'   print(ggbiplot(wine.pca, obs.scale = 1, var.scale = 1, groups = wine.class, ellipse = TRUE, circle = TRUE))
#'
ggtriplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE,
                     obs.scale = 1 - scale, var.scale = scale, 
                     groups = NULL, ellipse = FALSE, ellipse.level = 0.69,
                     ellipse.type = 'norm', ellipse.size = 0.5, ellipse.size.range = c(0.05,2),
                     ellipse.show.legend = NA, ellipse.alpha = 0.15,
                     labels = NULL, labels.size = 3, alpha = 1, 
                     circle = FALSE, circle.prob = 0.69, axes.lang = 'EN', name.adjust = 1.5,
                     var.arrows = TRUE, varname.size = 3, varname.abbrev = FALSE, arrow.alpha = 1,
                     var.arrow.scaling = 1, var.arrow.color = muted('red'),
                     bi.arrows = TRUE, biname.size = 3, biname.abbrev = F,
                     bi.arrow.scaling = 1, bi.arrow.color = muted('blue'), #bi.arrow.alpha = 1,
                     obs.size = 2, obs.color = 'black', obs.fill = NULL, obs.pal = NULL, #obs.binned = F, 
                     ...)
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  library(ggnewscale)
  
  stopifnot(length(choices) == 2)
  
  # Recover the SVD
  #if(inherits(pcobj, 'prcomp')){
  #  nobs.factor <- sqrt(nrow(pcobj$x) - 1)
  #  d <- pcobj$sdev
  #  u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
  #  v <- pcobj$rotation
  #  exp_var <- d[choices]^2/sum(d^2)
  #} else if(inherits(pcobj, 'princomp')) {
  #  nobs.factor <- sqrt(pcobj$n.obs)
  #  d <- pcobj$sdev
  #  u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
  #  v <- pcobj$loadings
  #  exp_var <- d[choices]^2/sum(d^2)
  #} else if(inherits(pcobj, 'PCA')) {
  #  nobs.factor <- sqrt(nrow(pcobj$call$X))
  #  d <- unlist(sqrt(pcobj$eig)[1])
  #  u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
  #  v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
  #  exp_var <- d[choices]^2/d.total
  #} else if(inherits(pcobj, "lda")) {
  #  nobs.factor <- sqrt(pcobj$N)
  #  d <- pcobj$svd
  #  u <- predict(pcobj)$x/nobs.factor
  #  v <- pcobj$scaling
  #  exp_var <- d[choices]^2/sum(d^2)
  #} else 
  if(inherits(pcobj, "rda") | inherits(pcobj, "cca")) {
    s = summary(pcobj)
    nobs.factor <- sqrt(nrow(s$sites) - 1)
    d <- sqrt(s$concont$importance[1, choices])
    u <- s$sites
    v <- s$species
    b <- s$biplot
    exp_var <- s$concont$importance[2, choices]
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, lda or rda')
  }
  # Observation
  choices <- pmin(choices, ncol(u))
  #u <- sweep(u, 2, d[choices]^obs.scale, FUN='*')
  df.u <- as.data.frame(u[,choices])
  
  # Variables
  #v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])
  
  # Biplot
  df.b <- as.data.frame(b[, choices])
  
  names(df.u) = names(df.v) = names(df.b) = c('xvar', 'yvar')
  
  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  
  # Scale the radius of the correlation circle so that it corresponds to 
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  
  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))
  
  # Scale biplot
  b.scale <- rowSums(v^2)
  df.b <- r * df.b / sqrt(max(b.scale))
  
  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- if(axes.lang == 'FR') paste0('ARD standardisée', choices) else paste0('standardized RDA', choices)
  } else {
    u.axis.labs <- if(axes.lang == 'FR') paste0('ARD', choices) else paste0('RDA', choices)
  }
  #change text for axis
  axis.txt = if(axes.lang == 'FR') 'var. expliquée)' else 'explained var.)'
  
  # Append the proportion of explained variance to the axis labels with text
  u.axis.labs <- paste(u.axis.labs, 
                       sprintf('(%0.1f%%', 
                               100 * exp_var),
                       axis.txt)
  
  # Score Labels
  if(!is.null(labels)) {
    df.u$labels <- labels
  }
  
  # Grouping variable
  if(!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable Names
  if(varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }
  
  # Biplot Names
  if(biname.abbrev) {
    df.b$biname <- abbreviate(rownames(b))
  } else {
    df.b$biname <- rownames(b)
  }
  
  # Variables for text label placement - species
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - name.adjust * sign(xvar)) / 2)
  
  # Variables for text label placement - biplot
  df.b$angle <- with(df.b, (180/pi) * atan(yvar / xvar))
  df.b$hjust = with(df.b, (1 - name.adjust * sign(xvar)) / 2)
  
  #= Base plot ==========================================================================================
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal() +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_vline(xintercept = 0, linetype="dashed") + 
    theme_minimal()
  
  if(!is.null(df.u$groups) & ellipse & is.null(obs.fill)) {
    g <- g + stat_ellipse(geom="polygon", aes(y = yvar, x = xvar, fill = groups), 
                          alpha = ellipse.alpha, 
                          show.legend = FALSE, 
                          level = ellipse.level,
                          type = ellipse.type)
  }

  if(var.arrows) {
    # Draw circle
    if(circle) 
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'), 
                         size = 1/2, alpha = 1/3)
    }
    
    # Get var.arrow color
    var.arrow.r = col2rgb(var.arrow.color)[1]/255 #
    var.arrow.g = col2rgb(var.arrow.color)[2]/255 #
    var.arrow.b = col2rgb(var.arrow.color)[3]/255 #
    var.arrow.labCol = rgb(var.arrow.r, var.arrow.g, var.arrow.b, arrow.alpha) #
    var.arrow.txtCol = rgb(var.arrow.r, var.arrow.g, var.arrow.b) #
    
    # Draw species arrows
    if(var.arrows){
      g <- g +
        geom_segment(data = df.v,
                     aes(x = 0, y = 0, xend = xvar*var.arrow.scaling, yend = yvar*var.arrow.scaling),
                     arrow = arrow(length = unit(1/2, 'picas')), 
                     color = var.arrow.labCol) #
    }
    # Get bi.arrow color
    bi.arrow.r = col2rgb(bi.arrow.color)[1]/255 #
    bi.arrow.g = col2rgb(bi.arrow.color)[2]/255 #
    bi.arrow.b = col2rgb(bi.arrow.color)[3]/255 #
    bi.arrow.labCol = rgb(bi.arrow.r, bi.arrow.g, bi.arrow.b, arrow.alpha) #
    bi.arrow.txtCol = rgb(bi.arrow.r, bi.arrow.g, bi.arrow.b) #
    
    # Draw biplot arrows
    if(bi.arrows){
      g <- g +
        geom_segment(data = df.b,
                     aes(x = 0, y = 0, xend = xvar*bi.arrow.scaling, yend = yvar*bi.arrow.scaling),
                     arrow = arrow(length = unit(1/2, 'picas')), 
                     color = bi.arrow.labCol)
    }
  }
  
  # Draw either labels or points
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)      
    }
  } else {
    fill = if(!is.null(obs.fill)) obs.fill else if(!is.null(groups)) groups else obs.color
    if(!is.null(obs.fill)){
      g <- g +  new_scale_fill() +
        scale_fill_manual(values = obs.pal) #grays
    }
    #shape = if(!is.null(obs.shape)) obs.shape else 21
    if(!is.null(df.u$groups)) {
      if(length(obs.size) == 1){ 
        if(length(fill) == 1) g <- g + geom_point(fill = fill, size = obs.size, color = obs.color, shape = 21, alpha = alpha)
        else g <- g + geom_point(aes(fill = fill), size = obs.size, color = obs.color, shape = 21, alpha = alpha)
      } else {
        if(length(fill) == 1) g <- g + geom_point(aes(size = obs.size), fill = fill, color = obs.color, shape = 21, alpha = alpha)
          else g <- g + geom_point(aes(fill = fill, size = obs.size), color = obs.color, shape = 21, alpha = alpha)
        #if(obs.binned){
        #  g <- g + scale_size_binned(guide = guide_bins())
        #}
      }
    } else {
      if(is.discrete(obs.size)){
        if(!is.null(fill)) g <- g + geom_point(fill = fill, size = obs.size, color = obs.color, shape = 21, alpha = alpha)
        else g <- g + geom_point(fill = obs.color, size = obs.size, color = obs.color, shape = 21, alpha = alpha)
      } else {
        if(!is.null(fill))  g <- g + geom_point(aes(size = obs.size), fill = fill, color = obs.color, shape = 21, alpha = alpha)
        else g <- g + geom_point(aes(size = obs.size), fill = obs.color, color = obs.color, shape = 21, alpha = alpha)
        #if(obs.binned){
        #  g <- g + scale_size_binned(guide = guide_bins())#, name = deparse(substitute(obs.size)))
        #}
      }
    }
  }
  ## Draw either labels or points
  #if(!is.null(df.u$labels)) {
  #  if(!is.null(df.u$groups)) {
  #    g <- g + geom_text(aes(label = labels, color = groups), size = labels.size)
  #  } else {
  #    g <- g + geom_text(aes(label = labels), size = labels.size)      
  #  }
  #} else {
  #  fill = if(!is.null(obs.fill)) obs.fill else groups
  #  shape = if(!is.null(obs.shape)) obs.shape else as.factor(21)
  #  print(fill)
  #  if(!is.null(df.u$groups)) {
  #    if(length(obs.size) == 1){ 
  #      if(length(fill) == 1 & length(shape) == 1) g <- g + geom_point(shape = shape, fill = fill, size = obs.size, color = obs.color, alpha = alpha)
  #      else if(length(shape) == 1) g <- g + geom_point(aes(fill = fill), shape = shape, size = obs.size, color = obs.color, alpha = alpha)
  #      else if(length(fill) == 1) g <- g + geom_point(aes(shape = shape), fill = fill, size = obs.size, color = obs.color, alpha = alpha)
  #      else g <- g + geom_point(aes(fill = fill, shape = shape), size = obs.size, color = obs.color, alpha = alpha)
  #    } else {
  #      if(length(fill) == 1 & length(shape) == 1) g <- g + geom_point(shape = shape, fill = fill, color = obs.color, alpha = alpha)
  #      else if(length(shape) == 1) g <- g + geom_point(aes(fill = fill), shape = shape, color = obs.color, alpha = alpha)
  #      else if(length(fill) == 1) g <- g + geom_point(aes(shape = shape), fill = fill, color = obs.color, alpha = alpha)
  #      else g <- g + geom_point(aes(fill = fill, shape = shape), color = obs.color, alpha = alpha)
  #    }
  #  } 
  #}
  # Draw either labels or points
  #if(!is.null(df.u$labels)) {
  #  if(!is.null(df.u$groups)) {
  #    g <- g + geom_text(aes(label = labels, color = groups), size = labels.size)
  #  } else {
  #    g <- g + geom_text(aes(label = labels), size = labels.size)      
  #  }
  #} else if(!is.null(df.u$groups)) {
  #  if(is.discrete(obs.size)){
  #    g <- g + geom_point(aes(fill = groups), size = obs.size, color = obs.color, shape = 21, alpha = alpha)
  #  } else if(obs.binned){
  #    g <- g + geom_point(aes(fill = groups, size = obs.size), color = obs.color, shape = 21, alpha = alpha) + 
  #      scale_size_binned(guide = guide_bins()) 
  #  } else {
  #    g <- g + geom_point(aes(fill = groups, size = obs.size), color = obs.color, shape = 21, alpha = alpha)
  #  }
  #} else if(is.discrete(obs.size)){
  #  g <- g + geom_point(fill = obs.color, size = obs.size, color = obs.color, shape = 21, alpha = alpha)
  #} else if(obs.binned){
  #  g <- g + geom_point(aes(size = obs.size), fill = obs.color, color = obs.color, shape = 21, alpha = alpha)
  #  g <- g + scale_size_binned(guide = guide_bins())#, name = deparse(substitute(obs.size)))
  #} else {
  #  g <- g + geom_point(aes(size = obs.size), fill = obs.color, color = obs.color, shape = 21, alpha = alpha)
  #}
  #g <- g + geom_point(aes(color = groups), alpha = alpha)
  #length.check = which(length(obs.size) == 1, length(obs.size) == 1) # if just one value, place outside aes
  #if(lenght(length.check) == 2) g <- g + geom_point(aes(fill = groups), color = obs.color, shape = 21, alpha = alpha, size = obs.size)
  #else if(length.check == 1)  g <- g + geom_point(aes(fill = groups, color = obs.color), size = obs.size, shape = 21, alpha = alpha) #only size length == 1, only size out
  #else if(length.check == 2) g <- g + geom_point(aes(fill = groups, size = obs.size), color = obs.color, shape = 21, alpha = alpha)
  #else{
  #  g <- g + geom_point(aes(fill = groups, size = obs.size, color = obs.color), shape = 21, alpha = alpha) #all in aes
  #  if(obs.binned){                                                                                        #and if binned = T; use binned scale
  #    discrete.check = which(!is.discrete(obs.size), !is.discrete(obs.size)) #check which is binnable and bin them
  #    if(length(discrete.check) == 2){
  #      g <- g + scale_size_binned(guide = guide_bins(), name = deparse(substitute(obs.size))) +
  #        scale_color_binned(guide = guide_bins(), name = deparse(substitute(obs.color)))
  #    } else if(discrete.check == 1){
  #      g <- g + scale_color_binned(guide = guide_bins(), name = deparse(substitute(obs.color))) #if binned = T; use binned scale
  #    } else if(discrete.check == 2){
  #      g <- g + scale_size_binned(guide = guide_bins(), name = deparse(substitute(obs.size))) #if binned = T; use binned scale
  #    }
  
  # Overlay a concentration ellipse if there are groups
  if(!is.null(df.u$groups) && ellipse) {
    #theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    #circle <- cbind(cos(theta), sin(theta))
    #
    #ell <- ddply(df.u, 'groups', function(x) {
    #  if(nrow(x) <= 2) {
    #    return(NULL)
    #  }
    #  sigma <- var(cbind(x$xvar, x$yvar))
    #  mu <- c(mean(x$xvar), mean(x$yvar))
    #  ed <- sqrt(qchisq(ellipse.prob, df = 2))
    #  data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'), 
    #             groups = x$groups[1])
    #})
    #names(ell)[1:2] <- c('xvar', 'yvar')
    #g <- g + geom_path(data = ell, aes(color = groups, fill = groups, group = groups), alpha = 0.5)
    if(is.null(obs.fill)){
      g <- g + stat_ellipse(geom="path", aes(y = yvar, x = xvar, color = groups), 
                            show.legend = FALSE, 
                            level = ellipse.level,
                            type = ellipse.type, 
                            lwd = ellipse.size)
    } else {
      if(ellipse.size == 'groups' | ellipse.size == 'gr'){
        print('gr')
        egroups = unlist(sapply(groups, function(x) paste0('ell', which(unique(groups) %in% x))))
        g <- g  + new_scale_size() + 
                  stat_ellipse(geom="path", 
                               aes(y = yvar, x = xvar, lty = groups, size = groups),
                               show.legend = ellipse.show.legend, 
                               level = ellipse.level,
                               type = ellipse.type
                               ) + 
                  scale_size_discrete(range = ellipse.size.range)
      } else {
        g <- g  + stat_ellipse(geom="path", 
                               aes(y = yvar, x = xvar), #, lty = groups),
                               show.legend = ellipse.show.legend, 
                               level = ellipse.level,
                               type = ellipse.type
                               )
      }
    }
    #if(!is.null(obs.fill)){
    #  g <- g +  scale_color_manual(values = c("#CCCCCC", "#969696", "#636363", "#252525")) #grays
    #}
  }
  
  # Label the variable arrows
  if(var.arrows){ 
    g <- g + 
      geom_text(data = df.v, 
                aes(label = varname, x = xvar*var.arrow.scaling, y = yvar*var.arrow.scaling, 
                    angle = angle, hjust = hjust), 
                color = var.arrow.txtCol, size = varname.size)
  }
  
  # Label the biplot arrows
  if(bi.arrows){ 
    g <- g + 
      geom_text(data = df.b, 
                aes(label = biname, x = xvar*bi.arrow.scaling, y = yvar*bi.arrow.scaling, 
                    angle = angle, hjust = hjust), 
                color = bi.arrow.txtCol, size = varname.size)
  }
  # Change the name of the legend
  # if(!is.null(groups)) {
  #   g <- g + scale_color_brewer(name = deparse(substitute(groups)), 
  #                               palette = 'Dark2')
  # }
  
  # TODO: Add a second set of axes
  
  return(g)
}
