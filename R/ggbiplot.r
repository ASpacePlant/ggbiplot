# 
#  ggbiplot.r
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
#'  Biplot for Principal Components using ggplot2
#'
#' @param pcobj               an object returned by prcomp() or princomp() or rda() (from the vegan package)
#' @param choices             which PCs to plot
#' @param scale               covariance biplot (scale = 1), form biplot (scale = 0). When scale = 1, the inner product between the variables approximates the covariance and the distance between the points approximates the Mahalanobis distance.
#' @param obs.scale           scale factor to apply to observations
#' @param var.scale           scale factor to apply to variables
#' @param pc.biplot           for compatibility with biplot.princomp()
#' @param groups              optional factor variable indicating the groups that the observations belong to. If provided the points will be colored according to groups
#' @param ellipse             draw a stat_ellipse() for each group?
#' @param ellipse.type        what type of stat_ellipse to draw : 't', 'norm' or 'euclid'
#' @param ellipse.level       the confidence level at which to draw an ellipse (default is 0.95), or, if type="euclid", the radius of the circle to be drawn.
#' @param ellipse.size        width of the path of the ellipse drawn (size aes of stat_ellipse)
#' @param ellipse.show.legend logical. Should this layer be included in the legends? NA, the default, includes if any aesthetics are mapped. FALSE never includes, and TRUE always includes.
#' @param ellipse.alpha       alpha transparency value for the fill aes of the ellipse
#' @param labels              optional vector of labels for the observations
#' @param labels.size         size of the text used for the labels
#' @param alpha               alpha transparency value for the points (0 = transparent, 1 = opaque)
#' @param circle              draw a correlation circle? (only applies when prcomp was called with scale = TRUE and when var.scale = 1)
#' @param var.axes            draw arrows for the variables?
#' @param axes.lang           in which language should the axis text be displayed : 'FR' for french, everything else is english for now
#' @param varname.size        size of the text for variable names
#' @param varname.adjust      adjustment factor the placement of the variable names, >= 1 means farther from the arrow
#' @param varname.abbrev      whether or not to abbreviate the variable names
#' @param obs.size            size aes of observation
#' @param obs.size.binned     weither to use scale_size_binned() or not
#' @param obs.size.name       to manually name the obsevation size legend
#' @param obs.color           color aes of observation; default 'black' help distinguishing the points imo
#' @param arrow.scaling       factor to which the length of the variables arrows will be scaled
#' @param arrow.color         color aes of the variables arrows
#' @param arrow.alpha         alpha transparency value for the variables arrows
#'
#' @return                    a ggplot2 plot
#' @export
#' @examples
#'   data(wine)
#'   wine.pca <- prcomp(wine, scale. = TRUE)
#'   print(ggbiplot(wine.pca, obs.scale = 1, var.scale = 1, groups = wine.class, ellipse = TRUE, circle = TRUE))
#'
ggbiplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE,
                     obs.scale = 1 - scale, var.scale = scale, 
                     groups = NULL, ellipse = FALSE, ellipse.level = 0.69, ellipse.show.legend = NA,
                     ellipse.type = 'norm', ellipse.size = 1, ellipse.alpha = 0.15,
                     labels = NULL, labels.size = 3, alpha = 1, 
                     var.axes = TRUE,  axes.lang = 'EN',
                     circle = FALSE, circle.prob = 0.69, 
                     varname.size = 3, varname.adjust = 1.5, 
                     varname.abbrev = FALSE, 
                     obs.size = 2, obs.color = 'black', obs.size.name = F, obs.fill = NULL, obs.pal = NULL,
                     arrow.scaling = 1, arrow.color = muted('red'), arrow.alpha = 1, 
                     ...)
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)

  stopifnot(length(choices) == 2)

  # Recover the SVD
 if(inherits(pcobj, 'prcomp')){
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
    exp_var <- d[choices]^2/sum(d^2)
  } else if(inherits(pcobj, 'princomp')) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$loadings
    exp_var <- d[choices]^2/sum(d^2)
  } else if(inherits(pcobj, 'PCA')) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
    exp_var <- d[choices]^2/d.total
  } else if(inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    exp_var <- d[choices]^2/sum(d^2)
  } else if(inherits(pcobj, "rda")) {
    nobs.factor <- sqrt(nrow(pcobj$CA$u) - 1)
    d <- sqrt(pcobj$CA$eig)
    u <- pcobj$CA$u
    v <- pcobj$CA$v
    exp_var <- summary(pcobj)$concont$importance[2, choices]
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, lda or rda')
  }

  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))

  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])

  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)

  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }

  # Scale the radius of the correlation circle so that it corresponds to 
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)

  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))

  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- if(axes.lang == 'FR') paste0('CP standardisée', choices) else paste0('standardized PC', choices)
  } else {
    u.axis.labs <- if(axes.lang == 'FR') paste0('CP', choices) else paste0('PC', choices)
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

  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)

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
                          type = ellipse.type, 
                          size = ellipse.size)
  }
  
  if(var.axes) {
    # Draw circle
    if(circle) 
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'), 
                         size = 1/2, alpha = 1/3)
    }
    
    # Get arrow color
    arrow.r = col2rgb(arrow.color)[1]/255 #
    arrow.g = col2rgb(arrow.color)[2]/255 #
    arrow.b = col2rgb(arrow.color)[3]/255 #
    arrow.labCol = rgb(arrow.r, arrow.g, arrow.b, arrow.alpha) #
    arrow.txtCol = rgb(arrow.r, arrow.g, arrow.b) #
    
    # Draw directions
    g <- g +
      geom_segment(data = df.v,
                   aes(x = 0, y = 0, xend = xvar*arrow.scaling, yend = yvar*arrow.scaling),
                   arrow = arrow(length = unit(1/2, 'picas')), 
                   color = arrow.labCol) #
  }
  # Draw either labels or points
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)      
    }
  } else {
    fill = if(!is.null(obs.fill)) obs.fill else groups
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
                            lwd = 1)
    } else {
      g <- g + stat_ellipse(geom="path", aes(y = yvar, x = xvar, lty = groups),
                            show.legend = FALSE, 
                            level = ellipse.level,
                            type = ellipse.type, 
                            lwd = 1)
    }
  }

  # Label the variable axes
  if(var.axes) {
    g <- g + 
    geom_text(data = df.v, 
              aes(label = varname, x = xvar*arrow.scaling, y = yvar*arrow.scaling, 
                  angle = angle, hjust = hjust), 
              color = arrow.txtCol, size = varname.size)
  }
  # Change the name of the legend for groups
  # if(!is.null(groups)) {
  #   g <- g + scale_color_brewer(name = deparse(substitute(groups)), 
  #                               palette = 'Dark2')
  # }

  # TODO: Add a second set of axes

  return(g)
}
