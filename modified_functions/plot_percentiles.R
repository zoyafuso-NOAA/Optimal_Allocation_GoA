plot_percentiles <- function(values = NULL,
                             xs = NULL,
                             inner_color = "cadetblue1",
                             outer_color = "dodgerblue",
                             pt.cex = 1,
                             pt.colors = rev(grey.colors(11,
                                                         start = 0, 
                                                         end = 0.9)),
                             plot = T){
        
        if(is.null(values)) stop("Must supply values")
        
        temp_quants <- apply(X = values,
                             MARGIN = 1,
                             FUN = quantile,
                             probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
                             na.rm = T)
        
        if(plot) {
                polygon(y = c(temp_quants["5%",], rev(temp_quants["95%",])),
                        x = c(xs, rev(xs)),
                        col = outer_color,
                        border = F)
                polygon(y = c(temp_quants["25%",], rev(temp_quants["75%",])),
                        x = c(xs, rev(xs)),
                        col = inner_color,
                        border = F)
                
                points(y = temp_quants["50%", ],
                       x = xs,
                       pch = 15,
                       cex = pt.cex,
                       col = pt.colors )   
        }
        return(temp_quants)
}