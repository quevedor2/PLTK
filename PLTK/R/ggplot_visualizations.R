#' Multiplot
#' @description As function states, plot multiple plots onto one figure 
#' 
#' @param ... ... (e.g., p and q) in this case are individual plots
#' @param plotlist  No idea?
#' @param cols  Number of columns
#' @param layout Matrix specificy layout (cols and rows)
#'
#' @return
#' @importFrom grid grid.newpage pushViewport viewport grid.layout
#' @export
#'
#' @examples 
#' multiplot(p, q) 
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' ggplot publication theme
#' @description ggplot theme for publication ready Plots: https://rpubs.com/Koundy/71792
#' 
#' @param base_size  Default=18
#' @param base_family Default='helvetica'
#'
#' @importFrom ggthemes theme_foundation
#' @return Theme for ggplot
#' @export
#'
#' @examples
#' data <- data.frame("A"=letters[1:10], "Val"=1:10)
#' ggplot(x, aes(x=A, y=Val)) +
#' geom_point() + 
#' theme_Publication()
theme_Publication <- function(base_size=18, base_family="helvetica") {
  (theme_foundation(base_size=base_size)
  # , base_family=base_family
  + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = NA),
          axis.title = element_text(size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(), 
          axis.line = element_line(colour="black", size = 0.5),
          axis.ticks = element_line(),
          # panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA, size = 5),
          # legend.position = "bottom",
          # legend.direction = "horizontal",
          # legend.key.size = unit(2, 'lines'),
          # legend.title=element_blank(),
#           legend.key.size= unit(0.2, "cm"),
#           legend.margin = unit(0, "cm"),
          # legend.title = element_text(face="italic"),
          strip.background=element_blank(),
          strip.text=element_text(size = rel(1)),
          plot.margin=unit(c(10,5,5,5),"mm")
#           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#           strip.text = element_text(face="bold")
  ))
}

# === Colour pallete ===
# Usage: ggplot(data, aes(x, y)) + scale_fill_Publication()
.scale_fill_Publication <- function(...){
  require(scales)
  discrete_scale("fill","Publication",
                 manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

.scale_colour_Publication <- function(...){
require(scales)
  discrete_scale("colour","Publication",
                 manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

# === Emulate GGPLOT colours ===
# Usage: gg_color_hue(num)
# Number will indicate how many colours you want to emulate
.gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
