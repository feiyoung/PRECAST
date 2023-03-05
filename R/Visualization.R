
## prepared-for-plot Tool function####################
coordinate_rotate <- function(pos, theta=0){# counter-clock rotation
  pos_new <- pos
  pos_new[,1] <- pos[,1]*cos(theta) - pos[,2]*sin(theta)
  pos_new[,2] <- pos[,1]*sin(theta) + pos[,2]*cos(theta)
  return(pos_new)
}  

get_indexList <- function(alist){
  nsample <- length(alist)
  nr <- 0
  indexList <- list()
  for(i in 1:nsample){
    indexList[[i]] <- (nr+1):(nrow(alist[[i]] )+nr)
    nr <- nr + nrow(alist[[i]] )
  }
  return(indexList)
}

get_sampleID <- function(XList){
  sampleID <- list()
  r_max <- length(XList)
  for(r in 1:r_max){
    sampleID[[r]] <- rep(r, nrow(XList[[r]]))
  }
  sampleID <- unlist(sampleID)
  return(sampleID)
}


firstup <- function(x) {
  ## First letter use upper capital
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


matlist2mat <- function(XList){
  # transfer a matrix list to a matrix stacked by rows.
  r_max <- length(XList)
  X0 <- XList[[1]]
  if(r_max>1){
    for(r in 2:r_max){
      X0 <- rbind(X0, XList[[r]])
    }
  }
  
  return(X0)
}




# Vector based plot function##############################
barPlot <- function(vec, ylabel='ARI', cols=NULL,...){
  # require(ggplot2)
  
  
  ## filter vec
  N <- length(vec)
  vec_use <-vec[!is.na(vec)]
  
  
  df_use <- data.frame(value=vec_use, 
                       Method=names(vec_use))
  df_use$Method <- factor(df_use$Method, levels=names(vec_use))
  
  
  
  
  ## CCor
  p1 <- ggplot(df_use, aes(x=aes_string('Method'), y=aes_string('value'), 
     fill=aes_string('Method'))) + 
    geom_bar(position = "dodge", stat="identity",width = 1, ...) + # , ...
    #geom_errorbar( aes(ymin=value-sd, ymax=value+sd), width=0.4, colour="orange",  size=1.3, position=position_dodge(.9)) + 
    #facet_grid(beta~Error , scales="fixed",labeller = label_bquote(beta == .(beta))) 
    labs(y=ylabel, x=NULL)+ 
    scale_x_discrete(breaks = NULL) 
  
  if(is.null(cols)){
    return(p1)
  }else{
    p1 + scale_fill_manual(values = cols)
  }
}


# Matrix based plot function -------------------------------------------------------
plot_RGB <- function(position, embed_3d, pointsize=2,textsize=15){
  
  # suppressMessages(require(ggplot2))
  
  info = as.data.frame(position)
  colnames(info) = c("sdimx","sdimy")
  
  
  r = (embed_3d[,1]-min(embed_3d[,1]))/(max(embed_3d[,1])-min(embed_3d[,1]))
  g = (embed_3d[,2]-min(embed_3d[,2]))/(max(embed_3d[,2])-min(embed_3d[,2]))
  b = (embed_3d[,3]-min(embed_3d[,3]))/(max(embed_3d[,3])-min(embed_3d[,3]))
  x =  info$sdimx
  y =  info$sdimy
  dat = data.frame(x,y,r,g,b)
  p1=ggplot(data=dat, aes(x=x, y=y, col=rgb(r,g,b))) +
    geom_point(size=pointsize) +
    scale_color_identity()+
    theme_void()+
    theme(plot.title = element_text(size = textsize),
          text = element_text(size = textsize),
          #axis.title = element_text(face="bold"),
          #axis.text.x=element_text(size = 22) ,
          legend.position = "bottom")
  
  p1
}

plot_scatter <- function (
    embed_use, meta_data, label_name, xy_names=c('tSNE1', 'tSNE2'), no_guides = FALSE, 
    cols = NULL, 
     point_size = 0.5, point_alpha=1, 
    base_size = 12, do_points = TRUE, do_density = FALSE, border_col='gray',
    legend_pos='right', legend_dir='vertical', nrow.legend=NULL) {
  # require(dplyr)
  # require(ggthemes)
  # require(ggrepel)
  # require(data.table)
  plt_df <- embed_use %>% data.frame() %>% cbind(meta_data) %>% 
    dplyr::sample_frac(1L)
  plt_df$given_name <- plt_df[[label_name]]
  
  if(is.null(cols)){
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    
    cluster <- as.vector(plt_df$given_name)
    ngrp <- length(unique(cluster))
    
      
    cols  <- gg_color_hue(ngrp)
  }
  
  
  plt <- plt_df %>% ggplot(aes_string(colnames(plt_df)[1],colnames(plt_df)[2], col = label_name, 
                                      fill = label_name)) + #  + theme_tufte(base_size = base_size, ticks= show_ticks)
    theme(axis.text.x=element_text(size=base_size, color=1),
          axis.text.y=element_text(size=base_size, color=1),
          axis.title.x = element_text(size=base_size+2, color='black'),
          axis.title.y = element_text(size=base_size+2, color='black'),
          strip.text =  element_text(size=base_size, color='black'),
          strip.background = element_rect(
            linetype = 'solid', color='gray3'
          ),
          legend.direction = legend_dir, legend.position = legend_pos,
          legend.text=element_text(size=base_size+1),
          legend.title=element_text(size=base_size+2),
          panel.background= element_rect(fill = 'white', color=border_col))+
    scale_color_manual(values = cols) + scale_fill_manual(values = cols) + 
    theme(plot.title = element_text(hjust = 0.5)) + labs(x = xy_names[1], 
                                                         y = xy_names[2])
  
  if (!is.null(nrow.legend)){
    plt <- plt + guides(color = guide_legend(nrow = nrow.legend,override.aes = list(stroke = 1,
                                                                                    alpha = 1, shape = 16, size = 4)))
  }else{
    plt <- plt + guides(color = guide_legend(override.aes = list(stroke = 1,
                                                                 alpha = 1, shape = 16, size = 4)),
                        alpha = "none")
  }
  
  
  if (do_points) 
    plt <- plt + geom_point( size = point_size, alpha=point_alpha)
  if (do_density) 
    plt <- plt + geom_density_2d()
  if (no_guides) 
    plt <- plt + guides(col = 'none', fill = 'none', alpha = 'none')
  
  return(plt)
}



volinPlot <- function(mat, ylabel='ARI', cols=NULL){
  # require(ggplot2)
  ## filter mat
  N <- nrow(mat)
  mat_use <- mat[,which(colSums(is.na(mat)) != N)]
  
  
  df_use <- data.frame(value=as.vector(mat_use), Method=rep(colnames(mat_use), each=N))
  df_use$Method <- factor(df_use$Method, levels=colnames(mat_use))
  
  p1 <- ggplot(df_use, aes_string(x = 'Method', y = 'value',
                           fill = 'Method')) + 
    geom_violin(aes_string(fill = 'Method' ), color = "transparent", alpha = 0.5) +
    geom_boxplot(outlier.alpha = 0, coef = 0, color = "gray40", width = 0.4) +
    labs(x = "", y = ylabel) +
    theme_classic() + theme( axis.text.x = element_blank() ) + theme(text = element_text(size=20))
  if(is.null(cols)){
    p1
  }else{
    pal <- cols
    p1+scale_fill_manual(values = pal, name = "") 
  }
  
}


boxPlot <- function(mat, ylabel='ARI', cols=NULL, ...){
  
  # require(ggplot2)
  ## filter mat
  N <- nrow(mat)
  mat_use <- mat[,which(colSums(is.na(mat)) != N)]
  
  
  df_use <- data.frame(value=as.vector(mat_use), Method=rep(colnames(mat_use), each=N))
  df_use$Method <- factor(df_use$Method, levels=colnames(mat_use))
  p1 <- ggplot(df_use, aes_string(x="Method", y="value",
                           fill="Method") ) +
    geom_boxplot(...) +  
    # facet_grid(beta~Error,scales= "fixed",labeller = label_bquote(beta == .(beta)) )+ 
    labs(x=NULL, y= ylabel) + scale_x_discrete(breaks = NULL)
  if(is.null(cols)){
    return(p1)
  }else{
    p1 + scale_fill_manual(values = cols)
  }
  
  
}



# Seurat based Plot function -----------------------------------------------------------

doHeatmap <- function(seu, features=NULL, cell_label='Cell type', grp_label = FALSE,
                      pt_size=4, grp_color=NULL, ...){
  # require(ggplot2)
  ngrp <- nlevels(Idents(seu))
  if(is.null(grp_color)){
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    grp_color <- gg_color_hue(ngrp)
  }
  
  
  Seurat::DoHeatmap(object = seu, features=features, group.colors = grp_color[1:ngrp], label = grp_label, ...) +
    guides(color = guide_legend(title=cell_label,override.aes = list(stroke = 1, 
                                                                     alpha = 1, shape = 16, size = pt_size, color=grp_color[1:ngrp])), 
           alpha =  "none")
}




featurePlot <- function(seu, feature=NULL, cols=NULL, pt_size=1, title_size=16, 
                        quant=0.5, assay='RNA', reduction="position"){
  dat <- as.data.frame(seu[[reduction]]@cell.embeddings)
  colnames(dat) <- c("Spatial_1", "Spatial_2")
  if(is.null(feature)) feature <- row.names(seu)[1]
  dat$Expression <- seu[[assay]]@scale.data[feature,]
  if(is.null(cols)) cols <- c("#0571B0",  "#CA0020")
  med <- quantile(seu[[assay]]@scale.data[feature,], quant)
  ggplot(data=dat, aes_string(x='Spatial_1', y='Spatial_2', color='Expression')) + geom_point(size=pt_size) +
    scale_colour_gradient2(
      low = cols[1],
      mid = "white",
      high = cols[2], midpoint = med) + mytheme_graybox() + 
    ggtitle(feature) + theme(title =element_text(size=title_size, color=1, face='italic'))
}



# Plot  Theme -------------------------------------------------------------------


doHeatmap.matrix <- function(corMat, cluster_orderd,legend_title='Cell type', grp_label = FALSE,
                             pt_size=4, grp_color=NULL){
  
  doHeatmap <- function(seu, features=NULL, cell_label='Cell type', grp_label = FALSE,
                        pt_size=4, grp_color=NULL, ...){
    # require(ggplot2)
    ngrp <- nlevels(Idents(seu))
    if(is.null(grp_color)){
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
      }
      grp_color <- gg_color_hue(ngrp)
    }
    
    
    Seurat::DoHeatmap(object = seu, features=features, group.colors = grp_color[1:ngrp], label = grp_label, ...) +
      guides(color = guide_legend(title=cell_label,override.aes = list(stroke = 1, 
                                                                       alpha = 1, shape = 16, size = pt_size, color=grp_color[1:ngrp])), 
             alpha =  "none")
  }
  
  # require(Seurat)
  seu <- CreateSeuratObject(counts=corMat)
  Idents(seu) <- factor(cluster_orderd, levels=1: max(cluster_orderd))
  
  seu[["RNA"]]@scale.data <- corMat
  
  doHeatmap(seu, features = row.names(seu), cell_label=legend_title, pt_size=pt_size,
            grp_color=grp_color, grp_label=grp_label)
  
}


# barPlot_enrich <- function(top_dat, source='Ont', term_name="Term", nlog10P='nlog10P',
#                            bar_width=0.8, base_size=20, font_family='serif', cols= ggthemes::canva_pal()(4)){
#   # source='Ont'; term_name="Term"; nlog10P='nlog10P'
#   require(ggplot2) # y=term_name,
#   order_idx <- order(top_dat[,nlog10P])
#   top_dat <- top_dat[order_idx,]
#   top_dat[, term_name] <- factor(top_dat[, term_name], levels=top_dat[order_idx,term_name])
#   p1 <- ggplot(data=top_dat, aes_string(x=term_name,y=nlog10P, fill=source)) +
#     scale_fill_manual(values=cols)+
#     geom_bar(position = "dodge", stat="identity",width =bar_width)+ coord_flip() +
#     theme_classic() + theme(text=element_text(size=base_size, family=font_family)) 
#   return(p1)
# }



mytheme_graybox <- function (base_size = 11, base_family = "", base_line_size = base_size/22,
                             base_rect_size = base_size/22, border_color = 'gray10', bg_fill='white')
{
  half_line <- base_size/2
  t <- theme(panel.background = element_rect(fill = bg_fill,
                                             colour = NA), panel.border = element_rect(fill = NA,
                                                                                       colour = border_color),
             #line = element_blank(), #rect = element_blank(),
             panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             text = element_text(family = base_family, face = "plain",
                                 colour = "black", size = base_size, lineheight = 0.9,
                                 hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(),
                                 debug = FALSE), axis.text = element_blank(), axis.title = element_blank(),
             axis.ticks.length = unit(0, "pt"), axis.ticks.length.x = NULL,
             axis.ticks.length.x.top = NULL, axis.ticks.length.x.bottom = NULL,
             axis.ticks.length.y = NULL, axis.ticks.length.y.left = NULL,
             axis.ticks.length.y.right = NULL, legend.box = NULL,
             legend.key.size = unit(1.2, "lines"), legend.position = "right",
             legend.text = element_text(size = rel(0.8)), legend.title = element_text(hjust = 0),
             strip.text = element_text(size = rel(0.8)), strip.switch.pad.grid = unit(half_line/2,
                                                                                      "pt"), strip.switch.pad.wrap = unit(half_line/2,
                                                                                                                          "pt"), panel.ontop = FALSE,
             panel.spacing = unit(half_line/2, "pt"), plot.margin = unit(rep(0.2,4), "lines"),
             plot.title = element_text(size = rel(1.2), hjust = 0,
                                       vjust = 1, margin = margin(t = half_line)), plot.title.position = "panel",
             plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(t = half_line)),
             plot.caption = element_text(size = rel(0.8), hjust = 1,
                                         vjust = 1, margin = margin(t = half_line)), plot.caption.position = "panel",
             plot.tag = element_text(size = rel(1.2), hjust = 0.5,
                                     vjust = 0.5), plot.tag.position = "topleft",
             complete = TRUE)
  #ggplot2:::ggplot_global$theme_all_null %+replace% t
  t
}


mytheme_void <- function (base_size = 11, base_family = "", base_line_size = base_size/22,
                          base_rect_size = base_size/22)
{
  half_line <- base_size/2
  t <- theme(panel.background = element_rect(fill = "white",
                                             colour = NA), panel.border = element_rect(fill = NA,
                                                                                       colour = "grey"),
             line = element_blank(), #rect = element_blank(),
             text = element_text(family = base_family, face = "plain",
                                 colour = "black", size = base_size, lineheight = 0.9,
                                 hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(),
                                 debug = FALSE), axis.text = element_blank(), axis.title = element_blank(),
             axis.ticks.length = unit(0, "pt"), axis.ticks.length.x = NULL,
             axis.ticks.length.x.top = NULL, axis.ticks.length.x.bottom = NULL,
             axis.ticks.length.y = NULL, axis.ticks.length.y.left = NULL,
             axis.ticks.length.y.right = NULL, legend.box = NULL,
             legend.key.size = unit(1.2, "lines"), legend.position = "right",
             legend.text = element_text(size = rel(0.8)), legend.title = element_text(hjust = 0),
             strip.text = element_text(size = rel(0.8)), strip.switch.pad.grid = unit(half_line/2,
                                                                                      "pt"), strip.switch.pad.wrap = unit(half_line/2,
                                                                                                                          "pt"), panel.ontop = FALSE, panel.spacing = unit(half_line,
                                                                                                                                                                           "pt"), plot.margin = unit(c(0, 0, 0, 0), "lines"),
             plot.title = element_text(size = rel(1.2), hjust = 0,
                                       vjust = 1, margin = margin(t = half_line)), plot.title.position = "panel",
             plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(t = half_line)),
             plot.caption = element_text(size = rel(0.8), hjust = 1,
                                         vjust = 1, margin = margin(t = half_line)), plot.caption.position = "panel",
             plot.tag = element_text(size = rel(1.2), hjust = 0.5,
                                     vjust = 0.5), plot.tag.position = "topleft",
             complete = TRUE)
  #ggplot2:::ggplot_global$theme_all_null %+replace% t
  t
}


mytheme <- function(legend.direction = "horizontal",
                    legend.position = "bottom", type = 'rect'){
  
  if(type=='rect'){
    th <- theme(axis.text.x=element_text(size=16, color=1, face='plain'),
                axis.text.y=element_text(size=16, color=1, face='plain'),
                axis.title.x = element_text(size=18, color='black', face='plain'),
                axis.title.y = element_text(size=18, color='black',face='plain'),
                strip.text =  element_text(size=16, color='black', face='plain'),
                strip.background = element_rect(
                  linetype = 'solid', color='gray3'
                ),
                legend.direction = legend.direction, legend.position = legend.position,
                legend.text=element_text(size=17, face='plain'),
                legend.title=element_text(size=18, face='bold'),
                panel.background= element_rect(fill = 'white', color='gray'))
  }
  return(th)
}



# Other plots -------------------------------------------------------------


drawFigs <- function(pList, layout.dim = NULL, common.legend=FALSE,legend.position='right',  ...){
  if(!is.list(pList)) stop('drawFigs: pList must be a list!')
  
  if(is.null(layout.dim) && length(pList)>1){
    layout.dim <- c(2, round(length(pList)/2) )
  }
  if(is.null(layout.dim) && length(pList) == 1){
    layout.dim <- c(1,1)
  }
  ggpubr::ggarrange(plotlist = pList, ncol = layout.dim[2],
                    nrow = layout.dim[1], common.legend = common.legend,
                    legend = legend.position, ...)
  
}

chooseColors <- function(palettes_name= c("Nature 10", "Light 13", "Classic 20", "Blink 23", "Hue n"), n_colors = 7,
                         alpha=1, plot_colors=FALSE){
  
  # require(colorspace)
  palettes_name <- match.arg(palettes_name)
  colors <- if(palettes_name == "Classic 20"){
    # require(ggthemes)
    # palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
    pal1 <- tableau_color_pal(palettes_name)
    pal1(n_colors)
  }else if(palettes_name == "Nature 10"){
    cols <- c("#E04D50", "#4374A5", "#F08A21","#2AB673", "#FCDDDE",
              "#70B5B0", "#DFE0EE" ,"#DFCDE4", "#FACB12", "#f9decf")
    cols[1:n_colors]
  }else if(palettes_name == "Blink 23"){
    cols <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#FE0092", "#bc9000","#4ffc00", "#00bcac", "#0099cc",
              "#D35400", "#00eefd", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8",
              "#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#007CC8")
    cols[1:n_colors]
  }else if(palettes_name == "Light 13"){
    cols <-c( "#FD7446" ,"#709AE1", "#31A354","#9EDAE5",
              "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,"#DADAEB" ,
              "#FF9896","#91D1C2", "#C7E9C0" ,
              "#6B6ECF", "#7B4173" )
    cols[1:n_colors]
  }else if(palettes_name == "Hue n"){
    gg_color_hue(n_colors)
  }else{
    stop(paste0("chooseColors: check palettes_name! Unsupported palettes_name: ", palettes_name))
  }
  #require(colorspace)
  colors_new = adjust_transparency(colors,   alpha = alpha)
  if(plot_colors){
    barplot(rep(1, length(colors_new)), axes = FALSE, space = 0, col = colors_new)
  }
  
  return(colors_new)
}

