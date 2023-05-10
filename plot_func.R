#!/usr/bin/env Rscript
#@uthor : Thomas NEFF

# lib ---------------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(edgeR)
library(ggnewscale)

# fun ---------------------------------------------------------------------
ggplotBCV <- function(y,col.raw="grey70",col.common="red",col.trend="blue",
                      show.genes=T,genes=NULL,pval=0.05,col.pval="FDR",
                      sel=list(name="Exp.",size=2,shape=21,colors=c("red","white"),labels=c("up","down")),
                      pts=list(size=2,alpha=0.1,shape=20),theme=theme_bw(),show.plot=T,
                      x.lab="Average log(CPM)",y.lab="Biological coefficient of variation",
                      save_plot=list(file="plot_bcv",fmt="png",width=600,height=550)){
  
  # Compute AveLogCPM if not found in y
  A <- eval(y)$AveLogCPM
  if(is.null(A)) A <- aveLogCPM(y$counts, offset=getOffset(y))
  
  # Points to determine y axis limits
  disp <- getDispersion(eval(y))
  if(is.null(disp)) stop("No dispersions to plot")
  if(attr(disp,"type")=="common") disp <- rep_len(disp, length(A))
  
  df <- data.frame(cbind(rownames(eval(y)@.Data[[1]]),A,eval(y)$tagwise.dispersion,
        eval(y)$common.dispersion,eval(y)$trended.dispersion),row.names = 1)
  colnames(df) <- c("A","tagwise","common","trended")
  df <- data.frame(merge(x = df, y = eval(y)$tt, by = 'row.names', all = TRUE),row.names = "Row.names")
  df <- df[,c("A","tagwise","common","trended","logFC","FDR")]

  df <- data.frame(apply(X = df, MARGIN = 2, FUN = function(x) as.double(as.character(x))))
  rownames(df) <- rownames(eval(y)@.Data[[1]])
  
  colors <- c("Tagwise"=col.raw,"Common"=col.common,"Trend"=col.trend)
  
  df$dir <- ifelse(sign(df$logFC)==1,"up","down")

  p <- ggplot() +
    geom_point(data = df,mapping = aes(x = A, y = sqrt(tagwise),color=col.raw),size=pts$size,alpha=pts$alpha,shape=pts$shape) +
    scale_color_manual(name="",values = c(col.raw),label="Tagwise",aesthetics = "color") + new_scale_color() +
    geom_line(data = df,mapping = aes(x = A, y = sqrt(trended), color="Trend")) +
    geom_line(data = df, mapping = aes(x = A, y = sqrt(common),color="Common")) +
    scale_color_manual(name="Mean",values = colors) + eval(theme) + ylab(y.lab) + xlab(x.lab)
    
  if (show.genes) {
    if (is.null(genes)) {
      df$test <- ifelse(test = df[eval(col.pval)] <= pval,"sign","no.sign")
    } else {
      df$test <- sapply(genes,function(x) ifelse(x %in% rownames(df),"sign","no.sign"))
    }
    p <- p +geom_point(data = df[df$test=="sign",],mapping = aes(x = A,y = sqrt(tagwise),fill=dir),size=sel$size,shape=sel$shape) +
      scale_fill_manual(name=sel$name,values = sel$colors,labels=sel$labels,aesthetics = "fill") +
      geom_text_repel(data = df[df$test=="sign",], mapping = aes(x = A, y = sqrt(tagwise),label=rownames(df[df$test=="sign",])))
  }
  
  if (as.character(save_plot$file)) {
    
    if (is.null(save_plot$fmt)) { save_plot$fmt <- "png" }
    
    if (is.null(save_plot$width)) { save_plot$width <- 600 }
    
    if (is.null(save_plot$height)) { save_plot$height <- 550 }
    
    png(filename = paste(save_plot$file,save_plot$fmt,sep = "."),
        width = save_plot$width, heigth = save_plot$heigth)
    p
    dev.off()
  }

  if (show.plot) { print(p) }
  
  return(p)
}

ggplotQLF <- function(y,theme=theme_bw(),show.plot=T,genes=NULL, show.genes=T, pval=0.05,
                      raw=list(size=1,alpha=0.3,shape=20,label="Raw",col="grey60"),
                      shr=list(size=1,alpha=1,shape=20,label="Squeezed",col="purple"),
                      sel=list(name="Exp.",col=c("red","red","white","white"),shape=c(22,21,21,22),label=c("Raw up","Squeezed up","Squeezed down","Raw down")),
                      trend=list(col="blue",label="Trended",size=0.5),
                      x.label="Average log2(CPM)",y.label="Quarter-Root Mean Deviance" ) {
  
  glmfit <- eval(y)$fit
  A  <-  glmfit$AveLogCPM 
  if ( is.null(A))  A  <-  aveLogCPM(glmfit) 
  s2  <-  glmfit$deviance  /  glmfit$df.residual.zeros
  if ( is.null (glmfit$var.post ))  {  stop ( "besoin d'exécuter glmQLFit avant plotQLDisp" )  }
  var.post <- glmfit$var.post
  var.prior <- glmfit$var.prior
  
  df <- data.frame(cbind(rownames(eval(y)@.Data[[1]]),A,s2,var.post,var.prior),row.names = 1)
  colnames(df) <- c("A","s2","var.post","var.prior")
  df <- merge(x = df, y = eval(y)$tt, by = 'row.names', all = TRUE)
  df <- df[,c("A","s2","var.post","var.prior","logFC","FDR")]
  
  df <- data.frame(apply(X = df, MARGIN = 2, FUN = function(x) as.double(as.character(x))))
  df$symbol <- rownames(eval(y)@.Data[[1]])
  df$dir <- ifelse(sign(df$logFC)==1,"up","down")
  
  df.s2 <- df[,!colnames(df) %in% "var.post"]
  colnames(df.s2) <- c("A","pos","var.prior","logFC","FDR","symbol","dir")
  df.s2$type.pos <- "s2"
  df.var.post <- df[,!colnames(df) %in% "s2"]
  colnames(df.var.post) <- c("A","pos","var.prior","logFC","FDR","symbol","dir")
  df.var.post$type.pos <- "var.post"
  
  df <- rbind(df.s2,df.var.post)
  df$type.pos <- as.factor(df$type.pos)
  df$type.pos <- factor(df$type.pos,levels = c("s2","var.post"))
  
  label.raw.shr <- c(raw$label,shr$label)
  size.raw.shr <- c(raw$size,shr$size)
  alpha.raw.shr <- c(raw$alpha,shr$alpha)
  col.raw.shr <- c(raw$col,shr$col)
  
  p <- ggplot() + eval(theme) +
    geom_point(data = df, 
               mapping = aes(x = A,y = sqrt(sqrt(pos)),colour=type.pos,
                             alpha=factor(type.pos),size=factor(type.pos)),shape=raw$shape) +
    scale_alpha_manual(name="",values = alpha.raw.shr,label=label.raw.shr) +
    scale_color_manual(name="",values = col.raw.shr,label=label.raw.shr) + 
    scale_size_manual(name="",values=size.raw.shr,label=label.raw.shr)
  
  p <- p + new_scale_color() + geom_line(data=df,aes(x = A, y = sqrt(sqrt(var.prior)),color=trend$col)) + 
    scale_color_manual(name="",values=trend$col,label=trend$label)
  
  p <- p + xlab(label = x.label) + ylab(label = y.label) 
  
  if (show.genes) {
    
    if (is.null(genes)) { genes <- df[df["FDR"] <= pval,]$symbol }
    
    if (any(genes %in% df$symbol)) { genes <- genes[genes %in% df$symbol] }
    
    if (length(genes)==0) { stop("ERROR: he p-value must be too lower") }
    
    df$type.pos.dir <- paste(df$type.pos,df$dir,sep = ".")
    df$type.pos.dir <- as.factor(df$type.pos.dir)
    df$type.pos.dir <- factor(df$type.pos.dir,levels = c("s2.up","var.post.up","var.post.down","s2.down"))
    
    p <- p + geom_point(data=df[df$symbol == genes,],
                        mapping = aes(x = A,y = sqrt(sqrt(pos)),fill=type.pos.dir,shape=type.pos.dir)) +
            scale_fill_manual(name=sel$name,values = sel$col, label=sel$label) +
            scale_shape_manual(name=sel$name,values = sel$shape,label=sel$label)
    
    p <- p + geom_text_repel(data = df[df$symbol == genes,], mapping = aes(x = A, y = sqrt(sqrt(pos)),label=genes))
  }

  if (show.plot) { print(p) }

  return(df)
  
}
