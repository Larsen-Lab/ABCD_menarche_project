## Load libraries
library(cowplot)
library(patchwork)
library(data.table)
library(gratia)
require(gratia)
library(scales)
library(purrr)
library(tidyverse)

## ---- globals ----
font_size <- 16
theme_set(theme_classic(base_family = "sans", base_size = font_size))
theme_replace(axis.text = element_text(colour = "black", size = font_size))
line_size <- 3
point_size <- 2

### function to extract derivative, confidence interval, significance, and plot for GAMs ###
get_derivs_and_plot <- function(modobj, smooth_var, low_color = NULL, hi_color = NULL, xlabel = NULL, border_colors = NULL, flat_fill = FALSE) {
  if(packageVersion("gratia") < "0.9.2") {
    stop("Need to gratia package version ‘0.9.2’ to get derivative plot. Please update the package.")
  }
  this_font_size = font_size
  if (is.null(low_color)){low_color = "white"}
  if (is.null(hi_color)){hi_color = "grey20"}
  
  derv <- derivatives(modobj, term = sprintf('s(%s)', smooth_var), partial_match = TRUE)
  derv <- derv %>%
    mutate(sig = !(0 > .lower_ci & 0 < .upper_ci)) %>%
    mutate(sig_deriv = ifelse(sig == TRUE, .derivative, 0))
  
  if ("unordered" %in% names(modobj$model)) {
    derv$fac_levels <- factor(
      sapply(derv$.smooth, FUN = function(x) {str_match(x, levels(modobj$model$unordered))[!is.na(str_match(x, levels(modobj$model$unordered)))]}),
      levels = levels(modobj$model$unordered)
    )
  } else {
    derv$fac_levels <- factor(x = smooth_var, levels = smooth_var)
  }
  
  if (!is.null(border_colors)) {
    outline <- data.frame(
      fac_levels = factor(levels(derv$fac_levels), levels = levels(derv$fac_levels)),
      outline_color = border_colors
    )
  }
  
  if (flat_fill == TRUE) {
    dplot <- ggplot(data=derv) +
      geom_tile(aes(x = !!sym(smooth_var), y = .5, fill = abs(sig_deriv))) +
      scale_fill_gradient(low = "white", high = "black") +
      facet_grid(rows = "fac_levels") +
      theme(panel.spacing = unit(-.01, "cm"))
  } else {
    dplot <- ggplot(data=derv) +
      geom_tile(aes(x = !!sym(smooth_var), y = .5, fill = sig_deriv)) +
      facet_grid(rows = "fac_levels") +
      theme(panel.spacing = unit(-.01, "cm")) +
      scale_fill_gradient2(low = "#a6cee3", midpoint = 0, mid = "white",
                           high = "#ea4791", limits = c(min(derv$sig_deriv), max(derv$sig_deriv))) #fb9a99
  }
  
  if (!is.null(xlabel)) {
    dplot <- dplot + labs(x = xlabel, fill = sprintf("\u0394%s", smooth_var))
  } else {
    dplot <- dplot + labs(x = smooth_var, fill = sprintf("\u0394%s", smooth_var))
  }
  
  dplot <- dplot + theme(axis.title.y = element_blank(),
                         axis.text.y = element_blank(),
                         axis.text.x = element_text(size = this_font_size),
                         axis.line = element_blank(),
                         axis.ticks.y = element_blank(),
                         strip.background = element_blank(),
                         strip.text.y = element_blank(),
                         text = element_text(size = this_font_size),
                         legend.text = element_text(size = this_font_size),
                         axis.title = element_text(size = this_font_size),
                         legend.key.width = unit(1, "cm"),
                         legend.position = "bottom",
                         legend.background = element_blank(),
                         plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    guides(fill = guide_colorbar(reverse = F, direction = "horizontal", title.position = "left"))
  
  if (!is.null(border_colors)) {
    dplot <- dplot +
      geom_rect(data = outline,
                aes(color = border_colors,
                    ymin = -0.001,
                    ymax = 1.001,
                    xmin = min(derv[[smooth_var]]) - max(derv[[smooth_var]])/1000,
                    xmax = max(derv[[smooth_var]]) + max(derv[[smooth_var]])/1000),
                fill = "white", alpha = 0, size = 2) +
      scale_color_identity()
  } else {
    dplot <- dplot + geom_rect(aes(ymin=0, ymax=1, xmin=min(data), xmax=max(data)), color="black", fill="white", alpha=0)
  }
  
  return(dplot)
}

## ---- visualize_model ----
# Func to visualize GAM model outputs
visualize_model <- function(modobj,
                            smooth_var, 
                            int_var = NULL ,
                            group_var = NULL, 
                            plabels = NULL,
                            xlabel=NULL,
                            ylabel=NULL,
                            ymax = NA,
                            xlimits = NULL, 
                            quantile_limits = FALSE,
                            quants = c(.05,.95),
                            check_diagnostics = FALSE,
                            derivative_plot = FALSE,
                            flat_fill = FALSE,
                            difference_plot = FALSE,
                            show.data=TRUE, 
                            line_color = "black",
                            # colors for interaction plot
                            high_color = "#decbe4", # "#91bfdb"
                            high_line =  "#440154FF", #"#4575b4"
                            low_color = "#D3D3FF", #c7eae5" #fed9a6"# "#fc8d59"
                            low_line = "#3E1D77", #24868EFF #"#ff7f00" # "#f46d43"
                            side_density=FALSE,
                            show.legend=TRUE,
                            # --- new arguments for second derivative
                            second_deriv_x = NULL,
                            second_deriv_lower = NULL,
                            second_deriv_upper = NULL,
                            v_line = NULL){
  this_font_size = font_size
  if (any(class(modobj)=="gam")) {
    model <- modobj
    is_gamm<-FALSE
  } else if (class(modobj$gam)=="gam") {
    model <- modobj$gam
    is_gamm<-TRUE
  } else {
    stop("Can't find a gam object to plot")
  }
  s<-summary(model)
  df <- model$model
  ## Generate custom line plot
  np <- 10000 #number of predicted values
  
  theseVars <- attr(model$terms,"term.labels")
  varClasses <- attr(model$terms,"dataClasses")
  thisResp <- as.character(model$terms[[2]])
  
  all_terms <- attr(terms(df),"term.labels")
  
  
  if (!is.null(int_var)) {
    # We will produce an interaction plot
    if (!any(grepl(x=as.character(model$formula),pattern = int_var))) {
      warning("int_var not recognized in model formula!")
      return()
    }
    if (any(class(df[,int_var])=="character")){
      warning("int_var is of class `character` and should be either factor or numeric")
      return()
    }
    newdata <- with(df,expand.grid(x=seq(from=min(df[,smooth_var]),to=max(df[,smooth_var]),length.out=np),
                                   int.var = (if(is.numeric(df[,int_var])) quantile(df[,int_var],quants) else levels(df[,int_var]))))
    grad_fill <- ifelse(is.numeric(df[, int_var]), TRUE, FALSE)
    labs <- levels(df[,int_var])
    colnames(newdata)=c(smooth_var,int_var)
    
    if(is.numeric(df[,int_var])) {
      newdata[,"int.varq"] = factor(newdata[,int_var])
      # newdata <- newdata %>% mutate(int.varq = factor(int.var))
      cbar_vals<-quantile(df[,int_var],c(quants[1],.5,quants[2]))}
    
  } else {
    newdata <- with(df,expand.grid(x=seq(from=min(df[,smooth_var]),to=max(df[,smooth_var]),length.out=np)))
    colnames(newdata)=smooth_var
    
  }
  covariates <- all_terms[!(all_terms%in%c(smooth_var,int_var))]
  if (!is_empty(covariates)) {
    covariate_values <- df %>% dplyr::select(all_of(covariates)) %>% summarise_all(.funs = function(x){ifelse(is.numeric(x),yes = mean(x,na.rm=T),no = levels(factor(x))[2])})
    for (cov in covariates) {
      newdata[,cov]=covariate_values[,cov]
    }
  }

  
  # In case the model is not gaussian!
  # 1. Get predictions on the link scale
  p_link <- predict(object = model, newdata = newdata, se.fit = TRUE, type = "link")
  
  # 2. Extract the inverse link function from your model object (if it is gaussian, this is just the identity function)
  ilink <- model$family$linkinv
  
  # 3. Calculate intervals on link scale, then transform to the response scale
  pred <- cbind(newdata, 
                fit  = ilink(p_link$fit),
                selo = ilink(p_link$fit - 2 * p_link$se.fit),
                sehi = ilink(p_link$fit + 2 * p_link$se.fit))
  if (!is.null(group_var)) {
    pred[,group_var] = NA #these columns have to exist in the dataframe for plotting
  }
  pred[,thisResp] = 1 #these columns have to exist in the dataframe for plotting
  
  # colors for interaction plot
  # high_color = "#decbe4" # "#91bfdb"
  # high_line =  "#440154FF" #"#4575b4"
  # low_color = "#D3D3FF" #c7eae5" #fed9a6"# "#fc8d59"
  # low_line = "#3E1D77" #24868EFF #"#ff7f00" # "#f46d43"
  
  if (!is.null(int_var)) {
    if (grad_fill == T) {
      p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp, color = int_var))
      if (!is.null(group_var)) {
        if (show.data==TRUE) {
          p1<- p1 + geom_line(aes(group = !!group_var),alpha = .5)
        }
      }
      p1 <- p1 +
        scale_color_gradientn(colors = c(low_line,"grey90",high_line), breaks = unname(cbar_vals),name = "") +
        geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi", fill = "int.varq"),alpha = .8, linetype = 0) +
        scale_fill_manual(values = c(low_color,high_color),name=int_var) +
        geom_line(data = pred,aes_string(x = smooth_var, y = "fit",group = "int.varq"),linewidth = line_size) +
        labs(title = plabels)
      if (show.legend == FALSE) {
        p1 <- p1 + theme(legend.position = "none")
      }
      if (!is.null(xlabel)) {
        p1<-p1+xlab(xlabel)+ylab(ylabel)
      }
      if (quantile_limits==TRUE){
        q_vals = quantile(df[,thisResp],probs = quants)
        p1 <- p1 + ylim(q_vals[1],q_vals[2])
      }
      if (side_density==TRUE) {
        p1.m = p1+theme(plot.margin = unit(c(0,0,1,1),"mm"))
        dy <- ggplot(data=df,aes_string(x=thisResp))+geom_density(alpha=.25,fill="black")+coord_flip()
        
        dy <- dy + theme(axis.line.x = element_blank(),axis.text.x = element_text(colour = "white"),axis.ticks = element_blank(),axis.title = element_blank(),
                         axis.line.y = element_blank(),axis.text.y = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm"))
        dx <-  ggplot(data=df,aes_string(x=smooth_var))+geom_density(alpha=.25,fill="black")
        dx <- dx + theme(axis.line.y = element_line(color = "white"),axis.text.y = element_text(colour = "white"),
                         axis.ticks = element_blank(),axis.title.y = element_text(colour = "white"),
                         axis.title.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank(),
                         plot.margin = unit(c(0,0,-0.5,0),"cm"))
        if (!is.null(plabels)) {
          dx=dx+labs(title=plabels)
          p1.m <- p1.m + labs(title=NULL)
        }
        p1<-dx+p1.m+dy + plot_layout(design = c(area(1,1),area(2,1),area(2,2)),widths = c(1,.1),heights = c(.1,1))
        # p1 <- p1 + geom_xsidedensity(aes(y=stat(density),alpha=.25),show.legend = FALSE) +
        #   geom_ysidedensity(aes(x=stat(density),alpha=.25,color="black"),show.legend = FALSE)
      }
    } else {
      # black_color = scales::muted('blue',l=60,c=80)#  "#4daf4a" # 'blue' #"#1c405eff"
      # green_color = scales::muted('red',l=60,c=80)# "#984ea3" #"#285c27ff"
      # black_color = brewer_pal("qual",palette = "Set1")(6)[3]
      # green_color = brewer_pal("qual",palette = "Set1")(6)[4]
      # border_colors = c(black_color,green_color) # this is passed to the derivative plotter if needed.
      border_colors = brewer_pal("qual",palette = "Set1")(length(labs)) # this is passed to the derivative plotter if needed.
      
      p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp, color = int_var,fill=int_var))
      if (show.data==TRUE) {
        p1<-p1+ geom_point(alpha = 0.75,stroke = 0, size = point_size) 
      }
      if (show.data==TRUE) {
        p1<- p1 +  
          geom_point(alpha = .35,stroke = 0, size = point_size,show.legend = FALSE)
        # geom_hex(color=NULL)
        if (!is.null(group_var)) {
          if (show.data==TRUE) {
            p1<- p1 + geom_line(aes_string(group = group_var),alpha = .3)
          }
        } 
      } 
      # col_labs <- c(black_color,green_color)
      col_labs <- c("#440154FF", #1
                    "#24868EFF", #2
                    "#FDE725FF", #3
                    "red", #4
                    "blue", #5
                    "green", #6
                    "purple", #7
                    "pink", #8
                    "#324", #9
                    "#124", #10
                    "#156", #11
                    "#495", #12
                    "#044234", #13
                    "#440", #14
                    "#248", #15
                    "#248212", #16
                    "lightpink", #17
                    "purple3", #18
                    "deeppink", #19
                    "deeppink3", #20
                    "darkmagenta", #21
                    "mediumpurple2" #22
      ) #brewer_pal("qual",palette = "Accent")(length(labs))
      names(col_labs)=labs
      if (difference_plot == TRUE) {
        col_labs<-c(col_labs,"sig"="black","non_sig"=NA)
      }
      p1 <- p1 +
        #scale_color_brewer(type = "qual",palette = "Accent",direction = 1) +
        scale_color_manual(values = col_labs) +
        geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi", fill = int_var),alpha = .3, linetype = 0,show.legend=FALSE) +
        # scale_fill_brewer(type = "qual",palette = "Set1",direction = 1) +
        scale_fill_manual(values = col_labs) +
        geom_line(data = pred,aes_string(x = smooth_var, y = "fit",color = int_var),linewidth = line_size,show.legend = TRUE) 
      p1<-p1+theme(legend.position = c(.3,.79),legend.background = element_blank()) 
      
      if (!is.na(ymax)){
        p1 <- p1 + ylim(NA,ymax)
      }
      if (show.legend==FALSE) {
        p1<-p1+theme(legend.position = "none")
      } else {
        #p1<-p1+theme(legend.position = c(1,1))
      }
      if (quantile_limits==TRUE){
        
        q_vals = quantile(df[,thisResp],probs = quants)
        #note for later: do not cut off ribbon vals. Set this so that it is whichever is larger, qvals or pred$selo/pred$sehi
        p1 <- p1 + ylim(q_vals[1],q_vals[2])
      }
      if (side_density==FALSE) {
        # p1 <- p1 + geom_rug(sides = "br",alpha = .1)
      }
      p1 <- p1 +labs(title = plabels)+theme(legend.title = element_blank())
      if (!is.null(xlabel)) {
        p1<-p1+xlab(xlabel)+ylab(ylabel)
      }
      if (!is.null(xlimits)) {
        p1 <- p1 + xlim(xlimits[1], xlimits[2])
      }
      if (side_density==TRUE) {
        p1.m = p1+theme(plot.margin = unit(c(0,0,1,1),"mm"))
        dy <- ggplot(data=df,aes_string(x=thisResp,fill = int_var,color = int_var))+geom_density(alpha=.25,show.legend = F) +
          scale_fill_manual(values = col_labs) + scale_color_manual(values = col_labs) +
          coord_flip()
        dy <- dy + theme(axis.line.x = element_blank(),axis.text.x = element_text(colour = "white"),axis.ticks = element_blank(),axis.title = element_blank(),
                         axis.line.y = element_blank(),axis.text.y = element_blank(),
                         plot.margin = unit(c(0,0,0,0),"cm"))
        dx <-  ggplot(data=df,aes_string(x=smooth_var,fill = int_var,color = int_var))+geom_density(alpha=.25,show.legend = F)+
          scale_fill_manual(values = col_labs) + scale_color_manual(values = col_labs)
        dx <- dx + theme(axis.line.y = element_line(color = "white"),axis.text.y = element_text(colour = "white"),
                         axis.ticks = element_blank(),axis.title.y = element_text(colour = "white"),
                         axis.title.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank(),
                         plot.margin = unit(c(0,0,-0.5,0),"cm"))
        if (!is.null(plabels)) {
          dx=dx+labs(title=plabels)
          p1.m <- p1.m + labs(title=NULL)
        }
        if (derivative_plot == F) {
          # we can go ahead and assemble the plot
          p1<-dx+p1.m+dy + patchwork::plot_layout(design = c(area(1,1),area(2,1),area(2,2)),widths = c(1,.1),heights = c(.1,1))
        }
      }
      
    }
  } else {
    
    # No interaction variable, just produce a single line plot
    # int_var = "" # This is a lazy solution to making the existing code workable with no int_var.
    # thisPred <- data.frame(init = rep(0,np))
    
    
    # df <- df %>%
    #   gratia::add_partial_residuals(model)
    # df$partial_resids <- unlist(df[,grep(x=names(df),pattern = "s(",fixed = T)])
    
    p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp))
    if (!is.null(v_line)) {
      p1 <- p1 +
        geom_vline(xintercept = v_line, 
                   linetype = "dashed",
                   color = "gray90")
    }
    if (show.data==TRUE) {
      p1<- p1 +  
        geom_point(alpha = .3,stroke = 0, size = point_size,color = line_color)
      # geom_hex(show.legend = TRUE) + scale_fill_gradient(low="white",high=line_color,limits = c(1, 9), oob = scales::squish)
    } 
    if (!is.null(group_var)) {
      if (show.data==TRUE) {
        #cat("adding lines")
        p1<- p1 + geom_line(aes_string(group = group_var),alpha = .5)
      }
    }
    p1 <- p1 + geom_ribbon(data = pred,aes_string(x = smooth_var ,y=thisResp, ymin = "selo",ymax = "sehi"),fill = "gray", alpha = 1, linetype = 0) +
      geom_line(data = pred,aes_string(x = smooth_var, y = "fit"),linewidth = line_size,color=line_color)
    
    # --- ADD SECOND DERIVATIVE POINT HERE ---
    if (!is.null(second_deriv_x) & !is.null(second_deriv_lower)) {
      newdata.sd <- newdata[1,]
      newdata.sd[1,smooth_var]=second_deriv_x
      pred.y <- predict(object=model, newdata = newdata.sd, type = "response")
      y <- as.numeric(pred.y[1])
      second_deriv_df <- data.frame(
        x = second_deriv_x,
        y = y,
        lower = second_deriv_lower,
        upper = second_deriv_upper
      )
      
      p1 <- p1 +
        geom_linerange(data = second_deriv_df,
                      aes(y = y, xmin = lower, xmax = upper),
                      color = "gray40", size = 1,
                      #linetype = "dashed",
                      inherit.aes = FALSE) +
        geom_point(data = second_deriv_df, aes(x = x, y = y), 
                   fill = "gray80",
                   colour = "black", 
                   size = 3,
                   shape = 21,
                   position = position_nudge(y = 0),
                   inherit.aes = FALSE)
    }
    if (show.data==TRUE) {
      p1<-p1+ geom_point(alpha = 0.75,stroke = 0, size = point_size) 
    }
    if (side_density==FALSE) {
      # p1 <- p1+geom_rug()
    }
    p1 <- p1 + labs(title = plabels)
    if (!is.na(ymax)){
      p1 <- p1 + ylim(NA,ymax)
    }
    if (quantile_limits==TRUE){
      
      q_vals = quantile(df[,thisResp],probs = quants)
      p1 <- p1 + ylim(q_vals[1],q_vals[2])
    }
    if (!is.null(xlabel)) {
      p1<-p1+xlab(xlabel)+ylab(ylabel)
    }
    if (side_density==TRUE) {
      p1.m = p1+theme(plot.margin = unit(c(0,0,1,1),"mm"))
      dy <- ggplot(data=df,aes_string(x=thisResp))+geom_density(alpha=.25,fill=line_color,color=line_color)+coord_flip()
      
      dy <- dy + theme(axis.line.x = element_blank(),axis.text.x = element_text(colour = "white"),axis.ticks = element_blank(),axis.title = element_blank(),
                       axis.line.y = element_blank(),axis.text.y = element_blank(),
                       plot.margin = unit(c(0,0,0,0),"cm"))
      dx <-  ggplot(data=df,aes_string(x=smooth_var))+geom_density(alpha=.25,fill=line_color,color=line_color)
      dx <- dx + theme(axis.line.y = element_line(color = "white"),axis.text.y = element_text(colour = "white"),
                       axis.ticks = element_blank(),axis.title.y = element_text(colour = "white"),
                       axis.title.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank(),
                       plot.margin = unit(c(0,0,-0.5,0),"cm"))
      if (!is.null(plabels)) {
        dx=dx+labs(title=plabels)
        p1.m <- p1.m + labs(title=NULL)
      }
      if (derivative_plot == F) {
        # we can go ahead and assemble the plot
        p1<-dx+p1.m+dy + patchwork::plot_layout(design = c(area(1,1),area(2,1),area(2,2)),widths = c(1,.1),heights = c(.1,1))
      }
    }
  }
  
  if (difference_plot == TRUE) {
    ## First attempt with gratia is not quite right since it does not incorporate the mean differences between groups, only the centered smooths.
    
    # diffs <-gratia::difference_smooths(model = model, smooth = sprintf("s(%s)",smooth_var),n=1000)
    # diffs<-diffs %>% mutate(significant = factor(!between(x = 0,lower = lower,upper = upper),levels = c("TRUE","FALSE"),labels = c("sig","non_sig")))
    # diffs$yvalue <- min(pred$selo)*.5
    # p1 <- p1 + geom_tile(data=diffs,aes_string(x=smooth_var,fill="significant",y="yvalue"),
    #                      show.legend = FALSE,
    #                      color=NA,
    #                      height=(pred$sehi[[1]]-pred$selo[[1]])/5)+
    #   guides(fill="none")
    
    p_data <- mgcv::plot.gam(model)
    plot_df <- data.frame(x=p_data[[2]]$x,se=p_data[[2]]$se,fit=p_data[[2]]$fit)
    plot_df <- plot_df %>% mutate(se_upper=fit+se,se_lower=fit-se)%>%
      mutate(significant = factor(!between(x=0,lower=se_lower,upper = se_upper),levels = c("TRUE","FALSE"),labels = c("sig","non_sig")))
    names(plot_df)[names(plot_df)=="x"] = smooth_var
    plot_df$yvalue <- min(pred$selo)*.5
    
    names(plot_df)[grep(names(plot_df),pattern = 'x')]=smooth_var
    
    p1 <- p1 + geom_tile(data=plot_df,aes_string(x=smooth_var,fill="significant",y="yvalue"),
                         show.legend = FALSE,
                         color=NA,
                         height=(pred$sehi[[1]]-pred$selo[[1]])/5)+
      guides(fill="none")
    
  }
  if (derivative_plot == T) {
    # We will add a bar that shows where the derivative is significant.
    
    if (side_density==TRUE) {
      # we need to combine the derivative plot with the density plots and scatter plot
      # use the version of the scatter plot from the side_density portion of the code
      p1 <- p1.m
    }
    
    # First make some adjustments to the line plot.
    p1<- p1+theme(text = element_text(size=this_font_size),
                  axis.text = element_text(size = this_font_size),
                  axis.title.y = element_text(size = this_font_size),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  legend.text = element_text(size = this_font_size),
                  # legend.title = element_text(size = this_font_size),
                  axis.title = element_text(size = this_font_size),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "transparent",colour = NA),
                  plot.margin = unit(c(.2, .2, 0, .2), "cm")) #Top, left,Bottom, right
    
    # Now add the plots using the derivative plotting function
    if (!is.null(int_var)){
      if (any(grepl(x = row.names(s$s.table),pattern =  ":") & grepl(x=row.names(s$s.table),pattern = as.character(int_var)))) {
        # Refit the model as by-factor smooths rather than a by-factor interaction (change from ordered factor to unordered)
        f<-formula(model) # current formula
        # DON'T call family(model) yet - it causes errors with gamm extended families
        df$unordered <- factor(df[,int_var],ordered = F) # Change the factor to unordered
        model$model$unordered <- df$unordered
        fterms <- terms(f)
        v<-attr(fterms,"variables")
        smooth_locations <- grep(pattern = smooth_var,x = v)
        int_locations <- grep(pattern = int_var,x=v)
        drop_idx <- smooth_locations[!(smooth_locations %in% int_locations)]
        if (length(drop_idx)>0) {
          new_terms <- drop.terms(fterms, dropx = drop_idx, keep.response = TRUE)
          drop_formula <- formula(new_terms) # Formula without the main effect smooth.
          drop_string <- deparse(drop_formula)
          new_formula <- formula(paste(gsub(x=drop_string,pattern = int_var,replacement = "unordered"),collapse = " "))
          
          if (is_gamm==TRUE) {
            # For gamm objects, just use the original model to avoid extended.family issues
            this_mod <- model
          } else {
            # Only extract family for regular gam objects
            fam <- family(model)
            this_mod <- gam(formula = new_formula, data = df, family = fam)
          }
        } else {
          this_mod <- model
        }
        
        num_levels <- length(levels(df$unordered))
        
        if (!exists("border_colors")) {
          border_colors=NULL
        }
        
        if(!is.null(xlabel)){
          this_d <- get_derivs_and_plot(modobj = this_mod,smooth_var = smooth_var,border_colors = border_colors,flat_fill = flat_fill,xlabel = xlabel)
        } else{
          this_d <- get_derivs_and_plot(modobj = this_mod,smooth_var = smooth_var,border_colors = border_colors,flat_fill = flat_fill)
        }
        legend <- get_legend(this_d)
        plotlist <- list(p1,this_d+theme(legend.position = "none"))
        
        if (side_density==TRUE) {
          # Assemble using `patchwork` rather than `cowplot` due to the complexity of the figure.
          this_d <- this_d + theme(legend.position = "none")
          final_plot <- dx+p1+dy+plot_spacer()+this_d + 
            plot_layout(design = c(area(1,1),area(2,1),area(2,2),area(3,1),area(4,1)),widths = c(1,.1),heights = c(.1,1+.1*num_levels,-.11,.1*num_levels))
        } else {
          pg<-cowplot::plot_grid(rel_heights = c(2*num_levels,num_levels),plotlist = plotlist,align = "v",axis = "lr",ncol = 1)
          final_plot <- cowplot::plot_grid(pg,legend,rel_heights = c(1,.1),ncol = 1)
        }
        
      }
    } else {
      # No need to split
      if (!is.null(xlabel)) {
        this_d <- get_derivs_and_plot(modobj = model,smooth_var = smooth_var,xlabel = xlabel,border_colors = line_color,flat_fill = flat_fill)
      } else{
        this_d <- get_derivs_and_plot(modobj = model,smooth_var = smooth_var,border_colors = line_color,flat_fill = flat_fill)
      }
      if (side_density==TRUE) {
        # Assemble using `patchwork` rather than `cowplot` due to the complexity of the figure.
        this_d <- this_d + theme(legend.position = "none")
        final_plot <- dx+p1+dy+plot_spacer()+this_d + 
          plot_layout(design = c(area(1,1),area(2,1),area(2,2),area(3,1),area(4,1)),widths = c(1,.1),heights = c(.1,1,-.11,.1))
      } else {
        scatter <- list(p1)
        bar <- list(this_d+theme(legend.position = "none"))
        legend <- get_legend(this_d+theme(legend.position = "bottom",legend.direction = "horizontal"))
        allplots <- c(scatter,bar)
        pg<-cowplot::plot_grid(rel_heights = c(1,.35),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
        final_plot <- cowplot::plot_grid(pg,legend,rel_heights = c(1,.15),ncol = 1)
      }
    }
    
  }    else {
    # No derivative plot
    p1<- p1+theme(text = element_text(size=font_size),
                  axis.text = element_text(size = font_size),
                  legend.text = element_text(size = font_size),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  plot.background = element_blank())
    final_plot<-p1
  }
  
  # print(final_plot)
  if (check_diagnostics == T) {
    cp <- appraise(model,method="simulate")
    print(cp)
  }
  
  return(final_plot)
}