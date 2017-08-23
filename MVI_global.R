require(magrittr)
# Missing at Random Function ----------------------------------------------
# MCAR_generate <- function(data, var_mis = 0.2, mis_prop = 0.2) {
#   mis_var_idx <- sample(1:ncol(data), round(ncol(data)*var_mis))
#   mis_idx_mat <- data.frame(row_idx = NULL, col_idx = NULL)
#   for (idx_temp in mis_var_idx) {
#     idx_df <- data.frame(row_idx = sample(1:nrow(data), round(ncol(data)*var_mis)), col_idx = idx_temp)
#     mis_idx_mat <- rbind(mis_idx_mat, idx_df)
#   }
#   mis_idx_mat <- as.matrix(mis_idx_mat)
#   data_mis <- data
#   data_mis[mis_idx_mat] <- NA
#   return (list(data_mis = data_mis, mis_idx = mis_idx_mat))
# }

MCAR_generate <- function(data, mis_prop = 0.5) {
  all_idx <- which(data != Inf, arr.ind = T)
  rdm_idx <- sample(1:nrow(all_idx), round(nrow(all_idx)*mis_prop))
  slc_idx <- all_idx[rdm_idx, ]
  data_res <- data
  data_res[slc_idx] <- NA
  return(list(data_res = data_res, mis_idx = slc_idx))
}

# RSR calculation ---------------------------------------------------------
# RSR_calculate <- function(data_c, data_i) {
#   RSR <- sqrt(sum((data_c - data_i)^2)/sum((data_c - mean(unlist(data_c)))^2))
#   return(RSR)
# }



# MNAR imputation compare -------------------------------------------------
MNAR_generate <- function (data_c, mis_var = 0.5, var_prop = seq(.3, .6, .1)) {
  data_mis <- data_c
  if (is.numeric(mis_var)) var_mis_list <- sample(1:ncol(data_c), round(ncol(data_c)*mis_var))
  else if (is.character(mis_var)) var_mis_list <- which(colnames(data_c) %in% mis_var)
  for (i in 1:length(var_mis_list)) {
    var_idx <- var_mis_list[i]
    cur_var <- data_mis[, var_idx]
    cutoff <- quantile(cur_var, sample(var_prop, 1))
    cur_var[cur_var < cutoff] <- NA
    data_mis[, var_idx] <- cur_var
  }
  mis_idx_df <- which(is.na(data_mis), arr.ind = T)
  return (list(data_mis = data_mis, mis_idx_df = mis_idx_df))
}

# Scale and recover -------------------------------------------------------
scale_recover <- function(data, method = 'scale', param_df = NULL) {
  results <- list()
  data_res <- data
  if (!is.null(param_df)) {
    if (method=='scale') {
      data_res[] <- scale(data, center=param_df$mean, scale=param_df$std)
    } else if (method=='recover') {
      data_res[] <- t(t(data)*param_df$std+param_df$mean)
    }
  } else {
    if (method=='scale') {
      param_df <- data.frame(mean=sapply(data, function(x) mean(x, na.rm=T)), 
                             std=sapply(data, function(x) sd(x, na.rm=T)))
      data_res[] <- scale(data, center=param_df$mean, scale=param_df$std)
    } else {stop('no param_df found for recover...')}
  }
  results[[1]] <- data_res
  results[[2]] <- param_df
  return(results)
}



# Multiplot 4 ggplot2 -----------------------------------------------------
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
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

