#Graph functions
plothistlog <- function(col_name,colname_log, df, ylabtext, color_1, breaks_colname = "Sturges",breaks_colname_log,  density_plot = TRUE) {
  
    par(mfrow=c(2,2))
    hist(df[[col_name]],
         main = paste("Histogram of", col_name),
         xlab = col_name,
         col = color_1,
         freq = FALSE,
         breaks = breaks_colname
    )
    boxplot(df[[col_name]],
            main= paste("Boxplot of", col_name),
            xlab= col_name,
            ylab=ylabtext,
            col=color_1,
            horizontal=TRUE)
    ## part 2: log plots
    hist(df[[colname_log]],
         main = paste("Histogram of", colname_log),
         xlab = col_name,
         col = color_1,
         freq = FALSE,
         breaks = breaks_colname_log
    )
    boxplot(df[[colname_log]],
            main= paste("Boxplot of", colname_log),
            xlab= colname_log,
            ylab=ylabtext,
            col=color_1,
            horizontal=TRUE)
    par(mfrow=c(1,1))
  }
















plothist <- function(col_name, df, ylabtext, color_1, breaks = "Sturges", density_plot = TRUE) {
  
  if(density_plot == TRUE){
  
  density_col <- density(df[[col_name]],
                         kernel="gaussian")
  
  # Min-Max x axis
  min_x_d <- min(density_col$x)
  max_x_d <- max(density_col$x)
  
  # Min-Max y axis 
  min_y_d <- min(density_col$y)
  max_y_d <- max(density_col$y)
  
  par(mfrow=c(1,2))
  hist(df[[col_name]],
       main = paste("Histogram of", col_name),
       xlab = col_name,
       col = color_1,
       freq = FALSE,
       xlim = c(min_x_d, max_x_d),
       ylim = c(min_y_d, max_y_d),
       breaks = breaks
  )
  lines(density_col,
        ylab="Density",
        col="blue",
        lwd=2)
  boxplot(df[[col_name]],
          main= paste("Boxplot of", col_name),
          xlab= col_name,
          ylab=ylabtext,
          col=color_1,
          horizontal=TRUE)
  par(mfrow=c(1,1))
  }
  # I donot want to plot the density plot
  else{
    par(mfrow=c(1,2))
    hist(df[[col_name]],
         main = paste("Histogram of", col_name),
         xlab = col_name,
         col = color_1,
         freq = FALSE,
         breaks = breaks
    )
    boxplot(df[[col_name]],
            main= paste("Boxplot of", col_name),
            xlab= col_name,
            ylab=ylabtext,
            col=color_1,
            horizontal=TRUE)
    par(mfrow=c(1,1))
  }
  
  
}



plothist_factor <- function(col_name, col_factor, df, ylabtext, color_2, color_3, density_plot = TRUE) {
  
  if(density_plot == TRUE){
  
  d_Yes <- density(df[[col_name]][df[[col_factor]]=="Yes"],
                   kernel="gaussian")
  
  min_x_d_Yes <- min(d_Yes$x)
  max_x_d_Yes <- max(d_Yes$x)
  min_y_d_Yes <- min(d_Yes$y)
  max_y_d_Yes <- max(d_Yes$y)
  
  d_No <- density(df[[col_name]][df[[col_factor]]=="No"],
                  kernel="gaussian")
  
  min_x_d_No <- min(d_No$x)
  max_x_d_No <- max(d_No$x)
  min_y_d_No <- min(d_No$y)
  max_y_d_No <- max(d_No$y)
  
  min_x <- min(c(min_x_d_Yes,min_x_d_No))
  max_x <- max(c(max_x_d_Yes,max_x_d_No))
  min_y <- min(c(min_y_d_Yes,min_y_d_No))
  max_y <- max(c(max_y_d_Yes,max_y_d_No))
  
  par(mfrow=c(1,2))
  boxplot(df[[col_name]]~df[[col_factor]],
          main= paste("Boxplot of ", col_name ,"in terms of Private College"),
          ylab= ylabtext,
          xlab=col_name,
          col=c(color_3,color_2),
          horizontal=TRUE)
  plot(c(min_x,max_x),
       c(min_y,max_y),
       xlab= col_name,
       ylab= "Density",
       main= paste("Kernel density of ", col_name, "in terms of \n Private College"),
       type= "n")
  lines(d_Yes$x,d_Yes$y,col=color_2,lwd=2)
  lines(d_No$x,d_No$y,col=color_3,lwd=2)
  legend("topright", inset=.05, title="Private College",
         c("Yes","No"), fill = c("seagreen2","orange2"), horiz=TRUE)
  par(mfrow=c(1,1))
  }else{
    
    d_Yes <- hist(df[[col_name]][df[[col_factor]]==1],
                    plot = FALSE)
    
    min_x_d_Yes <- min(d_Yes$breaks)
    max_x_d_Yes <- max(d_Yes$breaks)
    
    
    d_No <- hist(df[[col_name]][df[[col_factor]]==0],
                    plot=FALSE)
    
    min_x_d_No <- min(d_No$breaks)
    max_x_d_No <- max(d_No$breaks)
    
    min_x <- min(c(min_x_d_Yes,min_x_d_No))
    max_x <- max(c(max_x_d_Yes,max_x_d_No))
    
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
    boxplot(df[[col_name]]~df[[col_factor]],
            main= paste("Boxplot of ", col_name ,"in terms of ", ylabtext),
            ylab= ylabtext,
            xlab=col_name,
            col=c(color_3,color_2),
            horizontal=TRUE)
    plot(d_Yes,
         col=color_2,
         xlim=c(min_x,max_x),
         freq=FALSE,
         xlab=col_name,
         main= paste("Histogram of ", col_name ,"in terms of ", ylabtext, "= 1")
    )
    plot(d_No,
         col=color_3,
         xlim=c(min_x,max_x),
         freq=FALSE,
         xlab=col_name,
         main=paste("Histogram of ", col_name ,"in terms of ", ylabtext, "= 0")
    )
    par(mfrow=c(1,1))
  }
  
  
}
