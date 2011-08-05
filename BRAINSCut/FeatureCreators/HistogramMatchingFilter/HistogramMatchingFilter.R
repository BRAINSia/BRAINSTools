#!/usr/bin/r

# ------------------------------------------------------------------------------------ #
# * Plot three Histograms 
# 
#  USAGE::
#  $ R --slave --args  [data_file_name] [outputFigureFileName_Without_Extension] < plot_Bivariate.R
#
#
# ------------------------------------------------------------------------------------ #

myArg <- commandArgs();
dt_filename <- as.character( myArg[4] ); # first argument excluding R slave things
outputFilename<- as.character( myArg[5] );

print( paste( "dt_filename :: ", dt_filename ) );
# - Define Prefered Color
# - color name can be string such as 
#   
#   label_color <- rbind ( "blue", "red", ... etc.) if want to.
label_color <- rbind ( c( 1,0,0), 
                       c(0,1,0),  
                       c(0,0,1),  
                       c(1,1,0),
                       c(1,0,1),
                       c(0.5,0.5,0.5) ,
                       c(0.5,0.5,0.5) ,
                       c(0.5,0.5,0.5) ,
                       c(0.5,0.5,0.5) ,
                       c(0.5,0.5,0.5) ,
                       c(0.5,0.5,0.5) ,
                       c(0.5,0.5,0.5) ,
                       c(0.5,0.5,0.5) ,
                       c(0.5,0.5,0.5) 
                       );
label_color <- data.frame( label_color );

# give meaningful color name for better readability of the code

colnames( label_color ) <- c( "r","g","b"); 

#
# Command line argument Argument
#



# ------------------------------------------------------------------------------------ #
# Histogram 
# ------------------------------------------------------------------------------------ #

dt <- read.csv( dt_filename ,header=T);

colNames <- colnames(dt);
colnames(dt) <- c("bin","h1","h2","h3"); # give simple name

# ---------------------------------------------------------------
# Data Looks Like...
# ---------------------------------------------------------------
#  bin src_frequency ref_frequency adjusted_frequency
#  0         22506         31005              33084
#  1         22622         30173              34018
#  2         24968         29773              34257
#  3         29149         30749              35170
#  4         32983         31383              36961
#  5         33086         32824              38321
# ---------------------------------------------------------------
dev.set(10);
png( outputFilename, width=900, height=900 );
par(bg = "black", mfrow=c(3,1)) 

bar_colors = c("tan","violetred","gold")

# Decorating 



# Plot bar plot

barplot( dt$h1 , 
         col=bar_colors[1],
         border=bar_colors[1],cex.lab=2,
         )
x_label <- seq( dt$bin[1], dt$bin[nrow(dt)]*3, by = dt$bin[nrow(dt)]*3 /9 );
mtext ( text= x_label, side =1 , at = x_label ,  col = "grey" , line=1.5 );
axis( side =1, labels=T, tick=T, col="grey" , at = x_label, line=1);
axis( side =2, labels=T, tick=T, col="grey");
title(main=paste("Histogram Comparison: ", colNames[2]), col.main="white" ,
      cex.main=3);
barplot( dt$h2 , 
         col=bar_colors[2],
         border=bar_colors[2],cex.lab=2,
         )
x_label <- seq( dt$bin[1], dt$bin[nrow(dt)]*3, by = dt$bin[nrow(dt)]*3 /9 );
mtext ( text= x_label, side =1 , at = x_label ,  col = "grey" , line=1.5 );
axis( side =1, labels=T, tick=T, col="grey" , at = x_label, line=1);
axis( side =2, labels=T, tick=T, col="grey");
title(main=paste("Histogram Comparison: ", colNames[3]), col.main="white" ,
      cex.main=3);
barplot( dt$h3 , 
         col=bar_colors[3],
         border=bar_colors[3],cex.lab=2,
         )

x_label <- seq( dt$bin[1], dt$bin[nrow(dt)]*3, by = dt$bin[nrow(dt)]*3 /9 );
mtext ( text= x_label, side =1 , at = x_label ,  col = "grey" , line=1.5 );

axis( side =1, labels=T, tick=T, col="grey" , at = x_label, line=1);
axis( side =2, labels=T, tick=T, col="grey");
title(main=paste("Histogram Comparison: ", colNames[4]), col.main="white" ,
      cex.main=3);





temp <- dev.off(10);







