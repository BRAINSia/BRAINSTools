# ------------------------------------------------------------------------------------ #
# * Plot Bivariate Plot ( Bivariate Density Plot)
#
#  USAGE::
#  $ R --slave --args  \
#                       [data_file_name] \
#                       [label_mapping file name ] \
#                       [outputFigureFileName_Without_Extension] \
#                       < ImageRegionPlotter2.R
#
#  Input Description
#
#  - [data_file_name]
#
#    :: a text data file name having label, index_for_image1, index_for_image2,
#       and frequencies. First three line includes a file name with full path,
#       and real data starts from 4th line
#
#    :: data file ex) --> Comma separated file!
#
#       [Image1]: /BSITKBRAINSABC/CORRECTED_INDEX_0_LEVEL_4.nii.gz
#       [Image2]: /BSITKBRAINSABC/CORRECTED_INDEX_2_LEVEL_4.nii.gz
#       [LabelMap]: /2334_43971_3DSPGR_0_COR_ACPC_DirtyLabels_BRAINSABC.nii.gz
#       label, CORRECTED_INDEX_0_LEVEL_4.nii.gz, CORRECTED_INDEX_2_LEVEL_4.nii.gz, FrequenceOfQuantile
#       1,0,0,0
#       1,0,1,0
#       1,0,2,0
#       1,0,3,0
#       1,0,4,0
#
#  - [label_mapping file name]
#    :: label mapping file
#       label, name , color
#       0, AIR, black
#       3, BGM, grey
#       4, CSF, yellow
#       2, GM, blue
#       6, NOTCSF, grey
#       7, NOTGM , grey
#       9, NOTVB, grey
#       8, NOTWM, grey
#       5, VB, grey
#       1, WM, red
#
#
#
#  - [outputFigureFileName_Without_Extension]
#    :: output figure file name with out extension
#    :: for now we are creating png file image look for 'png' in the code
#
#
#
# ------------------------------------------------------------------------------------ #

# - Constant Used in the Code
# - Label 'String' Names to be displayed
# - Will automatically mapped from the first label
# - increasing order.
# - ex) 1 -> "WM, 2->"GM"
# - no effect on results, but for displaying. look for 'legend'


#
# Command line argument Argument
#

myArg <- commandArgs();
dt_filename <- as.character( myArg[4] ); # first argument excluding R slave things
lb_filename <- as.character( myArg[5] ); # second argmumet for label mapping table
outputFilename<- as.character( myArg[6] );
mode<- as.character( myArg[7] ); # mode for s, bw, bs

# ------------------------------------------------------------------------------------ #
# - min/max line for each label
# ------------------------------------------------------------------------------------ #
DrawBoundingBoxBasedOnStd<- function( dt, nsd , index_color_rgb ){
  # percent <- 0.00001;
  img1.totalFreq <- sum( dt$freq );
  img1.mean <- sum( dt$freq * dt$img1 ) / img1.totalFreq ;
  img1.std  <- sqrt (sum( ((dt$img1 - img1.mean)^2)* dt$freq) / img1.totalFreq ) ;

  print( paste(" image1 mean::", img1.mean,  ", var:: ", img1.std ) );

  img2.totalFreq <- sum( dt$freq );
  img2.mean <- sum( dt$freq * dt$img2 ) / img1.totalFreq ;
  img2.std  <- sqrt (sum( ((dt$img2 - img2.mean)^2) *dt$freq ) / img2.totalFreq ) ;

  print( paste(" image2 mean::", img2.mean,  ", var:: ", img2.std ) );


  img1.min = img1.mean - nsd * img1.std; if( img1.min < 0 ) { img1.min = 0; }
  img1.max = img1.mean + nsd * img1.std; if( img1.max > 100 ) { img1.max = 100; }
  img2.min = img2.mean - nsd * img2.std; if( img2.min < 0 ) { img2.min = 0; }
  img2.max = img2.mean + nsd * img2.std; if( img2.max > 100 ) { img2.max = 100; }

  rect( img1.min, img2.min, img1.max, img2.max ,
        border = rgb( index_color_rgb[1],
                      index_color_rgb[2],
                      index_color_rgb[3] )
       );

}

# ------------------------------------------------------------------------------------ #
# - min/mix line for each label
# ------------------------------------------------------------------------------------ #
DrawBoundingBoxBasedOnPercent <- function( label_dt, percent, index_color_rgb ){
  # percent <- 0.00001;
  threshold <- sum( label_dt$freq ) * percent;

  sub_label_dt<- subset( label_dt, label_dt$freq > threshold  ) ;

  print(head( sub_label_dt));

  img1_min = min( sub_label_dt$img1 ); print( img1_min );
  img1_max = max( sub_label_dt$img1 ); print( img1_max );
  img2_min = min( sub_label_dt$img2 ); print( img2_min );
  img2_max = max( sub_label_dt$img2 ); print( img2_max );

  rect( img1_min, img2_min, img1_max, img2_max ,
       border = rgb( index_color_rgb[1],
                   index_color_rgb[2],
                   index_color_rgb[3] )
       );

}
# ------------------------------------------------------------------------------------ #
# Contour plot Function: Put White Lines to the plot.
# ------------------------------------------------------------------------------------ #

ContourPlot <- function( data , colorname) {

  sumFreq <- sum( data$freq );

  contourFreq <- matrix( data$freq / sumFreq,
                         nrow = 101,
                         ncol = 101,
                         byrow = TRUE);

  contour_levels <- pretty( contourFreq, 3 );
  print(contour_levels);
  contour( seq( 0, 100, by=1 ),
         seq( 0, 100, by=1),
         contourFreq,
         col = rgb(colorname[1], colorname[2],colorname[3],),
         add = TRUE,
         levels = contour_levels ,
#nlevel=5,
#col = rgb(colorname),
       );
  print("********************************");

}

# ------------------------------------------------------------------------------------ #
# Plot Joint Histogram
# ------------------------------------------------------------------------------------ #

PlotJointHistogram <- function( dt_filename,
                                lb_filename,
                                outputFilename ,
                                numberOfLineForSkip,
                                mode='s')
{
   #
   # Read in Data File
   #
   # -- Ex of Data File )------------------
   # label img1 img2 freq
   #  1     1         0         0        20
   #  2     1         0         1        84
   #  3     1         0         2       275
   #  4     1         0         3       482
   # --------------------------------------

   #
   # - Read comma separated file, skipping first 3 lines as 4th line column names
   print("read dt");
   dt <- read.csv( dt_filename, skip = numberOfLineForSkip  , header=T);
   print( head( dt ) )

   XImageFilename <- colnames( dt ) [2];
   YImageFilename <- colnames( dt ) [3];
   if( mode=='s') #individual subject
   {
      #
      # - retrieve some information from column name and then rename columes
      # -  for coding convenience.

      FreqeuncyName <- colnames( dt ) [4];
   }
   else if( mode=='bw') # batch mode weighted sum
   {
     dt <- subset( dt, select=c("label","bin1","bin2", "normalizedFreqSum") ) ;
     FreqeuncyName <- " weighted sum of frequency"
   }
   else if( mode=='bs') # batch mode simple sum
   {
     dt <- subset( dt, select=c("label","bin1","bin2", "freqSum") ) ;
     FreqeuncyName <- " sum of frequency"
   }

      # - rename columes

   colnames(dt)<-c("label","img1", "img2", "freq");

   #
   # - get existing labels in the data set
   #
   listOfLevel <- levels(factor(dt$label)) ;

   #
   # read label table

   print("read lb");
   lb_table <- read.csv(  lb_filename, header=T ) ;

   lb_used_table <-matrix(data=NA,ncol=3) ;
   colnames(lb_used_table) <- c("label","name", "color");

   for( i in seq( 1, length( listOfLevel) ) )
   {
      lb_used_table<-rbind( lb_used_table,
                            subset( lb_table, lb_table$label == listOfLevel[ i] )) ;
   }
   lb_used_table <- lb_used_table[-1,]; # delete first row with NAs in it

   lb_colors <- matrix( col2rgb( as.character(lb_used_table$color)), ncol=3 ,byrow=T);
   lb_colors <- lb_colors/255;


   #
   # - plot each label
   #

   #


   #
   # - Get Subset for the first label we have
   Nlb <- length( listOfLevel );

   #label_dt<- subset( dt, dt$label == as.numeric(listOfLevel[ Nlb ]) ) ;
   label_dt <- subset( dt, ( dt$label == as.numeric(listOfLevel[ Nlb]) &
                               dt$freq != 0 )) ;

   #
   # - Plot the first label
   frequency_devider <-  max ( label_dt$freq );

   #
   # - write figure to the device
   #

   png( paste(outputFilename, "png", sep="."),
        width = 800, height = 800,
       ) ; #Output Figure


   print( paste( "Processing Label : ", listOfLevel[ Nlb ]  ) );

   #
   # First plot
   #

   # make black background
   par(bg = "white")
   par(mar=c(5, 5, 2, 5) + 0.1);

   if( mode == 's' ){
       # three  figures.

       layout( matrix( c(3 ,1,  2, 0 ) , 2,2,byrow=T), width=c(4,1) , height=c(4,1));

       #
       # plot mapping
       #

       print( paste( dt_filename, "map1.txt", sep="") );
       #
       #figure locations
       #
       #mid.x <- 0.8;  mid.y<-0.8;
       #
       # figure 1:: Volume 2's Intensity vs. Quantiles
       #
       map1 <- read.csv( paste( dt_filename, "map1.txt", sep=""), header=T );
       map1$quantile <- map1$quantile * 100;
       #par( fig = c( mid.x, 1, 0, mid.y ) )
       plot( map1$intensity, map1$quantile  , type="l" ,col="grey",
             axes=F,
             xlab="",
             ylab="",
             xlim = c( 0.0, max( map1$intensity ) * 1.001 ),
             ylim = c( -0.5, 105 )
           );
       # - ticks at 0,0.25,0.5, 0.75 and 1.00 Quantiles
       q_ticks <- c( 0, 25, 50, 75, 100 );
       i_ticks <- c( map1$intensity[which( map1$quantile == 0) ],
                     map1$intensity[which( map1$quantile == 25)],
                     map1$intensity[which( map1$quantile == 50) ],
                     map1$intensity[which( map1$quantile == 75)],
                     map1$intensity[which( map1$quantile == 100) ] );

       axis( 4 , at=q_ticks, col="grey", col.ticks="grey", labels=F);
       mtext(q_ticks, at=q_ticks,col="grey", side=4, line=1.5);
       mtext("Volume 2 Quantile", side=4, line=3.5 , col="grey");

       axis( 1 , col="grey", col.ticks="grey");
       mtext("Volume 2 Intensity", side=1, line=2.5 , col="grey");

       #
       # - pointsat quantiles with values are specified
       position_points <- ( rbind( i_ticks[2:5], q_ticks[2:5] )  ) ;
       print( position_points )
       points( position_points,  col="grey" );
       mtext( c( "I:",i_ticks[2:5]),  at=q_ticks,col="grey", side=4, line=-3.5);


       #
       # figure 2
       map2 <- read.csv( paste( dt_filename, "map2.txt", sep=""), header=T );
       map2$quantile <- map2$quantile * 100;

       #par( fig = c( mid.x, 1, 0, mid.y) );
       plot( map2$quantile, map2$intensity, type="l" ,
             bty='n',
             col="grey",
             axes=F,
             xlab="",
             ylab="",
             xlim = c( -0.5, 105 ),
             ylim =  c( 0.0, max( map2$intensity ) ),
       #     new=T
             );
       # - ticks at 0,0.25,0.5, 0.75 and 1.00 Quantiles
       i_ticks <- c( map2$intensity[which( map1$quantile == 0) ],
                     map2$intensity[which( map1$quantile == 25)],
                     map2$intensity[which( map1$quantile == 50) ],
                     map2$intensity[which( map1$quantile == 75)],
                     map2$intensity[which( map1$quantile == 100) ] );

       axis( 1 , at=q_ticks, col="grey", col.ticks="grey", labels=F);
       mtext(q_ticks, at=q_ticks,col="grey", side=1, line=1.5);
       mtext("Volume 2 Quantile", side=1, line=3.5 , col="grey");

       axis( 2 , col="grey", col.ticks="grey", );
       mtext("Volume 2 Intensity", side=2, line=2.5 , col="grey");

       #
       # - pointsat quantiles with values are specified
       print( t( rbind( i_ticks[2:5], q_ticks[2:5] )  ) );
       position_points <- t( rbind( q_ticks[2:5], i_ticks[2:5] )  ) ;
       points( position_points,  col="grey" );
       mtext( c( "I:",i_ticks[2:5]),  at=q_ticks,col="grey", side=1, line=-3.5);

   }

   #
   # figure 3
   #
   # find proper color from label_color table

   print(lb_used_table);

   index_color <- which( lb_used_table$label == listOfLevel[ Nlb ] ) ;

   print( paste( "index_color::" , index_color ) );

   index_color_rgb <- lb_colors[ index_color, ];

   print( paste( "index_color_rgb::" ,index_color_rgb) );

   #par( fig = c( 0, mid.x, 0, mid.y) )
   plot( label_dt$img1,
         label_dt$img2,
         col = rgb( index_color_rgb[1],index_color_rgb[2],index_color_rgb[3] ,
               round( label_dt$freq/ frequency_devider + 0.002, 2)  ),
         pch=16,
         axes = FALSE, # no axis at this point
         ann = FALSE,  # so that we can give color for the axis
         xlim=c( -0.05, 105),
         ylim=c( -0.05, 105),
         );
   #
   # - contour plot
   #
   #ContourPlot ( label_dt , index_color_rgb);

   #
   # - min/mix line for each label
   #
   percent <- 0.00001;
   numberOfStd <- 2.58;
   #DrawBoundingBoxBasedOnPercent( label_dt, percent, index_color_rgb);
   DrawBoundingBoxBasedOnStd( label_dt, numberOfStd , index_color_rgb )

   #
   # - legend
   #

   legend("topright",
          lb_used_table$name,
          pch = 15,
          col= lb_used_table$color,
          bty="n" , # no border for the legend
          text.col = lb_used_table$color
          );
   #
   # - Title and Axis Specification
   #

   title( main= FreqeuncyName ,
          col.main= "gray",
          xlab= paste( "Volume 1 :", XImageFilename ) ,
          ylab= paste( "Volume 2 :", YImageFilename ) ,
          col.lab = "gray"
        );
   ticks <- seq( 0,100,by=10 );
   mtext ( text=ticks,side =1 , at = ticks ,
           col = "grey" );
   axis( 1 , # x axis
         col = "gray",
         labels = FALSE,
         at = ticks
       );

   mtext ( text=ticks,side = 2 , at = ticks ,
           col = "grey" );
   axis( 2 , # y axis
         col = "gray",
         labels = FALSE,
         at = ticks
       );

   for( lb_index in seq( Nlb, 1, by=-1 ) ){

     print( paste( "Processing Label : ", listOfLevel[ lb_index]  ) );

     label_dt <- subset( dt, ( dt$label == as.numeric(listOfLevel[ lb_index ]) &
                               dt$freq != 0 )) ;
     frequency_devider <- max( label_dt$freq ) ;

     #
     # following plot added to the above plot
     #
   index_color <- which( lb_used_table$label == listOfLevel[ lb_index ] );
   index_color_rgb <- lb_colors[ index_color ,];

     points ( label_dt$img1,
              label_dt$img2,
              col = rgb( index_color_rgb[1],
                         index_color_rgb[2],
                         index_color_rgb[3] ,
                    round(  label_dt$freq/ frequency_devider +0.002  , 2)  ),
              pch=16);

   #ContourPlot ( label_dt , index_color_rgb )
     #
     # - min/mix line for each label
     #
   #DrawBoundingBoxBasedOnPercent( label_dt, percent, index_color_rgb);
     DrawBoundingBoxBasedOnStd( label_dt, numberOfStd  , index_color_rgb );

   }



   #
   # TODO Add small plots
   #

   dev.off();
}

#
# actual calling
#
print( paste( "mode:::: ", mode) );
if( mode == 's' )
{
  PlotJointHistogram( dt_filename, lb_filename, outputFilename , 3, mode );
} else
{
  PlotJointHistogram( dt_filename, lb_filename, outputFilename , 0, mode );
}


