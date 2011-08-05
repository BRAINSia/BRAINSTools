plot_RO_ROI <-function(ListFile,OutputDir,ROI){
#### Deal with list file
  listOfFile <- read.table(ListFile,as.is=TRUE)
  names(listOfFile)<-c("directory","Type")
  ##
  ## Read in the table
  ##
  count<-1;
  Nsample <- 10;
  subjects<-matrix(data="NA",nrow=Nsample,ncol=1)
  for( i in 1:nrow(listOfFile) ){
    print( listOfFile$Type[i] )
    if( listOfFile$Type[i] == "Apply" ) {
        subjects[count] <- unlist(strsplit( listOfFile$directory[i], "/"))[7]
        print( paste(OutputDir,"/",ROI,"RO_Vs_Threshold",subjects[count],".txt", sep="" ) )
        assign( paste( "tb",count,sep=""), 
               read.table(  paste(OutputDir,"/",ROI,"RO_Vs_Threshold",subjects[count],".txt", sep="" ),
                                  sep=":",header=T ) 
               )
        count<-count+1;
        print(count)
    }
  }

require( "fields" ) ## for xline
##
## write out the plot
##
Filename <- paste(OutputDir,"/",ROI,"_RO_plot.png",sep="");
png(Filename, width = 1000, height = 480  );
par(mfrow=c(1,2))
##
## plotting
##
my_color <- c( "firebrick4","darkorange","lemonchiffon4","seagreen2","darkgreen", "deepskyblue", "cadetblue4","mediumpurple4","pink","purple" );
y_max <- max( c ( tb1$RO , tb2$RO ,tb3$RO , tb4$RO ,tb5$RO , tb6$RO ,tb7$RO,tb8$RO,tb9$RO,tb10$RO) )
plot( tb1$Th, tb1$RO, col=my_color[1], type="l" ,ylim=c(0,1.0) ,
       ylab="Relative Overlap",
       xlab="Threshold Value"
       )
str_len <- length( strsplit( OutputDir, "/")[[1]] )
title( paste( unlist( strsplit( OutputDir, "/")[[1]][str_len] ) , ROI) )
#print( str_len )
#print( unlist( strsplit( OutputDir, "/")[[1]][str_len] ) )

points( tb2$Th, tb2$RO, col=my_color[2], type="l")
points( tb3$Th, tb3$RO, col=my_color[3], type="l")
points( tb4$Th, tb4$RO, col=my_color[4], type="l")
points( tb5$Th, tb5$RO, col=my_color[5], type="l")
points( tb6$Th, tb6$RO, col=my_color[6], type="l")
points( tb7$Th, tb7$RO, col=my_color[7], type="l")
points( tb8$Th, tb8$RO, col=my_color[8], type="l")
points( tb9$Th, tb9$RO, col=my_color[9], type="l")
points( tb10$Th, tb10$RO, col=my_color[10], type="l")

##
## plotting of average
##
#print(
#     rowMeans( 
#       matrix(
#         c( tb1$RO , tb2$RO ,tb3$RO , tb4$RO ,tb5$RO , tb6$RO ,tb7$RO ),
#         ncol=7),
#     ) 
#)
RO_means <- rowMeans( 
               matrix(
                 c( tb1$RO , tb2$RO ,tb3$RO , tb4$RO ,tb5$RO , 
                    tb6$RO ,tb7$RO , tb8$RO
                    , tb9$RO ,tb10$RO
                    ),
                 ncol=Nsample),
            ) 

points( tb1$Th,
        RO_means,
        type="l" ,
        lwd="3",
        lty=4,
        col="darkblue" )
###
### find max of average
###
max_avg <- max(RO_means)
max_avg_loc <- which.max(  RO_means )

xline( tb1$Th[max_avg_loc], col="mistyrose4", lty=2 ,lwd=5)

###
### Intraclass correlation (ICC)
###
require(psy)
max_avg_ref_vol <- rbind( tb1$Reference_volume[max_avg_loc] ,
                         tb2$Reference_volume[max_avg_loc] ,
                         tb3$Reference_volume[max_avg_loc] ,
                         tb4$Reference_volume[max_avg_loc] ,
                         tb5$Reference_volume[max_avg_loc] ,
                         tb6$Reference_volume[max_avg_loc] ,
                         tb7$Reference_volume[max_avg_loc] ,
                         tb8$Reference_volume[max_avg_loc] ,
                         tb9$Reference_volume[max_avg_loc] ,
                         tb10$Reference_volume[max_avg_loc]
                        )
max_avg_man_vol <- rbind( tb1$ANN_volume[max_avg_loc] ,
                         tb2$ANN_volume[max_avg_loc] ,
                         tb3$ANN_volume[max_avg_loc] ,
                         tb4$ANN_volume[max_avg_loc] ,
                         tb5$ANN_volume[max_avg_loc] ,
                         tb6$ANN_volume[max_avg_loc] ,
                         tb7$ANN_volume[max_avg_loc] ,
                         tb8$ANN_volume[max_avg_loc] ,
                         tb9$ANN_volume[max_avg_loc] ,
                         tb10$ANN_volume[max_avg_loc]
                        )

print( cbind( max_avg_ref_vol, max_avg_man_vol ) )
icc_result <- icc( cbind( max_avg_ref_vol, max_avg_man_vol ) )
###
### arrow.plot( max_avg_loc, 0.6 )
line_1 <-  paste("Maximum Average at " ,  tb1$Th[max_avg_loc] );
line_2 <- paste("value of ", round(max_avg,3) );
line_3 <- paste( "icc( A, C ): (",
                 round( icc_result$icc.consistency, 3 ) ,
                 " , ",
                 round( icc_result$icc.agreement, 3) ,
                 " )"
               )

print( "============================================" );
print( line_1);
print( line_2);
print( line_3);
print( "============================================" );
legend( min( 0.4 , tb1$Th[ max_avg_loc] ), 0.4 , line_1 , bty="n", col="red", cex=1.3)
legend( min( 0.4 , tb1$Th[ max_avg_loc]), 0.35, line_2 , bty="n", col="red", cex=1.3)


vol_min <- min( rbind( max_avg_ref_vol, max_avg_man_vol ) ) 
vol_max <- max( rbind( max_avg_ref_vol, max_avg_man_vol ) ) 
range <- ( vol_max - vol_min )

plot( max_avg_ref_vol ,  max_avg_man_vol , 
      xlim=c( vol_min - 0.1*range , vol_max + 0.1*range),
      ylim=c( vol_min - 0.1*range , vol_max + 0.1*range),
      col=my_color[1:Nsample], 
      pch=16,
      xlab="Reference Volume (1000cc)",
      ylab="ANN Volume (1000cc)")
title( ROI )
legend( "bottomright", 
        c(subjects[1:Nsample], "Average of RO"),
        col=c( my_color[1:Nsample], "darkblue" ), 
        pch=15,
        bty="n",
        cex=1.3)
legend( "topright" , line_3 , bty="n", col="red", cex=1.3)
abline( 0,1,col="red",lty=2,lwd=0.5)

dev.off()

}


plot_batch_RO <- function( OutputDir ){
  plot_RO_ROI( paste( OutputDir , "Accumben/Output",sep="/" ) ) ;
  plot_RO_ROI( paste( OutputDir , "Caudate/Output",sep="/" ) ) ;
  plot_RO_ROI( paste( OutputDir , "Globus/Output",sep="/" ) ) ;
  plot_RO_ROI( paste( OutputDir , "Hippocampus/Output",sep="/" ) ) ;
  plot_RO_ROI( paste( OutputDir , "Thalamus/Output",sep="/" ) ) ;
  plot_RO_ROI( paste( OutputDir , "Putamen/Output",sep="/" ) ) ;
}
