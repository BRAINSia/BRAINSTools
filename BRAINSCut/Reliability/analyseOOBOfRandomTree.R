drawLines<-function( dt, lineType, count,
                     error.min, error.max, NTree.min, NTree.max, myColors, numDepth) 
{
## start counting
  currCnt <- 1;
  currDepth <- numDepth[ currCnt];

  subDt <- subset( dt, dt$depth==currDepth );
  plot( subDt$NTree, subDt$error,
        type=lineType,
        pch=20, 
        xlim=c(NTree.min,NTree.max),
        ylim=c(error.min,error.max),
        ann = FALSE,
        frame.plot=F,
        col=myColors[ 1+count]);


## add plots
  for( currCnt in 2:length(numDepth) ) 
    {
    currDepth <- numDepth[ currCnt ];
    subDt <- subset( dt, dt$depth==currDepth );
    points( subDt$NTree, subDt$error,
            type=lineType, 
            pch=20,
            col=myColors[ currCnt+count ]);
    }
}
plotOOB <- function( filename, titleTxt, detailDepth)
{
## read data
  dt<-read.csv( filename, header=T);
## sort order
  dt<- dt[ order(dt$NTree, dt$depth ), ];

## color scheme
  numDepth = as.numeric(levels(factor( dt$depth)) );
  myColors=rainbow( length( numDepth) , alpha=0.3);

## range
  error.range <- max(dt$error) - min( dt$error ) ;
  error.min   <- min(dt$error ) - error.range * 0.05;
  error.max   <- max(dt$error ) + error.range * 0.05;

  NTree.range <- max(dt$NTree) - min( dt$NTree ) ;
  NTree.min   <- min(dt$NTree ) - NTree.range * 0.05;
  NTree.max   <- max(dt$NTree ) + NTree.range * 0.1;



pdf( paste(titleTxt,".pdf",sep=""));##, width=500, height=800 );

## 
  par( fig=c(0,1,0,1), new=T );
  drawLines( dt , "o", 1,
             error.min, error.max, NTree.min, NTree.max, myColors, numDepth);
  title( main=titleTxt , cex.main=1.6,
         xlab="number of trees", cex.lab=1.3,
         ylab="out of bag training error" );

##
  legend( "topright", 
          "depth",
          bty="n",
          cex=1.2
          );
  require(fields);

  image.plot( legend.only=T, 
              zlim=c( numDepth[1], numDepth[length(numDepth)]), 
              legend.shrink = 0.8,
              col=rainbow(length( numDepth)) );

## magnified plot
  magDT <- subset( dt, dt$depth > detailDepth );
  magNumDepth = as.numeric(levels(factor( magDT$depth)) );
  
  rect( 0,min( magDT$error ), 101,max( magDT$error ),
        border="red");

  par( fig=c(0.3,0.9,0.4,0.9), new=T );

  drawLines( magDT,  "l",  length(numDepth) - length(magNumDepth),
             min(magDT$error), max(magDT$error), NTree.min, NTree.max, myColors, magNumDepth);
  legend( "top",paste("detail ",detailDepth, "~") , bty="n" );


  dev.off();
}
