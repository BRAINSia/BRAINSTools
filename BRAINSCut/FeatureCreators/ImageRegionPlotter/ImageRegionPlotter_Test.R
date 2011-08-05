# ------------------------------------------------------------------------------------ #
# * Plot Countour Plot ( Bivariate Density Plot) 
#  USAGE::
#  $ R --slave --args  filename1 filename2 outputFigurename < plot_Bivariate.R
# ------------------------------------------------------------------------------------ #

#
# Argument
#

myArg <- commandArgs();
imageFilename1<- as.character( myArg[4] );
imageFilename2<- as.character( myArg[5] );
labelMapFilename<- as.character( myArg[6] );
outputFilename <- as.character( myArg[7] );


#
# Read in Images
#

#
# * R Nifti Library
# * See http://cran.r-project.org/web/packages/Rniftilib/Rniftilib.pdf
#
require("Rniftilib");

#
# * Read in Images
#
Img1 <- nifti.image.read( imageFilename1 );
Img2 <- nifti.image.read( imageFilename2 );
LabelMap <- nifti.image.read( labelMapFilename );

#
# * Get Labels in the LabelMap
#
labelList <- levels( factor(LabelMap[,,]) ); 

#
# * Compute Histogram and Plotting
#

require(ks);
label_color <- cbind( rgb( 1,0,0,0.01),  
                      rgb( 0,1,0,0.01),
                      rgb( 0,0,1,0.01),
                      rgb( 0,1,1,0.01),
                      rgb( 1,0,1,0.01),
                      rgb( 1,1,0,0.01) );

png( outputFilename );
# 
# First Basic Plot [Skipping 0]
#
print( paste(" Plot Label :: " , labelList[4] ));
ImgLab1 <- Img1[,,][ which( LabelMap[,,] == labelList[4] ) ] ;
ImgLab2 <- Img2[,,][ which( LabelMap[,,] == labelList[4] ) ] ;
H.pi <- Hpi.diag( x=cbind( ImgLab1, ImgLab2 ), binned=TRUE);
fhat <- kde( cbind( ImgLab1, ImgLab2 ),
             H=H.pi, 
             binned=TRUE);

plot( fhat, 
      cont=seq(10,90,by=20) ,col="red");
     #drawpoints=TRUE,
#      pch=".",
#      col=label_color[1]    );
#
# * Plot Rest of Labels
#
#for label_index in 3:length(labelList){
#  print( paste(" Plot Label :: " , labelList[label_index] );
#
#  ImgLab1 <- Img1[,,][ which( LabelMap[,,] == labelList[ label_index ] ) ] ;
#  ImgLab2 <- Img2[,,][ which( LabelMap[,,] == labelList[ label_index ] ) ] ;
#  H.pi <- Hpi.diag( x=cbind( ImgLab1, ImgLab2 ), binned=TRUE);
#  fhat <- kde( cbind( ImgLab1, ImgLab2 ),
#               H=H.pi, 
#               binned=TRUE);
#  plot( fhat, 
#        cont=seq(10,90,by=20),
#        #drawpoints=TRUE,
#        pch=".",
#        col=label_color[label_index] ,
#        add=TRUE
#        );

#}
dev.off();


