install.packages("Rniftilib")
library("Rniftilib")

http://cran.r-project.org/web/packages/Rniftilib/Rniftilib.pdf
nifti.labelmap <- nifti.image.read("/paulsen/IPIG/predict_3T_MR/site-073/2334/43971/10_AUTO.NN3Tv20101112/43971_ACPC_Discrete_Class.nii.gz")
nifti.inputImg <- nifti.image.read("/paulsen/IPIG/predict_3T_MR/site-073/2334/43971/10_AUTO.NN3Tv20101112/43971_AVG_T1.nii.gz")


label_list<- levels( factor(nifti.labelmap[,,]) )
merged_image_mask <- merge( nifti.inputImg[,,], nifti.labelmap[,,])

hist(nifti.inputImg[,,][ which( nifti.labelmap[,,] == as.numeric(label_list[4])  ) ])



## T1 & T2 Merge and Compute Hist and Plot
  require(kd)

  firstImgLbs  <- T1[,,][which( labelMap[,,] == 10 ) ]
  secondImgLbs <- T2[,,][which( labelMap[,,] == 10 ) ]


  H.pi <- Hpi.diag( x=cbind(t1_label,t2_label), binned=TRUE)
  fhat <- kde( cbind(t1_label,t2_label), H=H.pi, binned=TRUE)
  plot(fhat, cont=seq(10,90,by=20), drawpoints=TRUE )
  plot(fhat, cont=seq(10,90,by=20),drawpoints=TRUE,pch=".",col=rgb(0,1,0,0.1) )


