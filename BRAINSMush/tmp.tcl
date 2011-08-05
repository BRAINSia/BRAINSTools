OptimizeMushMixture $T1AlignedToACPCBfcImage $T2AlignedToACPCBfcImage $AnyAlignedToACPCBrainMask $MushAlignedToACPCBfcImage
      } 


      # set EdgeFormingRadius 1
      set Temp_MushThresholdBrainMask ${ScratchDirectory}/[file rootname [file rootname [file tail ${MushAlignedToACPCBfcBrainMask}]]]_Threshold.nii.gz
      if {[CheckOutputsNewer [list $MushAlignedToACPCBfcBrainMask $Temp_MushThresholdBrainMask] [list $MushAlignedToACPCBfcImage $AnyAlignedToACPCBrainMask] ] == false} {


        # Threshold mean-gm +|- C * stddev-gm;   
        # erode(1)-LargestFilledMask(5)-erode(3)-dilate(4)-LargestFilledMask(0)
        
        if {[CheckOutputsNewer [list $Temp_MushThresholdBrainMask] [list $MushAlignedToACPCBfcImage ${AnyAlignedToACPCBrainMask}] ] == false} {

          set mushImage [Brains::itk::LoadImage ${MushAlignedToACPCBfcImage} Float-Single]
          set Temp_MushEdgeImage ${ScratchDirectory}/[file rootname [file rootname [file tail ${MushAlignedToACPCBfcBrainMask}]]]_Edges.nii.gz
          set GreaterEdgeFormingRadius [expr ${LesserEdgeFormingRadius} + 1]
          puts "grayscale erode by ${GreaterEdgeFormingRadius}"
          set minImage [Brains::itk::ApplyStructuringElementToGrayscaleImage $mushImage Erode Ball ${GreaterEdgeFormingRadius}]
          puts "grayscale dilate by ${LesserEdgeFormingRadius}"
          set maxImage [Brains::itk::ApplyStructuringElementToGrayscaleImage $mushImage Dilate Ball ${LesserEdgeFormingRadius}]
          puts "subtract"
          set edgeImage [Brains::itk::ImageMath $maxImage $minImage Subtract Float-Single]
          ${minImage} Delete
          ${maxImage} Delete
          Brains::itk::SaveImage $edgeImage $Temp_MushEdgeImage

          set ROIMask [Brains::itk::LoadImage ${AnyAlignedToACPCBrainMask} Signed-16bit]
          set caution 7
          puts "erode by ${caution}"
          set CautiousROIMask [Brains::itk::ApplyStructuringElementToMaskImage ${ROIMask} Erode Ball ${caution}]
          ${ROIMask} Delete
          set ROI_Mush_Table [Brains::itk::measureLabelImageStatistics $CautiousROIMask $mushImage]
          set ROI_Edge_Table [Brains::itk::measureLabelImageStatistics $CautiousROIMask $edgeImage]
          ${CautiousROIMask} Delete

          set MushROIMean [Brains::Utils::ltraceSafe 0.0 $ROI_Mush_Table 0 1 2]
          set MushROIStdDev [Brains::Utils::ltraceSafe 0.0 $ROI_Mush_Table 0 4 2]          
          set lower [expr round( ${MushROIMean} * 0.5 ) ]
          set upper [expr round( ${MushROIMean} + ${MushROIStdDev} * 100.0 ) ]
          puts "\nMushROIMean ${MushROIMean} MushROIStdDev ${MushROIStdDev} lower ${lower} upper ${upper} "
          set threshToHeadMask [Brains::itk::BinaryThresholdImage ${mushImage} ${lower} ${upper} 0]
          ${mushImage} Delete
          set clipToHeadMask [Brains::itk::LargestRegionFilledMask ${threshToHeadMask} 1 1 10]
          ${threshToHeadMask} Delete

          set EdgeROIMean [Brains::Utils::ltraceSafe 0.0 $ROI_Edge_Table 0 1 2]
          set EdgeROIStdDev [Brains::Utils::ltraceSafe 0.0 $ROI_Edge_Table 0 4 2]          
          set lwindow [expr 0.1 * ${EdgeROIStdDev} - 0.5 ]
          set lower [expr round( 0.0 + $lwindow ) ]
          # round the lower bound down but the upper bound up.
          set uwindow [expr 1.0 * ${EdgeROIStdDev} + 0.5 ]
          set upper [expr round( ${EdgeROIMean} + $uwindow ) ]
          puts "\nEdgeROIMean ${EdgeROIMean} EdgeROIStdDev ${EdgeROIStdDev} lower ${lower} upper ${upper} "
          set threshToBrainMask [Brains::itk::BinaryThresholdImage ${edgeImage} ${lower} ${upper} 0]
          ${edgeImage} Delete

          set threshMask  [Brains::itk::ImageMath $threshToBrainMask $clipToHeadMask Minimum "Signed-16bit"]
          Brains::itk::SaveImage ${threshMask} ${Temp_MushThresholdBrainMask}
          
          ${threshMask} Delete
          ${threshToBrainMask} Delete
          ${clipToHeadMask} Delete
        }
        
