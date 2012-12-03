inputDir="/hjohnson/HDNI/PREDICT_TRAINING/regina_ann/TrainingModels/BAW2012Dec//listFiles/"
maskDir="/hjohnson/HDNI/PREDICT_TRAINING/regina_ann/TrainingModels/BAW2012Dec/Manuals/Masks/"
DefDir="/hjohnson/HDNI/PREDICT_TRAINING/regina_ann/TrainingModels/BAW2012Dec/listFiles/Deformations/"
echo "siteID, subjectID, sessionID, imageList, roiList, deformationList"
while read site subject session
do
    roiDict=""
    for roi in accumben caudate putamen globus thalamus hippocampus
    do
        for side in l r
        do
            filename=(`ls $maskDir/${session}_${side}_$roi.nii.gz`)
            if [ "$roiDict" == "" ]; then
                roiDict="'${side}_$roi':'$filename'"
            else
              roiDict="$roiDict,'${side}_$roi':'$filename'"
            fi
        done
    done
    roiDict="{$roiDict}"

    ################################################################################
    
    imageDict=""
    for type in t1 
    do
        filename=(`ls $inputDir/${session}_${type}_average_BRAINSABC.nii.gz`)
        if [ "$imageDict" == "" ]; then
            imageDict="'$type':'$filename'"
        else
            imageDict="$imageDict, '$type':'$filename'"
        fi
    done
    imageDict="{$imageDict}"

    ################################################################################

    atlasToSubject=(`ls  $DefDir/${session}_AtlasToSubject_Composite.h5`)
    subjectToAtlas=(`ls  $DefDir/${session}_AtlasToSubject_InverseComposite.h5`)
    deformationDict="{'atlasToSubject':'$atlasToSubject','subjectToAtlas':'$subjectToAtlas'}"
    echo $site, $subject, $session, \"$imageDict\", \"$roiDict\", \"$deformationDict\"
done < subjectTest.list

while read site subject session
do
    imageDict=""
    for type in  t2
    do
        filename=(`ls $inputDir/${session}_${type}_average_BRAINSABC.nii.gz`)
        if [ "$imageDict" == "" ]; then
            imageDict="'$type':'$filename'"
        else
            imageDict="$imageDict, '$type':'$filename'"
        fi
    done
    imageDict="{$imageDict}"

    deformationDict="{'atlasToSubject':'$atlasToSubject','subjectToAtlas':'$subjectToAtlas'}"
    echo $site, $subject, $session, \"$imageDict\"
done < subjectTest.list 
