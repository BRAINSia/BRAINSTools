Hi Hans and Regina,

I’m just back from a summer vacation with my parents, so the summer is great here, and I hope it’s going pretty well for you too.

To help on your data I reviewed all super-resolution repositories, so I document all of my findings here for future reference.

First, I found this private repository:
https://github.uiowa.edu/SINAPSE/BRAINSSuperResolution
Above repository has an informative README file, and it seems that Matlab codes are fully customized for 4D DWI processing.

It is good to consider above repository as a general reference; however, I suggest to follow the following steps to work on 3D super-resolution (as the problem here).
Some minor changes may be needed to be applied to the Matlab codes that already exist in BRAINSTools repository (https://github.com/BRAINSia/BRAINSTools/tree/master/BRAINSSuperResolution).

First, the T1 image should be passed to “GenerateEdgeMap” tool to generate the edge map needed for super-resolution reconstruction algorithm.
You can find all the options to this program here: https://github.com/BRAINSia/BRAINSTools/blob/master/BRAINSSuperResolution/GenerateEdgeMap/GenerateEdgeMapImage.xml

Then, the output edge map should be passed to the super-resolution reconstruction algorithm (currently implemented in Matlab) along with the low-resolution FLAIR image.
Following directory was generated to show an example on DWI-b0 component (a 3D image):
https://github.com/BRAINSia/BRAINSTools/tree/master/BRAINSSuperResolution/SuperResolutionMatlabExample

Please look at the following M-file that takes the input low-resolution image and the edge map and generates the super-resolution reconstructed image:
https://github.com/BRAINSia/BRAINSTools/blob/master/BRAINSSuperResolution/SuperResolutionMatlabExample/SuperResolution.m

If the images are in NRRD formats, they can be read/written using nrrdLoadWithMetadata and nrrdSaveWithMetadata as shown in following example:
https://github.uiowa.edu/SINAPSE/BRAINSSuperResolution/blob/master/MATLAB_SCRIPTS/run_sr.m

Hope this helps, and please let me know if you have further questions.

Best,
Ali


From: Johnson, Hans J [mailto:hans-johnson@uiowa.edu]
Sent: Monday, July 3, 2017 7:01 PM
To: Ali Ghayoor <alii.ghayoor@gmail.com>
Cc: Eun Young (Regina) Kim <reginaeunyoungkim@gmail.com>
Subject: Help with super-resolution

Ali,

I hope your summer is going well.

Regina and I have a set of FLAIR images that are .5x.5x7mm  We have T1 images that are 1x1x1mm.   We would like to use your super resolution code to optimally resample the FLAIR images to the T1 image space.

Could you please help us find the tools and documentation to achieve this?

Thanks,
Hans

``` bash
#!/bin/bash

time python /scratch/SuperResolution/BRAINSSuperResolution/HCPWorkflows/HCPWF.py \
--inputDWIScan /Shared/johnsonhj/HDNI/ReferenceData/HCP_DATA/${SUBJID}/T1w/Diffusion/data.nii.gz \
--inputT1Scan /Shared/sinapse/CACHE/20160610_HCP_base_Results/HCP_DATA/${SUBJID}/HCP_${SUBJID}_01/TissueClassify/t1_average_BRAINSABC.nii.gz \
--inputT2Scan /Shared/sinapse/CACHE/20160610_HCP_base_Results/HCP_DATA/${SUBJID}/HCP_${SUBJID}_01/TissueClassify/t2_average_BRAINSABC.nii.gz \
--inputStandardLabels /Shared/sinapse/CACHE/20160610_HCP_base_Results/HCP_DATA/${SUBJID}/HCP_${SUBJID}_01/JointFusion/JointFusion_HDAtlas20_2015_fs_standard_label.nii.gz \
--inputLobeLabels /Shared/sinapse/CACHE/20160610_HCP_base_Results/HCP_DATA/${SUBJID}/HCP_${SUBJID}_01/JointFusion/JointFusion_HDAtlas20_2015_lobe_label.nii.gz \
--program_paths /scratch/NAMICExternalProjects/release-20160523/bin \
--python_aux_paths '/scratch/SuperResolution/BRAINSSuperResolution/HCPWorkflows:/scratch/wmql/tract_querier/tract_querier/' \
--labelsConfigFile /scratch/BS/BRAINSTools/AutoWorkup/DWIProcessingWorkflows/FS_Extended_Labels_Config.csv \
--workflowCacheDir /Shared/johnsonhj/HDNI/20160804_HCP_Processing_Pipeline/Results/${SUBJID} \
--resultDir /Shared/johnsonhj/HDNI/20160804_HCP_Processing_Pipeline/Results/${SUBJID}
_${SUBJID}

```
