# restingstate-reliability-spinalcord

This repository is associated with the following [manuscript](https://www.sciencedirect.com/science/article/pii/S1053811923003038?via%3Dihub#sec0036) and the corresponding [dataset](https://openneuro.org/datasets/ds004386).

This repository contains the preprocessing and analysis code needed to assess the reliability of spinal cord resting-state networks as presented in above-mentioned manuscript. Note that the functions are organized according to the organization of the manuscript.  

## Required software
- [MATLAB](https://www.mathworks.com/products/matlab.html), version 2018a or higher

- [Spinal Cord Toolbox(SCT)](https://spinalcordtoolbox.com/en/latest/), version 4.2.2 or higher

- [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki), version 5.0 or higher

- [AFNI](https://afni.nimh.nih.gov/),  version 22.1.09n

### Github repos
- https://github.com/NYU-DiffusionMRI/mppca_denoise

- https://github.com/andersonwinkler/PALM

- https://github.com/ekmolloy/fmri_test-retest

### Matlab repo
[distributionPlot](https://github.com/eippertlab/zshim-spinalcord/tree/main/ZShim_Results/Step2_CalculateResults/Helper_Code/Figure_Code/distributionPlot/distributionPlot) :
Jonas (2021). Violin Plots for plotting multiple distributions (distributionPlot.m) (https://www.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m), MATLAB Central File Exchange.


## Preprocessing & Analysis
For preprocessing and calculation of results, [SCT version 4.2.2](https://github.com/spinalcordtoolbox/spinalcordtoolbox/releases/tag/4.2.2), FSL version 6.0.3, and MATLAB version 2021a were used.

## References
- Kaptan, M., Horn, U., Vannesjo, S. J., Mildner, T., Weiskopf, N., Finsterbusch, J., Brooks, J. C. W., & Eippert, F. (2023). Reliability of resting-state functional connectivity in the human spinal cord: Assessing the impact of distinct noise sources. NeuroImage, 120152. https://doi.org/10.1016/j.neuroimage.2023.120152
- https://openneuro.org/datasets/ds004386
