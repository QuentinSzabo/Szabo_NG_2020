# Szabo_NG_2020
Codes used to analyze microscopy images in Szabo et al., Nature Genetics 2020.

FISH_3D_analysis_QS.m

Uses as input files 2.2 x 2.2 µm (x and y) z-stacks of images (channel-separated, '.tif' format) surrounding FISH loci. Analyses the intermingling (overlap fraction and 3D distances between centroids) between segmented FISH probes (Otsu’s method) if 2 channels are selected for the analysis, or the structure (volume, principal axis length, sphericity, number of subdomains and subdomain volumes) of the segmented probe (Otsu’s method and watershed segmentation for subdomains) if one channel is selected for the analysis.

FISH_DAPI_2D_analysis_QS.m

Uses as input files 2.2 x 2.2 µm (x and y) z-stacks of images (channel-separated, '.tif' format) surrounding FISH loci or individual z-slices containing a single nucleus. Measures areas and extrapolated diameters of the channel-segmented objects using Otsu and watershed segmentations (if FISH analysis is selected, the analysis of areas and diameters are performed within individual z-slices randomly selected between the minimum and maximum z-coordinate of the 3D-segmented FISH object).

Codes were run on Windows 10 using MATLAB R2019b and its Image Processing Toolbox.
