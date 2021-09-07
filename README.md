# SOAX
A software for 2D/3D biopolymer network extraction and quantification.

About
-----
SOAX is an open source software tool to extract the centerlines, junctions and filament lengths of biopolymer networks in 2D and 3D images. It facilitates quantitative, reproducible and objective analysis of the image data. The underlying method of SOAX uses multiple Stretching Open Active Contours (SOACs) that are automatically initialized at image intensity ridges and then stretch along the centerlines of filaments in the network. SOACs can merge, stop at junctions, and reconfigure with others to allow smooth crossing at junctions of filaments.

SOAX provides 3D visualization for exploring image data and visually checking results against the image. Quantitative analysis functions based on extracted networks are also implemented in SOAX, including spatial distribution, orientation, and curvature of filamentous structures. SOAX also provides interactive manual editing to furthure improve the extraction results, which can be saved in a file for archiving or further analysis.

**\*Note\*** [TSOAX](https://github.com/tix209/TSOAX) (based on Qt 5) is an extension of SOAX to support tracking filaments and networks over multiple frames.

Reference
---------
T. Xu, D. Vavylonis, F. Tsai, G. H. Koenderink, W. Nie, E. Yusuf, I-Ju Lee, J.-Q. Wu, and X. Huang, "[SOAX: A Software for Quantification of 3D Biopolymer Networks](http://www.nature.com/srep/2015/150313/srep09081/full/srep09081.html)", In Scientific Reports, 5, 9081; DOI:10.1038/srep09081, 2015.

T. Xu, D. Vavylonis, X. Huang, "[3D Actin Network Centerline Extraction with Multiple Active Contours](http://www.sciencedirect.com/science/article/pii/S136184151300159X)", In Medical Image Analysis, 18(2):272-84, 2014.

Acknowledgement
---------------
This work has been supported by NIH grants R01GM098430, R01GM114201, and R35GM136372.
