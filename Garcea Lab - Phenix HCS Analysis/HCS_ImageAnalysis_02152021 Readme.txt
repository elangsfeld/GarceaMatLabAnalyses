File: HCS_ImageAnalysis_02092021

Changelog:
--02.09.2021:
	-Fixed some new bugs that popped up in the previous version.
	-Added a field to the Results_Sorted variable: aSMAGradientMeansNucAdjusted.
		-This shows the aSMAGradientMean value for each [adjusted number] nucleus. Basically, it creates repeats for each nucleus in each cell object.
--02.09.2021:
	-Added functionality to InFocusImage function:
		-If there are <4 z-planes in the input image, then this function creates a MIP (rather than finding a focal plane). Additionally, this disables
		 the filter that removes images whose ZFocus value equate to the 1st or last z-plane.
		-If there are any *completely* empty images (as can happen in error occasionally), this function substitutes a nearby z-plane of the same
		 channel for the missing image. This avoid issues with later analysis, but could theoretically result in the incorrect focal plane being used.
	-Added a filter to HCS_ImageAnalysis that removes "hidden" files alongside real images. These files seem to be created by Mac machines when transferring
	 data to/from Windows machine. Their names start with "._" and end with ".tiff", which made them a problem for this script. These unwanted files are now
	 ignored.

Editables (Variables that can be adjusted to suit the data set):
--Folder: Where are the images located? IMPORTANT: Use format 'FILEPATH/' for Mac/Apple computers, and use format 'FILEPATH\' for Windows computers.
	  Also, the apostrophes and slash are important.
--FigShow: Do you want to show the Figure with segmentation overlay and demographic information? (1=yes, 0=no)
--FigAnnotate: Do you want to show the number designation for nuclei and cell objects in the Figure? (1=yes, 0=no)
--FigSave: Do you want to save the Figure that is generated during analysis? (1=yes, 0=no)
--InFocus_ImageSave: Do you want to save the in-focus image hyperstack (XYC)? (1=yes, 0=no)
--Channels: How many fluorescent channels are in the image? (This script is compatible with 2-4 channels.)
--CH_DAPI: Which channel corresponds to DAPI (or other nuclear dye) signal?
--CH_CellMask: Which channel corresponds to CellMask (or other cell membrane dye) signal?
--CM_LowSizeFilt: What is the minimum area (pixels^2) of cell objects that you want to include in the analysis? Don't touch this unless you are seeing a problem.
--TF_Analysis: Do you want to perform TF translocalization analysis? (1=yes, 0=no)
--CH_TF: Which channel corresponds to TF signal?
--aSMA_Analysis: Do you want to perform gradient-based aSMA analysis? (1=yes, 0=no)
--CH_aSMA: Which channel corresponds to aSMA signal?
--aSMA_ImageSave: Do you want to save an image of the aSMA channel and gradient intensity image? (1=yes, 0=no)
--SplitClusterResults: Do you want to split the results for each well into "clusters" (i.e., cell objects containing equal or more nuclei than the value of the
	  	       "MinClusterSize" variable) and "non-clusters"? (1=yes, 0=no)
--MinClusterSize: This dictates how "clusters" are identified for result sorting. Minimum value is inclusive.

Functions Called (in order they are called):
--InFocusImage:
	-Uses a gradient-based system to identify the z-plane that is most likely to be in focus. The rationale is that pixel-pixel gradients will be brighter
	 in z-planes with in-focus signal. If the "highest scoring" z-plane is the first or the last z-plane, then the image is thrown out of all subsequent
	 analyses, since the real focal plane is almost always outside of the image stack in these situations.
	-This uses the DAPI channel to find the focal plane.
	-Focal images are created for each channel, and they are max intensity projections of the focal plane and the two adjacent planes. This method tends to
	 more accurately capture focal signal in samples that aren't perfectly flat.
	-If there are <4 Z-planes in the image (as in TCPS data), a max intensity projection will be made for all available Z-planes (and the first/last ZFocus filter is disabled).
--CellMaskSegmentation:
	-Uses a marker-based watershedding segmentation method to identify the boundaries of Cell Mask signal and split them into "cell objects," which are
	 clusters of one or more cells.
	-Could technically work with any good cytoplasmic/cell membrane marker, but it was optimized for Cell Mask.
--NuclearSegmentation:
	-Uses a marker-based watershedding segmentation method to identify the boundaries of DAPI signal and split touching nuclei into individual objects,
	 when applicable. Also filters out objects that don't lie within "cell objects."
--CellularAnalysis:
	-This is the main feature-extraction function of the script. It extracts the intensity information for each channel, as well as nearest neighbor stats
	 for each nucleus. Importantly, nearest neighbor calculations are limited to neighbors WITHIN THE SAME CELL OBJECT.
--NucTranslocation:
	-This function will only be called if the "TF_Analysis" variable is set to 1.
	-Calculates the ratio of mean intensities between each nucleus and the "nucleus-free" cytoplasm in which the nucleus resides.
--aSMAActivation:
	-This function will only be called if the "aSMA_Analysis" variable is set to 1.
	-Uses a gradient-based system to quantify aSMA "activation" level in each cell object. In theory, cells with more aSMA fibers should yield higher
	 mean gradient values, but this theory hasn't been thoroughly vetted yet.
	**As far as I'm concerned (circa 12.23.2020), the rationale for this module still needs to be tested thoroughly. The activation level appears
	 to correlate really well with mean aSMA intensity, which I'm not sure it always should.
--WellSort:
	-This function sorts all the information from individual images in the "Results" variable into individual wells. In other words, it combines all
	 fields of view in each well into a single cell of results.
--ClusterSort:
	-This function will only be called if the "SplitClusterResults" variable is set to 1.
	-Uses the "MinClusterSize" variable value to split the results in each well into "in cluster" and "out cluster" subgroups.

Folders/Files Saved (some depend on options selected at top of script):
--FigureImages Folder, 'FILENAME Segmentation.tiff' files
	-This folder (nested under Analysis folder) and files will only be saved if the "FigSave" variable is set to 1.
	-An image of the figure showing the segmentation of DAPI and CellMask channels, as well as some demographic data.
		-If the "FigAnnotate" variable is set to 1, this saved image will also show the numbered designation of each nucleus and cell object.
--InFocusImages Folder, 'FILENAME In-Focus Image.tiff' files
	-This folder (nested under Analysis folder) and files will only be saved if the "InFocus_ImageSave" variable is set to 1.
	-A hyperstack of the in-focus image (as determined by InFocusImage function) will be saved. These are the images used for analysis.
--aSMAImages Folder, 'FILENAME aSMA Gradient Image.tiff' files
	-This folder (nested under Analysis folder) and files will only be saved if the "aSMA_Analysis" and "aSMA_ImageSave" variables are BOTH set to 1.
	-A hyperstack will be saved, with two channels. The first is the in-focus aSMA image. The second is the gradient intensity image. The values of the
	 second channel are used to calculate the mean gradient intensity in the Results files.
--AnalysisResults.mat
	-This file contains the "Results" variable. The analysis results are recorded for each image.
--AnalysisResultsSortedByWell.mat
	-This file contains the "Results_Sorted" variable. The analysis results are combined and sorted by well to make subsequent analysis easier.
--AnalysisResultsSortedByWellandOutsideClusters.mat
	-This file will only be saved if the "SplitClusterResults" variable is set to 1.
	-This file contains the "Results_Sorted_OutClusters" variable. The analysis results for each well (from "Results_Sorted") are separated by
	 how many nuclei were detected in that "cell object." The objects in this variable contain FEWER nuclei than the "MinClusterSize" variable value.
--AnalysisResultsSortedByWellandInsideClusters.mat
	-This file will only be saved if the "SplitClusterResults" variable is set to 1.
	-This file contains the "Results_Sorted_InClustesr" variable. The analysis results for each well (from "Results_Sorted") are separated by
	 how many nuclei were detected in that "cell object." The objects in this variable contain EQUAL/MORE nuclei than the "MinClusterSize" variable value.

Variables in each file:
--AnalysisResults.mat / "Results" variable:
	-FileName: The name of the file created by the FileSort script. (r = row, c = column, f = field).
	-Well: The first 6 characters of FileName, indicating the row and column of the image (i.e., the well).
	-FieldofView: The last 2 characters of FileName, indicating which field of view of the image.
	-ZFocus: The in-focus z-plane, as determined by a gradient-based scoring system in the InFocusImage function file.
	-NucNearestNeighborDistance: The distances between each nucleus and the nucleus nearest it (Euclidean by centroid-centroid distance). These distances are 
		WITHIN cell objects ONLY, so if there is only one nucleus detected in a cell object, there will be no nearest neighbor value. Duplicates are also
		removed, so a cell object with two nuclei will only report one nearest neighbor value.
	-TotalNuclei: The total number of nuclei detected in the image.
	-NucleiPerGroup: The number of nuclei detected in each cell object in the image.
	-AdjustedTotalNuclei: The total number of nuclei detected in the image, and adjusted using the median nuclear area value in the image. This theoretically
		helps correct for improperly segmented clumps of nuclei. As far as I know, this value isn't used for any additional analyses in the script unless
		it explicitly says otherwise.
	-AdjustedNucGroupSizes: The number of nuclei detected in each cell object in the image, but using the AdjustedTotalNuclei value.
	-MeanIntensities: The mean intensity values for each cell object in the image. Each column is a different channel (Col #1 = Channel 1, etc.).
	-SumIntensities: The sum (integrated) intensity values for each cell object in the image. Each column is a different channel (Col #1 = Channel 1, etc.).
	-CMNumber: The total number of cell objects detected in the image (based on Cell Mask (CM) signal segmentation).
	-CMAreas: The area (pixels) of each cell object in the image (based on Cell Mask (CM) signal segmentation).
	-CMTotalArea: The sum (integrated) area of all cell objects in the image (based on Cell Mask (CM) signal segmentation).
	-TFNucCytoRatios: The ratio of mean intensities between each nucleus and the cell object in which it resides. Will only be calculated if TF_Analysis = 1.
	-aSMAGradientMeans: The mean gradient value for each cell object in the image. Will only be calculated if aSMA_Analysis = 1.
--AnalysisResultsSortedByWell.mat / "Results_Sorted" variable:
	-Well: Same as in Results variable.
	-NumFields: The number of images that are analyzed for each well (images can be removed from analysis for various QC reasons).
	-TotalNumberofCMObjects: The sum of all "CMNumber" fields (from Results variable) for each well.
	-TotalNumberofNuclei: The sum of all "TotalNuclei" fields (from Results variable) for each well.
	-TotalAdjustedNumberofNuclei: The sum of all "AdjustedTotalNuclei" fields (from Results variable) for each well.
	-CMNumberObjects: The "CMNumber" field (from Results variable) for each image in each well.
	-CMTotalArea: The "CMAreas" field (from Results variable) for each image in each well.
	-NearestNucDistance: The "NucNearestNeighborDistance" field (from Results variable) for each image in each well.
	-TotalNuclei: The "TotalNuclei" field (from Results variable) for each image in each well.
	-NucleiPerGroup: The "NucleiPerGroup" fields (from Results variable) for each image in each well. The field values from each image are appended in order.
	-AdjustedTotalNuclei: The "AdjustedTotalNuclei" field (from Results variable) for each image in each well.
	-AdjustedNucGroupSizes: The "AdjustedNucGroupSizes" fields (from Results variable) for each image in each well. The field values from each image are appended in order.
	-MeanCh1: The 1st column of the "MeanIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-MeanCh2: The 2nd column of the "MeanIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-MeanCh3: The 3rd column of the "MeanIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-MeanCh4: The 4th column of the "MeanIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-SumCh1: The 1st column of the "SumIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-SumCh2: The 2nd column of the "SumIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-SumCh3: The 3rd column of the "SumIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-SumCh4: The 4th column of the "SumIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-TFNucCytoRatios: The "TFNucCytoRatios" field (from Results variable) for each image in each well. The field values from each image are appended in order.
		-This field will only be included if TF_Analysis = 1.
	-aSMAGradientMeans: The "aSMAGradientMeans" field (from Results variable) for each image in each well. The field values from each image are appended in order.
		-This field will only be included i aSMA_Analysis = 1.
	-aSMAGradientMeansNucAdjusted: The "aSMAGradientMeans" field (above) for EACH [ADJUSTED] NUCLEUS in each cell object.
--AnalysisResultsSortedByWellandOutsideClusters.mat / "Results_Sorted_OutClusters" variable:
	-Well: Same as in Results variable.
	-CMAreas: The "CMAreas" field (from Results variable) for the cell objects in each image in each well that were determined to be outside of clusters.
	-NumberCMObjects: The number of cell objects in each image in each well that were determined to be outside of clusters.
	-PercentofAllCMObjects: The percentage of all cell objects in all images in the well that were determined to be outside of clusters.
	-TotalCMArea: The sum of "CMAreas" field (from Results variable) for the cell objects in each image in each well that were determined to be outside of clusters.
	-PercentofAllCMArea: The percentage of all cell areas in all images in the well that were determined to be outside of clusters.
	-NucleiPerGroup: The "NucleiPerGroup" field (from Results variable) for cell objects in all images in the well that were determined to be outside of clusters.
	-AdjustedNucGroupSizes: The "AdjustedNucGroupSizes" field (from Results variable) for cell objects in all images in each well that were determined to be outside of clusters.
	-TotalNuclei: The sum of "TotalNuclei" field (from Results variable) for the cell objects in all images in the well that were determined to be outside of clusters.
	-PercentofAllNuclei: The percentage of all nuclei in all images in the well that were determined to be outside of clusters.
	-NearestNucDistanceFiltered: The "NucNearestNeighborDistance" field (from Results variable) for cell objects that were determined to be outside of clusters.
	-MeanCh1: The 1st column of the "MeanIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-MeanCh2: The 2nd column of the "MeanIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-MeanCh3: The 3rd column of the "MeanIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-MeanCh4: The 4th column of the "MeanIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-SumCh1: The 1st column of the "SumIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-SumCh2: The 2nd column of the "SumIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-SumCh3: The 3rd column of the "SumIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-SumCh4: The 4th column of the "SumIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-TFNucCytoRatios: The "TFNucCytoRatios" field (from Results variable) for nuclei with cell objects that were determined to be outside of clusters.
		-This field will only be included if TF_Analysis = 1.
	-aSMAGradientMeans: The "aSMAGradientMeans" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
		-This field will only be included if aSMA_Analysis = 1.
--AnalysisResultsSortedByWellandInsideClusters.mat / "Results_Sorted_InClusters" variable:
	-All fields are the same as in the "Results_Sorted_OutClusters" variable, but only cell objects containing a number of nuclei equal or greater than the MinClusterSize
		variable value are included for analysis.