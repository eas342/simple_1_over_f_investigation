This README describes some of the data within this OTIS Dark Frame Experiment.

The 108 frame dark ramp is "sliced and diced" into 54 subtractive pairs.

These 54 pairs are then used to pretend we have a dark frame time series.
This experiment allows us to test out the contributions of read noise to exoplanet time series.

## Original Raw Data:
The original raw data for the ALONG (A5) detector is here
/surtrdata1/External/OTIS_unzipped/NRCNRCALONG-DARK-7235074213_1_1_924_JW1_JLAB40_20170823T074323.606_20170823T080253.912/NRCNRCALONG-DARK-72350742131_1_485_SE_2017-08-23T16h49m51.fits

For the A3 detector, it is here:
/fenrirdata1/External/OTIS_unzipped/NRCNRCA3-DARK-7255203548_1_1_3697_JW1_JLAB40_20170912T203659.807_20170912T205630.113/NRCNRCA3-DARK-72552035481_1_483_SE_2017-09-12T23h40m37.fits

## Directories:
pairwise_sub - the pairwise subtracted 54 pairs from the raw ALONG (A5) data. No reference pixel or other NCDHAS pipeline processing is applied.
pairwise_sub_A3 - the pairwise subtracted 54 frames from the raw A3 data. No reference pixel or other NCDHAS pipeline processing is applied.
pairwise_sub_red - the pairwise subtracted 54 pairs from the NCDHAS-processed ALONG (A5) data with default parameters
pairwise_sub_red_no_side_refcor - the pairwise subtracted 54 pairs from NCDHAS-processed ALONG (A5) data but with side reference pixels subtraction fully turned off. See directory "proc2" as well.
pairwise_sub_red_rowSub - the pairwise subtracted 54 pairs from reduced ALONG (A5) data. The pairs were made from a red file where each group had a row-by-row median subtraction after reduction by ncdhas.
pairwise_sub_red_rowColSub - the pairwise subtracted 54 pairs from the reduced ALONG (A5) data. The pairs were made from a red file where each group had a row-by-row and column-by-column median subtraction after reduction by ncdhas.
pairwise_sub_red_smoothedRowKernel - pairwise subtracted 54 pairs from the reduced ALONG (A5) data. The pairs were made from a red file where each group had a subtraction from a smoothed image where the smoothing was a 40 pix moving average kernel along the row direction.
grp_split_red_file - directory where the processed ALONG red file is split into one file per group
grp_split_red_A3 - directory where the processed A3 red file is spit into one file per group
grp_split_red_smoothed_RowKernel_file - directory where a ALONG red file is spit into one file per group. Here the processed red file is one where each group is subtracted by smoothing kernel along the rows.

## Old pairwise files:
old_1_over_f_after_pairwise_sub/pairwise_sub_red_additional_sub - the pairwise subtracted 54 pairs from NCDHAS-processed ALONG (A5) data with an extra median subtraction along each row only. The row-by-row subtraction is done on the slope image instead of the red image.
old_1_over_f_after_pairwise_sub/pairwise_sub_red_rowcol_sub - the pairwise subtracted 54 pairs from NCDHAS-processed ALONG (A5) data with an extra median subtraction along each row and column. The row-by-row and column-by-column subtraction is done on the slope image instead of the red image.

## Proc Red files:
proc - the NCDHAS output from processing the raw ALONG (A5) and raw A3 dark exposures with default NCDHAS parameters.
proc2 - the NCDHAS output of ALONG (A5) data from processing the raw A5 dark exposure with side reference pixel subtraction fully turned off. (see "run_ncdhas.sh")
proc_red_additional_rowSub - the same as proc but with the median of each row in each read subtracted from that row. This should be better than subtracting after the fact
proc_red_additional_rowcol_sub - the same as proc but with the median of each row and median of each column subtracted for each group/read.
proc_red_smoothedRowKernel - the same as proc, but each group is subtracted by a smoothing kernel applied along the rows
proc_red_refAmpFlip0 - uses one reference amplifier (the first) that is copied and flipped over the image to subtract other amplifiers.
proc_red_pcaEach - uses PCA with 10 principal components on each image. The first 10 components are mulitplied by their eigenvectors to create an model for each group. Each image will hav its own principal component vectors.

## Scripts:
sub_pairwise.py - takes a 108 frame dark frame and creates read pairs
additional_1_over_f_correction.py - applies median subtractions to reduce 1/f noise
split_groups.py - split groups into separate files



