# Identification of differential gene expression responses in pulmonary cells from bulk RNA-Seq using deconvolution
University of Minnesota Supercomputing Institute

RI Analyst: Natalia Calixto Mancipe

This project aims to determine the cell-specific differential gene expression and cell number differences in pulmonary cells extracted from rats’ lungs from three experimental groups: control, pulmonary hypertension (induced by monocrotaline injection) treated with vehicle and monocrotaline treated with ferrostatin-1, a ferroptosis inhibitor.  RNAseq was extracted from bulk tissue from 4 Male Sprague-Dawley rats per treatment (12 libraries total).

Illumina reads were aligned to the reference genome and used to estimate gene counts with the CHURP-PURR pipeline. The resulting count matrixes are stored in the bulk_data folder. 

The reference single cell data used for deconvolution was taken from [Hong et al. 2020](https://www.atsjournals.org/doi/full/10.1164/rccm.202006-2169OC?role=tab). Only the control and monocrotaline data was used.

These matrixes were further analysed with the scripts found in the 'code' folder as follows:

1. scPrep.r : get reference single cell data and cell type markers using Seurat's algorithm and the author's cell type labels.
2. sigMatrix.MAST.r : fix bug in DWLS buildSignatureMatrixMAST function and run it separately for control and MCT.
3. solveDWLS.r : deconvolve each experimental treatment using the signature matrixes from the previous step.
4. TOAST.r : define cell types of interest and test for differential gene expression. Perform GSEA using the differential expression tests results.
