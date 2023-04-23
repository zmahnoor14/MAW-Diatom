# MAW-Diatom

Metabolome Annotation Workflow (MAW) was used to perform metaboli profiling of Skeletonema marinoi using the spectral data from Liquid Chromtohtaphy Tandem Mass Spectrometry (LCMS2). The data of the RAW and mzML files is submitted to MetaboLights with the accession number: [MTBLS2892](www.ebi.ac.uk/metabolights/MTBLS2892) 

MAW is divided into MAW-R and MAW-Py. MAW-R was performed with two scripts:
1. Ponder-Smarinoi-R-Spectra.ipynb performed spectral database dereplication using GNPS, HMDB and MassBank [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6528931.svg)](https://doi.org/10.5281/zenodo.6528931). 
2. Ponder-Smarinoi-R-SIRIUS.ipynb performed dereplication using the structures present in the "ALL" database from SIRIUS4.

MAW-Py used Workflow_Python_Script_all.py to perform the candidate selection from spectral databases and SIRIUS, while the SMILES_postprocess.py was used to finalize and curate the list of annotated metabolites.

# Suspect List
The suspect list curation was performed using the suspectlits_curation.py. The comments in the script provide detail on all steps involved in the curation process and can be applicable to other suspect lists or any list with SMILES.

## Citation
[1] Zulfiqar, M., Gadelha, L., Steinbeck, C. et al. MAW: the reproducible Metabolome Annotation Workflow for untargeted tandem mass spectrometry. J Cheminform 15, 32 (2023). https://doi.org/10.1186/s13321-023-00695-y
