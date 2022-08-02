# Effect of concurrent precipitation on P texture in Al-Mn alloys

Zenodo DOI coming after first `release' of repository.

This repository contains Jupyter notebook files, MATLAB scripts, and other files, apart from the raw BSE images and EBSD datasets, which are necessary to reproduce the results and figures in the paper X, which was recently submitted to X. The paper preprint is available on [arXiv](). The raw EBSD and BSE data used in the paper is available on [Zenodo]().

The contents in this repository is licensed under the GPLv3+, since many of the softwares used have the same license.

## Contents

### notebooks

Static views of the notebooks are available via [nbviewer](https://nbviewer.org/github/hakonanes/p-texture-al-mn-alloys/tree/main/notebooks/).

Python packages used in the notebooks are listed in `requirements.txt` and can be installed into a virtual or conda environment:

```bash
pip install -r requirements.txt
```

* `ebsd1_dewrap.ipynb`: EBSD datasets acquired with a NORDIF detector are sometimes written with the last column of patterns as the first column. This notebook checks if this is the case, places the patterns correctly in the map, and writes them to an HDF5 file in the kikuchipy h5ebsd format.
* `ebsd2_preprocess.ipynb`: Generate indexing-independent views of an EBSD dataset (mean intensity map, virtual backscatter electron images, image quality map, and average neighbour dot product map) and calibrate the detector-sample geometry via projection center (PC) optimization with the [PyEBSDIndex Python package](https://github.com/USNavalResearchLaboratory/PyEBSDIndex) (cubic materials only!). An average PC is used in dictionary indexing.
* `ebsd3_dictionary_indexing.ipynb`: Obtain crystal orientations from the EBSD patterns via dictionary indexing (DI) as implemented in kikuchipy. Requires an Al master pattern generated with EMsoft.
* `ebsd4_refinement.ipynb`: Refine crystal orientations obtained from DI.
* `bse1_crop_bse_images.ipynb`: Crop BSE images prior to and after stitching with BigStitcher in ImageJ.
* `generate_figures_in_paper.ipynb`: Generate final figures used in the paper.
* `image_registration.ipynb`: Manually identify control points in BSE images and EBSD intensity maps, perform image registration of BSE and EBSD datasets, and insert particles detected in BSE images and their sizes into EBSD datasets, making up the multimodal datasets.
* `particle_detection.ipynb`: Detect particles in BSE images and EBSD intensity maps.

Installation of packages has been tested to work on Linux (Ubuntu 22.04) and Windows 10. All notebooks have been tested to work on Linux, while the two indexing and refinement notebooks have also been tested to work on Windows 10.
