# Scripts for the paper *Orientation dependent pinning of (sub)grains by dispersoids during recovery and recrystallization in an Al-Mn alloy*

A *Zenodo* DOI is added here after the first 'release' of this repository.

This repository contains *Jupyter* notebooks, *MTEX* (*MATLAB*) and *ImageJ*scripts, and other files, apart from the raw BSE images and EBSD datasets, which are necessary to reproduce the results and figures in the paper *Orientation dependent pinning of (sub)grains by dispersoids during recovery and recrystallization in an Al-Mn alloy* (2022), submitted to X:

```bibtex
@article{entry,
  author   = {Håkon W. Ånes and Antonius T. J. {van Helvoort} and Knut Marthinsen},
  year     = {2022},
}

```

The preprint is available on *arXiv* (doi link):

```bibtex
@article{entry,
  author   = {Håkon W. Ånes and Antonius T. J. {van Helvoort} and Knut Marthinsen},
  year    = {2022},
}
```

The raw EBSD and BSE data is available on *Zenodo* ([doi](https://doi.org/10.5281/zenodo.7383087)):

```bibtex
@dataset{entry,
  year      = {2022},
}
```

The contents in this repository is licensed under the GPLv3+, since many of the softwares used have the same license.

## Contents

### notebooks

Static views of the notebooks are available via [nbviewer](https://nbviewer.org/github/hakonanes/p-texture-al-mn-alloys/tree/main/notebooks/).

*Python* packages used in the notebooks are listed in `requirements.txt` and can be installed into a virtual or *conda* environment:

```bash
pip install -r requirements.txt
```

Files:

* `ebsd1_dewrap.ipynb`: EBSD datasets acquired with a *NORDIF* detector are sometimes written with the last column of patterns as the first column. This notebook checks if this is the case, places the patterns correctly in the map, and writes them to an HDF5 file in the *kikuchipy* h5ebsd format.
* `ebsd2_preprocess.ipynb`: Generate indexing-independent views of an EBSD dataset (mean intensity map, virtual backscatter electron images, image quality map, and average neighbour dot product map) and calibrate the detector-sample geometry via projection center (PC) optimization with the [*PyEBSDIndex* *Python* package](https://github.com/USNavalResearchLaboratory/PyEBSDIndex). An average PC is used in dictionary indexing.
* `ebsd3_dictionary_indexing.ipynb`: Obtain crystal orientations from the EBSD patterns via dictionary indexing (DI) as implemented in kikuchipy. Requires an Al master pattern generated with *EMsoft*.
* `ebsd4_refinement.ipynb`: Refine crystal orientations obtained from DI.
* `bse1_crop_bse_images.ipynb`: Crop BSE images prior to and after stitching with *BigStitcher* in *ImageJ*.
* `bse2_particle_detection.ipynb`: Detect particles in BSE images and EBSD intensity maps.
* `bse3_image_registration.ipynb`: Manually identify control points in BSE images and EBSD intensity maps, perform image registration of BSE and EBSD datasets, and insert particles detected in BSE images and their sizes into EBSD datasets, making up the multimodal datasets.
* `bse4_density_disp_by_const_part.ipynb`: Calculate the density and average size of dispersoids within deformation zones of large constituent particles relative to elsewhere in the microstructure.
* `bse4_density_disp_by_const_part.py`: The above notebook as a *Python* script accepting input parameters.
* `run_bse4_density_disp_by_const_part.sh`: Shell script to run the above *Python* script per dataset with varying parameter values. Useful when doing a parameter search or just automating the calculations.
* `figures_in_paper.ipynb`: Generate most of the figures used in the paper (and some for the supplementary).
* `figures_in_supplementary.ipynb`: Generate most of the figures in the supplementary.

Installation of the packages has been tested to work on *Linux* (*Ubuntu* 22.04) and *Windows* 10. All notebooks have been tested to work on *Linux*, while the two indexing and refinement notebooks have also been tested to work on *Windows* 10.

### matlab_scripts

*MATLAB* packages used in the scripts are [*MTEX*](https://mtex-toolbox.github.io/) and [*export_fig*](https://mathworks.com/matlabcentral/fileexchange/23629-export_fig).
Note that *MATLAB* requires a license.

Files:

* `macrotexture.m`: Estimation of the orientation distribution function from four pole figures collected with an XRD diffractometer.
* `orientation_analysis.m`: Correlated analysis of subgrains and particles based on the multimodal dataset, specifically dispersoids at subgrain boundaries and subgrains at constituent particles.

### imagej_scripts

Files:

* `hist_match_subtract_bg.bsh`: *ImageJ* *BeanShell* script automating histogram matching and creation of images to detect large and small particles in a BSE image.
