# PSI
Python implementation of Phase Sorting Interferometry (PSI) for estimating non-common path errors in High Contrast Imaging (HCI) instruments.

This PSI package is developed primarily for ELT/METIS, but also for VLT/ERIS.


## Legacy contribution
This work is based on the initial work of Emiel Por, followed-up by the work of Matthew Willson†:
- [fepsi](https://github.com/mkenworthy/fepsi) : written by Emiel Por and modified by Matthew Kenworthy (Leiden University)
- [psi](https://github.com/mwillson-astro/PSI/tree/master) : preliminary package developed by Matthew Willson† for METIS (Liège University)

The initial commit of this repository is a copy of [psi](https://github.com/mwillson-astro/PSI/tree/master).

## References
- ["Focal Plane Wavefront Sensing Using Residual Adaptive Optics Speckles" by Codona and Kenworthy (2013)](https://iopscience.iop.org/article/10.1088/0004-637X/767/2/100),  ApJ, 767, 100.

## Dependencies

Code requires the following Python package : `hcipy` and all its dependencies.

The ffmpeg tools are required for generating movies.

[HCIPy](https://github.com/ehpor/hcipy) is a Python software package written and developed by Emiel Por for performing end-to-end simulations of high contrast imaging instruments for astronomy.

**HCIPy version currently required: '0.4.1.dev137+gc1ea9e5' (master version of 5/12/2021).**

**Do not pip install but do a git clone instead**
