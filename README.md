# data-fusion
This repository contains the MATLAB package implementing the data fusion algorithm on a computational example and experimental cross-sections of developing fly embryos. It is associated to the article "Synthesizing developmental trajectories by semi-supervised learning", Villoutreix P., Andén J., et al., 2017.

The code is written in MATLAB. The datasets can be found at <https://github.com/paulvill/data-fusion-images>.

The code uses the ScatNet toolbox implementing the scattering transform, all the necessary functions are stored in the folder scatnet-0.2. More details can be found at <http://www.di.ens.fr/data/software/scatnet/>.

Some of the code in this package has been published previously:
- Dsilva, Carmeline J., et al. "Temporal ordering and registration of images in studies of developmental dynamics." Development 142.9 (2015): 1717-1724.
    Scripts are stored in the folder DLL15.
- Lederman, Roy R., and Ronen Talmon. "Common manifold learning using alternating-diffusion."  Tech. Report YALEU/DCS/TR-1497 (2014).
    Scripts are stored in the folder LT14.

To test the functionality of this package, open MATLAB and type
```matlab
run_examples;
```
This will initialize the package and run a toy example reconstructing a one-dimensional trajectory embedded in three dimensions. To reproduce the results presented in the paper, type
```matlab
run_experiments;
```
This script downloads an an experimental dataset consisting of live imaging movies of nuclei morphology and fixed snapshots containing both morphology and the spatial distribution of five chemical species. It then applies the harmonic extension method to color the movie, overlaying the spatiotemporal dynamics of those chemical species on top of the morphological changes.
