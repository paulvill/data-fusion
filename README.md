# data-fusion
This repository contains the MATLAB library implementing the data fusion algorithm on a computational example and experimental cross-sections of developing embryos. It is associated to the article "Synthesizing developmental trajectories by data fusion", Villoutreix P., Andén J., et al., 2017.

The code is written in Matlab. The datasets can be found at this URL. 

The code uses the ScatNet library implementing the Scattering Transformation which can be obtained at this adress http://www.di.ens.fr/data/software/scatnet/ .

Some of the scripts used have been published previously:
- Lim, Bomyi, et al. "Dynamics of inductive ERK signaling in the Drosophila embryo." Current Biology 25.13 (2015): 1784-1790.
    Scripts are stored in the folder LDL15.
- Dsilva, Carmeline J., et al. "Temporal ordering and registration of images in studies of developmental dynamics." Development 142.9 (2015): 1717-1724.
    Scripts are stored in the folder DLL15.
- Lederman, Roy R., and Ronen Talmon. "Common manifold learning using alternating-diffusion."  Tech. Report YALEU/DCS/TR-1497 (2014).
    Scripts are stored in the folder LT14.
    
The main scripts are 
- spiral_example.m : implements the data fusion algorithm on a non-linear 1-dimensional trajectory in a 3-dimensional space, knicknamed "spiral".
- spiral_example_cross_validation.m : implements the K-fold cross-validation on the "spiral".
- experimental_dataset.m : implements the data fusion on the experimental datasets leading to multimodal movie contain the spatio-temporal dynamics of 5 chemical species
- experimental_dataset.m : implements the K-fold cross validation on the "experimental datasets".
