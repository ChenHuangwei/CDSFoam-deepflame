# CDSFoam: An detonation solver based on deepflame
1. The fluxschemes is form [blastFoam](https://github.com/synthetik-technologies/blastfoam) and [detonationFoam](https://github.com/JieSun-pku/detonationFoam)
2. The HLLC-LM scheme is introduced
3. The explicit third-order SSP Runge-Kutta method for time integration
4. The CVODE solver of Cantera for chemical reaction rate evaluation
5. The dynamic load balancing (DLB) is used for chemical solution

## How to install
1. Install [OpenFOAM-7](https://openfoam.org/version/7/)
2. Install [deepflame](https://github.com/deepflameCFD/deepflame-dev)
3. Compile dfCDSFoam:```./Allwmake```

## Derived cases
1. ***Riemann problem***

   ![Riemann problem](https://github.com/ChenHuangwei/CDSFoam-deepflame/blob/master/Figs/Riemann.png)
2. ***Deflagration to detonation transition***

   ![DDT](https://github.com/ChenHuangwei/CDSFoam-deepflame/blob/master/Figs/DDT.png)
4. ***Rotating detonation engine***

   ![RDE](https://github.com/ChenHuangwei/CDSFoam-deepflame/blob/master/Figs/RDE.png)
6. ***Oblique detonation***

   ![ODW](https://github.com/ChenHuangwei/CDSFoam-deepflame/blob/master/Figs/ODW.png)

## Citation
>**Implementation and verification of an OpenFOAM solver for gas-droplet two-phase detonation combustion(Revised)**
