# Implementation of synthetic-eddy-method in Nek5000

This contains a NEK case for turbulent straight pipe flow. The turbulence is generated with the [SEM by Jarrin et al.](http://cfd.mace.manchester.ac.uk/desider/symposium/symp05/Session_2/S2-B.pdf) and verified with the [DNS results of El Khoury et al.](http://link.springer.com/article/10.1007%2Fs10494-013-9482-8).

The SEM is implemented to generate nonuniform, isotropic turbulence, like in case `SEM3` in Jarrin's paper. 
#### Parameters/Limitations
* An input file `sem_input.txt`is required, that gives the mean-velocity, turbulent kinetic energy and dissipation rate and the wall distance of these scalars. This can be generated with `statistics/matlab_scripts/generate_k_eps.m`
* The following parameters are taken from `pipe.upar` 
    * `neddy` the number of eddies. Should *(probably??)* suffice A_inlet/(πσ_min²). However, it was observed that too many eddies (at too small sigma_min) seem not to contribute to the generated signal properly -> underestimation of turbulent kinetic energy, because of factor sqrt(Vb/neddy) 
    * `nElInlet` the number of elements in the inlet plane
    * `sigma_max` The maximal eddy size, i.e. max(k^1.5/eps). Can be taken from above Matlab script
    * `yplus_cutoff` The wall distance below which no synthetic eddies will be created (to satisfy divergence-free condition). Should be == 10 in plus units
    * `bbox_max` The additional size (wall-normal) of the box/cylinder in which synthetic eddies are convected, i.e. such that compact support of all eddies is included. Can be taken from above Matlab script
    * `u0` The bulk mean velocity
    * `Vb`, `z_inlet`, and `ybmax` are be set in `usrdat()`. They are the volume of the convective box/cylinder, the `z` coordinate of the inlet and the _radius_
 
As of now, the SEM is only implemented to generate a turbulent signal in z-positive direction. It should however be straightforward to change the relevant code (x-positive was before [commit c820c3](../../commit/c820c3d9f9ae82491efa70bcbb80dae23970e9b5) ). The first elements in the `pipe.rea` should be the face, at which the SEM inflow is placed. This enables us to precompute eddy-size once for the whole simulation (using linear interpolation from `sem_input.txt`) instead of calculate the eddy-size at each time step. An unoptmized version without this limitation shoud be in [commit aaa095](../../commit/aaa095)  

**Future outlook**
In the [SEM as implemented in Code Saturne by Jarrin](http://code-saturne.org/viewvc/saturne/trunk/src/turb/cs_les_inflow.c?view=markup), there is a `rescale_flowrate`, which normalizes the mass-flux of the fluctuations s.t. they add zero bulk-net to the mean flow. This is obviously especially relevant for a small number of eddies.

## Workflow
* Create a cartesian pipe-mesh with positive z as downstream-direction. Typically, the mesh generator [pipeMeshNek by Jacopo](https://bitbucket.org/jacopo-canton/pipemeshnek) was used  
This will yield files `base.rea, base2d.rea`. The first one contains the actual mesh, the latter one a single pipe crosssectionial face, i.e. it contains `nElInlet` 2d elements. It is required for the averaging/statistics code.

* run `reatore2, genmap, makenek`
* run `nek5000` 

###Statistics 
* generate a (2D) crossectional pipe-face in polar coordinates using Matlab script `statistics/matlab_scripts/generate_polar_mesh.m`. Nek in postproc-mode will interpolate the time-averaged statistics onto this  
* compile and run `stats.usr`-nek5000 in subfolder `statistics`, using `base2d.rea` as `stats.rea` 
* Run Matlab script `statistics/matlab_scripts/plot_stats.m` to evaluate stastics

Lorenz Hufnagel, hufnagel@kth.se
