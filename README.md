# Implementation of synthetic-eddy-method in Nek5000

This contains a NEK case for turbulent straight pipe flow. The turbulence is generated with the [DFSEM by Poletto et al.](https://www.escholar.manchester.ac.uk/uk-ac-man-scw:262700) and verified with the [DNS results of El Khoury et al.](http://link.springer.com/article/10.1007%2Fs10494-013-9482-8).

The DFSEM is implemented to generate divergence-free, nonuniform, isotropic turbulence like the DFSEMiso in Poletto's thesis. 

#### Parameters/Limitations
* An input file `sem_input.txt`is required, that gives the mean-velocity, turbulent kinetic energy and dissipation rate and the wall distance of these scalars. This can be generated with `statistics/matlab_scripts/generate_k_eps.m`
* The following parameters are taken from `pipe.upar` 
    * `neddy` the number of eddies. Should *(probably??)* suffice A_inlet/(πσ_min²). However, it was observed that too many eddies (at too small sigma_min) seem not to contribute to the generated signal properly -> underestimation of turbulent kinetic energy, because of factor sqrt(Vb/neddy) 
    * `nElInlet` the number of elements in the inlet plane
    * `sigma_max` The maximal eddy size, i.e. max(k^1.5/eps). Can be taken from above Matlab script
    * `u0` The bulk mean velocity
    * `Vb` and `z_inlet` are being calculated in `usrdat()`. They are the volume of the convective box/cylinder and the `z` coordinate of the inlet plane respectively. Futher, the wall distance (or here pipe radius `ybmax`) is required to bound eddy sizes.
 
As of now, the DFSEM is only implemented to generate a turbulent signal in z-positive direction. It should however be straightforward to change the relevant code (x-positive SEM was before [commit c820c3](../../commit/c820c3d9f9ae82491efa70bcbb80dae23970e9b5) ). The first elements in the mesh need to be the face, at which the DFSEM inflow is placed. This enables us to precompute eddy-size once for the whole simulation (using linear interpolation from `sem_input.txt`) instead of recomputing the eddy-sizes at each time step. An unoptmized version without this limitation is found in [commit aaa095](../../commit/aaa095)  

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
