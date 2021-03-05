# mfs-fds: 3D case 

MFS via FDS.

authors: Alex Barnett, Ken Ho, and Yuxiang (Larry) Liu.

## Usage

### Dependencies:

* MATLAB, any recent version.
* [FLAM](https://github.com/klho/FLAM)
* [FMMLIB2D](https://github.com/zgimbutas/fmmlib2d) v 1.2
* [FMMLIB3D](https://github.com/zgimbutas/fmmlib3d) v 1.2
* optionally: [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html)


### Example:

* scat_sphere: scattering off sphere
* scat_smooth: arbitrary smooth shape, a BOR shape is chosen, but no FFT trick is used. 
* scat_ellipsoid: scattering off ellipsoid 
* Results are saved in correspoinding .txt files

## Params

MFS-
* M: total number of boundary points 
* N: total number of source points (M = 1.2 * N)
* P: number of disretization in \phi direction 
* a, b, w: params for the boundary curve
* d: controls the location of the source points 

FLAM-
* rank_or_tol: tolerance 
* rp: radius of the proxy sphere 
* p: p-by-p/2 points on the proxy sphere, used quasi-uniform scheme from 
     https://github.com/ahbarnett/BIE3D/blob/master/utils/setupspherequad.m
* occ: another FDS param, controls max pts per box for quadtree
