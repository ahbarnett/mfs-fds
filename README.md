# mfs-fds

MFS via FDS.

authors: Alex Barnett, Ken Ho, and Yuxiang (Larry) Liu.

## Usage

### Dependencies:

* MATLAB, any recent version.
* [FLAM](https://github.com/klho/FLAM)
* [FMMLIB2D](https://github.com/zgimbutas/fmmlib2d)
* optionally: [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html)

If you place (eg soft-links to)
each of the directories for the above packages alongside this `mfs-fds`
directory,
there is no adjustment to do.  Otherwise adjust the location of your
distributions in the file `startup.m`.

### Example:

Open MATLAB from the `mfs-fds` directory

type `startup` (from this directory)

type `cd 2D`

type `scat_driver`, which will run a small example for a couple of
seconds, print some (hopefully small) solution error value, and
produce a 2D wave scattering plot.

## Contents

* `2D/scat_driver` : small smooth Dirichet scattering example  
* `2D/scat_amoeba_conv` : medium smooth Dirichlet scattering convergence tests  
* `2D/scat_amoeba_conv_nF100_k100' : medium-large Dirichlet smooth scattering + movie  
* `2D/scat_panel_onecorner_conv_demos` : small high-accuracy corner Dirichlet scattering examples
* `2D/run_all_onecorner_conv_demos` : sweeps through 5 {0,1}-corner small shapes as above
* `2D/scat_panel_polygon` : choice of two large-scale polygon examples  


