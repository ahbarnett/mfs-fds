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
produce a 2D scattering plot.

## Contents

* `2D/scat_driver` : small smooth scattering example in 2D
* `2D/scat_panel_corner_dr` : small corner scattering example in 2D


