## README - A-KAZE Features

This library is an implementation of the AKAZE feature detector and descriptor
algorithms with the OpenCV dependency removed. This allows for more flexible use
of the AKAZE features without requiring the large dependency of OpenCV. The only
dependency for this library is the Eigen matrix library. We made made an effort
to keep the code relatively in tact to the original implementation so that it
would be easier to follow and compare.

The major difference are:
- everything is in the namespace libAKAZE
- the keypoints and descriptors are a minimal, custom struct instead of OpenCV's structs
- images are held in a row-major Eigen matrix with floating point pixels
- the number of threads can be set at runtime

In my experiments, this implementation is less than 1.5x slower than the OpenCV
implementation. This is because OpenCV has many hand-written SIMD optimizations
for image processing. However, this implementation does give the same keypoints
and descriptors (within 0.1% error) as OpenCV.

**NOTE:** If you have suggestions for how to improve performance (namely the image
convolution and image half-sampling) please email me or make a pull request!

I would also like to thank Pablo Alcantarilla and Pierre Moulon for help and encouragement in creating this fork of akaze.

## Contact Info

If you have questions about this particular implementation, find a bug, or have an improvement suggested please email Chris Sweeney at cmsweeney@cs.ucsb.edu

If you have questions about AKAZE, Pablo Alcantarilla (the primary author of the original work) would like to know! If you work in a research institution, university, company or you are a freelance and you are using KAZE or A-KAZE in your work, please send him an email!! Here is his contact information:

Pablo F. Alcantarilla
email: pablofdezalc@gmail.com

## Citation

If you use this code as part of your work, please cite the following papers:

1. **Fast Explicit Diffusion for Accelerated Features in Nonlinear Scale Spaces**. Pablo F. Alcantarilla, J. Nuevo and Adrien Bartoli. _In British Machine Vision Conference (BMVC), Bristol, UK, September 2013_

2. **KAZE Features**. Pablo F. Alcantarilla, Adrien Bartoli and Andrew J. Davison. _In European Conference on Computer Vision (ECCV), Fiorenze, Italy, October 2012_

The original AKAZE implementation is housed at:
`https://github.com/pablofdezalc/akaze`

## Library Dependencies

The code only relies on the Eigen matrix library (version 3.2.0 or higher)

If you want to use OpenMP parallelization you will need to install OpenMP in your system
In Linux you can do this by installing the gomp library

Tested compilers:
- GCC 4.8
- Apple Clang 3.5

Tested systems:
- Ubuntu 14.04
- Mac OS X Yosemite

## Getting Started

Compiling:

1. `$ mkdir build`
2. `$ cd build>`
3. `$ cmake ..`
4. `$ make`

If the compilation is successful you should see two executables in the folder bin:
- `akaze_features`
- `akaze_match`

Additionally, the library `libAKAZE[.a, .lib]` will be created in the folder `lib`.

If there is any error in the compilation, perhaps some libraries are missing.
Please check the Library dependencies section.

Examples:
To see how the code works, examine the two examples provided.

## Computing A-KAZE Features

For running the program you need to type in the command line the following arguments:
`./akaze_features img.jpg [options]`

The `[options]` are not mandatory. In case you do not specify additional options, default arguments will be
used. Here is a description of the additional options:

- `--verbose`: if verbosity is required
- `--help`: for showing the command line options
- `--num_threads`: number of threads AKAZE will use (default is 1)
- `--soffset`: the base scale offset (sigma units)
- `--omax`: the coarsest nonlinear scale space level (sigma units)
- `--nsublevels`: number of sublevels per octave
- `--diffusivity`: diffusivity function `0` -> Perona-Malik 1, `1` -> Perona-Malik 2, `2` -> Weickert
- `--dthreshold`: Feature detector threshold response for accepting points
- `--descriptor`: Descriptor Type, 0-> SURF_UPRIGHT, 1->SURF
                                   2-> M-SURF_UPRIGHT, 3->M-SURF
                                   4-> M-LDB_UPRIGHT, 5->M-LDB
- `--descriptor_channels`: Descriptor Channels for M-LDB. Valid values: 1, 2 (intensity+gradient magnitude), 3(intensity + X and Y gradients)
- `--show_results`: `1` in case we want to show detection results. `0` otherwise

## Important Things:

* Check `AKAZEConfig.h` in case you would like to change the value of some default settings
* The **k** constrast factor is computed as the 70% percentile of the gradient histogram of a
smoothed version of the original image. Normally, this empirical value gives good results, but
depending on the input image the diffusion will not be good enough. Therefore I highly
recommend you to visualize the output images from save_scale_space and test with other k
factors if the results are not satisfactory

## Image Matching Example with A-KAZE Features

The code contains one program to perform image matching between two images. The
ground truth homography files are given in the datasets directory and must be
provided at runtime

For running the program you need to type in the command line the following arguments:
`./akaze_match img1.jpg img2.pgm homography.txt [options]`

The `[options]` are not mandatory. In case you do not specify additional
options, default arguments will be used.

The datasets folder contains the **Iguazu** dataset described in the paper and
additional datasets from Mikolajczyk et al. evaluation.  The **Iguazu** dataset
was generated by adding Gaussian noise of increasing standard deviation.

For example, with the default configuration parameters used in the current code
version you should get the following results:

```
./akaze_match ../../datasets/iguazu/img1.pgm
              ../../datasets/iguazu/img4.pgm
              ../../datasets/iguazu/H1to4p
              --descriptor 4
```

```
Number of Keypoints Image 1: 1823
Number of Keypoints Image 2: 2373
A-KAZE Features Extraction Time (ms): 304.796
Matching Descriptors Time (ms): 54.1619
Number of Matches: 1283
Number of Inliers: 1047
Number of Outliers: 236
Inliers Ratio: 81.6056
```
