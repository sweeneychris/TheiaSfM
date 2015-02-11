=============================================

            data_loader directory

=============================================

This directly contains utility files for dealing with the Big SfM localization
datasets Dubrovnik6k and Rome16k:
http://www.cs.cornell.edu/projects/bigsfm/#data

The datasets provide visual information after an SfM reconstruction (with bundle
adjustment). This includes camera pose, camera intrinsic parameters, 3D points,
2D-3D matches, and 2D descriptors. The datasets are split into two parts: db and
query. The db contains a reconstruction of *only* the db images with no links to
the query images. The query images provide 2D descriptors along with ground truth
pose provided for comparison.

NOTE: The "orig" component of the datasets ar not consistent with the query
descriptors. The orig component uses query images with higher resolution than
the query images provided with the dataset. As such, more descriptors per query
image were used in the reconstruction and the true 2D-3D correspondences are
different than what is provided (SIFT keys for the smaller query images are
given, which do not align with the SIFT keys from the larger resolution query
images used for the orig reconstruction).


=============================================
	      What is provided
=============================================

This directory contains programs to read the "db" component of the Dubrovnik6k
and Rome16k datasets. The "orig" datasets are not included at this time for the
reasons discussed above. Additionally, a script is provided to convert the
(massive) text files provided by the datasets to binary files. In my experience,
this reduces the size of the dataset by roughly 4x, and is much faster to load.

Additionally, we provide methods to load the binary versions of the files and a
very rudimentary application to view the reconstruction in OpenGL. The viewer is
provided mostly as a verification that the reconstructions were loaded properly.


=============================================
            Questions or Issues
=============================================

If you have any questions or issues contact the author of these files, Chris
Sweeney, at cmsweeney@cs.ucsb.edu
