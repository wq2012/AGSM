# AGSM: Active Geometric Shape Models

## Copyright

Copyright (C) 2012 Quan Wang <wangq10@rpi.edu>,
Signal Analysis and Machine Perception Laboratory,
Department of Electrical, Computer, and Systems Engineering,
Rensselaer Polytechnic Institute, Troy, NY 12180, USA

You are free to use this software for academic purposes if you cite our paper:
```
Quan Wang, Kim L. Boyer,
The active geometric shape model: A new robust deformable shape model and its applications,
Computer Vision and Image Understanding, Volume 116, Issue 12, December 2012, Pages 1178-1194,
ISSN 1077-3142, 10.1016/j.cviu.2012.08.004.
```

For commercial use, please contact the authors.

## Overview

This software implements active geometric shape models for line fitting,
circle fitting, ellipse fitting, and cubic spline fitting.

This software requires Matlab and Matlab compiler.

Before using, you must compile. Just run the compile.m file under the code
directory.

Four examples of using this software:
```
code/line_fitting/demo_line_fitting.m
code/circle_fitting/demo_circle_fitting.m
code/arbitrary_ellipse_fitting/demo_ellipse_fitting.m
code/spline_contour_fitting/demo_spline_contour_fitting.m
```

You can also play with the AGSM Canvas app. It runs on 64-bit Windows or 64-bit Mac.

Current version: 2.3

## Useful Links

For patient privacy, the CSF PC-MR images used in our paper cannot be provided in this package.


The Matlab version of the GVF algorithm can be downloaded at:
* http://www.iacl.ece.jhu.edu/static/gvf/


The s-t graph cuts algorithm can be downloaded at:
* http://vision.csd.uwo.ca/code/


The active contours algorithm that we used for comparison can be downloaded at:
* http://www.engr.uconn.edu/~cmli/DRLSE/

Please also visit our project website for updates:

* https://sites.google.com/site/agsmwiki/

## History

### New features in AGSM version 2.3

1. Robustness of closed contour fitting is dramatically improved by smart initialization (see InitialCircle.cpp).

2. Line fitting is compared to Hough transform and RANSAC least squares.

3. Circle fitting is compared to circle Hough transform.

### New features in AGSM version 2.2

1. We made the AGSM Canvas app more robust by smart initialization and image padding.

2. Only one app for both 64-bit Windows and 64-bit Mac.

### New features in AGSM version 2.1

1. We added the AGSM Canvas app, which is an interactive software.

### New features in AGSM version 2.0

1. We added the line fitting and circle fitting packages.

2. We re-implemented GVF in C++/MEX to make it much faster.

3. We implemented Bresenham algorithm in C++/MEX to generate points on a line.

4. All force/torque computations are re-implemented in C++/MEX, which significantly improves efficiency.

