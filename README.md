# openip - The open source image processing library

## Introduction

This is a lightweight image processing library developed to support academic, research and industrial applications.

## How to build

Building the library requires:
* <code>cmake</code> - to manage the build process,
* <code>libgsl</code> - for numerical methods,
* <code>libpng</code> - for PNG IO,
* <code>libjpeg</code> - for JPG IO,
* <code>libtiff</code> - for TIFF IO.
    
## Components

The library has 5 main components:
* <code>openipDS</code> - the core data structures like Vector, Image, Volume, Filter, etc.;
* <code>openipIO</code> - the IO interfaces;
* <code>openipSC</code> - interfaces to GSL;
* <code>openipLL</code> - high level image processing algorithms, like skeletonization, segmentation, etc.;
* <code>openipML</code> - machine learning techniques, like decision trees, Support Vector Machines, etc.
    
## Applications

There is a sample application to demonstrate the use of the main components for template matching. Once the application is built, one can apply template matching using various (dis)similarity measures by calling the application as <code>templateMatching --pearsoncc input.png template.png output.png</code>

