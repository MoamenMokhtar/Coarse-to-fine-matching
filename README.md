# DESCRIPTION
-----------------------

We present a new approach to wide baseline matching. We propose to use a hierarchical decomposition of the image domain and coarse-to-fine selection of regions to match. In contrast to interest point matching methods, which sample salient regions to reduce the cost of comparing all regions in two images, our method eliminates regions sys​systematically to achieve efficiency. One advantage of our approach is that it is not restricted to covariant salient regions, which is too restrictive under large viewpoint and leads to few corresponding regions. Affine invariant matching of regions in the hierarchy is achieved efficiently by a coarse-to-fine search of the affine space. experiments on two benchmark datasets shows that our method finds more correct correspondence of the image with fewer false alarms than other wide baseline methods on large viewpoint change.



This code is the implementation of the algorithm provided in the following paper:
 > Yang, Yanchao, Zhaojin Lu, and Ganesh Sundaramoorthi. "Coarse-to-Fine Region Selection and Matching." Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition. 2015.

--------------------
 USAGE
=================

***NOTE***: This code was tested only on Windows 64-bits on Visual Studio 2010.
Requirements: **Opencv**

In order to setup the code on Windows, please follow the following steps:

> 1. Make sure you have OpenCV downloaded on your machine.
> 2. Clone the project on your machine, and create an empty Win 32 console application on Visual Studio.
> 3. Add the cloned source code and headers to your project.
> 4. Make sure to change your project settings to use Opencv.
> 5. This code uses Open MP for running in parallel, so make sure to set the "**Open MP support**" flag to "**true**".

> In order to run the code, you have to create **Params.txt** and fill it as shown in the Params.template file.


