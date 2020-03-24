# nifti_histogram_matching
This repository includes different methods for performing histogram matching between 2 (reference and source) medical images. 

# Setup
1. Once you downloaded the repository, you can run setup.py (python setup.py) to install dependencies and unpack resources.
2. Once installed, you can import nifti_histogram_matching as a Python library.

# Example

import nifti_histogram_matching as nhm
nhm.exact_histogram_matching('reference_image.nii','source_image.nii','test_ehm.nii')

Example data can be found in resources/examples

# Methods included:

Conventional HM (HM): Source and reference images (IS and IR) are histogrammed mapping the voxel intensity values to discrete bins. The number of bins was 230, calculated by the Rice Rule. Following, the discrete cumulative distribution functions (CDFs) for each of the histograms are calculated. Let us denote the CDF of the template image as CR, while CS is the CDF of the source image. For each value in the source histogram Vs, a value VR in the template histogram must be found, by simple interpolation, such that CS (VS) = CR (VR). This will produce a mapping of VR and VS values, which is then applied to modify the source image to the reference histogram space, providing the intensity-harmonized image.

Exact HM (EHM): Because the calculated CDFs are based on discrete bins, the CDFs used in HM are not exactly invertible, making histogram matching an ill-posed problem.  Different solutions have been proposed for this. Here, we have modified a previous Python implementation of the exact HM method proposed by Coltuc et al. (Coltuc et al., 2006) (https://github.com/StefanoD/ExactHistogramSpecification) in order to accept 3D medical images as an input. The implemented EHM method tries to obtain invertible cumulative distribution functions by translating the problem into a k-dimensional space and further inducing a strict ordering among image pixels values. Thus, Let IS (IR) be the source (reference) image, with (M x N x L) pixels, and let HS (HR) = {(V0, h0), (V1, h1), (V2, h2)…, (V229, h229)} the image histogram, where h0 is the number of voxels in the first histogram bin and V0 is the histogram center value. Let us denote a strict ordering relation of order k among the voxels on each image: I(x1,y1,z1) → I(x2,y2,z2) → … → I(xM,xN,xL). Such an ordering can be obtained in different ways. Then, each of the ordered strings (source and reference) is splitted from left to right in groups, so that group j has hj voxels. At last, for each group, the mapping is performed as VS (IS (j)) = VR (IR (j)).

Bi-Histogram HM (BiHM): One of the potential problems of HM and EHM is that they can produce artificially enhanced contrast when applied to images with reduced value ranges, for example, for pathological cases where both striatums have reduced bilaterally the binding. This is a problem inherent to HM, and several authors have tried to develop brightness-preserving versions of HM (Kim, 1997; Ooi & Mat Isa, 2010; Tang & Mat Isa, 2017). Here we have developed a custom implementation of the bi-histogram method proposed by Kim (Kim, 1997). The idea behind the algorithm is to utilize independent histogram matching separately over two or more sub-images.  In our implementation, the input image values are splitted in half based on the distribution mean. HM is performed over the resulting image segments, with a constraint that the resulting matched sub-images are bounded to each other around the input mean. Finally, the resulting HM matched images are added up to generate the final harmonized image.

Logpow HM (LpHM): Another implementation of HM that provides brightness preservation features is the log-pow HM implementation proposed by Toet and Wu (Toet & Wu, 2014). The method consists on applying a logarithm operator that reduces the effect of spikes and noise, smoothing the image histogram while retaining the relative size ordering of the original bins. Then, a power operator is applied to obtain a smoothed version of the original histogram. After the log-pow process, conventional HM is applied. The method proposes the following modification:

H^' [i]=〖(log⁡(H[i]+ α)〗^β

where h and h’ represent the original and modified histograms, α is a constant larger than 1, introduced to avoid taking the algorithm of zero on empty bins, and β is the exponent of the power function. For our application, we have fixed α=1 and β=2.5, values found as conservative by the original authors (Toet & Wu, 2014).
