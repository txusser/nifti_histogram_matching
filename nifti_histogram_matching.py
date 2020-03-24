import nibabel as nib
import numpy as np
from exact_hm.histogram_matching import ExactHistogramMatcher

def histogram_matching(reference_nii, input_nii, output_nii):

    # Load the template image
    template = nib.load(reference_nii)
    nt_data = template.get_data()[:,:,:]
    
    # Load the patient image
    patient = nib.load(input_nii)
    pt_data = patient.get_data()[:,:,:]

    # Stores the image data shape that will be used later
    oldshape = pt_data.shape

    # Converts the data arrays to single dimension and normalizes by the maximum
    nt_data_array = nt_data.ravel()
    pt_data_array = pt_data.ravel()

    # get the set of unique pixel values and their corresponding indices and counts
    s_values, bin_idx, s_counts = np.unique(pt_data_array, return_inverse=True, return_counts=True)
    t_values, t_counts = np.unique(nt_data_array, return_counts=True)

    # take the cumsum of the counts and normalize by the number of pixels to
    # get the empirical cumulative distribution functions for the source and
    # template images (maps pixel value --> quantile)
    s_quantiles = np.cumsum(s_counts).astype(np.float64)
    s_quantiles /= s_quantiles[-1]
    t_quantiles = np.cumsum(t_counts).astype(np.float64)
    t_quantiles /= t_quantiles[-1]

    # interpolate linearly to find the pixel values in the template image
    # that correspond most closely to the quantiles in the source image
    interp_t_values = np.interp(s_quantiles, t_quantiles, t_values)

    #Reshapes the corresponding values to the indexes and reshapes the array to input
    final_image_data = interp_t_values[bin_idx].reshape(oldshape)
    #final_image_data[indx] = 0

    #Saves the output data
    img = nib.Nifti1Image(final_image_data, patient.affine, patient.header)
    nib.save(img, output_nii)

    return output_nii


def logpow_histogram_matching(reference_nii, input_nii, output_nii, alpha = 1, beta = 3):

    # Load the template image
    template = nib.load(reference_nii)
    nt_data = template.get_data()[:,:,:]
    
    # Load the patient image
    patient = nib.load(input_nii)
    pt_data = patient.get_data()[:,:,:]

    # Stores the image data shape that will be used later
    oldshape = pt_data.shape

    # Converts the data arrays to single dimension and normalizes by the maximum
    nt_data_array = nt_data.ravel()
    pt_data_array = pt_data.ravel()

    # get the set of unique pixel values and their corresponding indices and counts
    s_values, bin_idx, s_counts = np.unique(pt_data_array, return_inverse=True, return_counts=True)
    t_values, t_counts = np.unique(nt_data_array, return_counts=True)

    s_counts = np.power(np.log10(s_counts + alpha),beta)
    t_counts = np.power(np.log10(t_counts + alpha),beta)

    # take the cumsum of the counts and normalize by the number of pixels to
    # get the empirical cumulative distribution functions for the source and
    # template images (maps pixel value --> quantile)
    s_quantiles = np.cumsum(s_counts).astype(np.float64)
    s_quantiles /= s_quantiles[-1]
    t_quantiles = np.cumsum(t_counts).astype(np.float64)
    t_quantiles /= t_quantiles[-1]

    # interpolate linearly to find the pixel values in the template image
    # that correspond most closely to the quantiles in the source image
    interp_t_values = np.interp(s_quantiles, t_quantiles, t_values)

    #Reshapes the corresponding values to the indexes and reshapes the array to input
    final_image_data = interp_t_values[bin_idx].reshape(oldshape)
    #final_image_data[indx] = 0

    #Saves the output data
    img = nib.Nifti1Image(final_image_data, patient.affine, patient.header)
    nib.save(img, output_nii)

    return output_nii


def bi_histogram_matching(reference_nii, input_nii, output_nii, nbins = 256):

    # Load and prepare the template image
    template = nib.load(reference_nii)
    nt_data = template.get_data()[:,:,:]
    indx = np.where(nt_data < 0)
    # Set array values in the interval to val
    nt_data[indx] = 0
    nt_max = np.amax(nt_data)
    nt_data = ((nbins-1)/nt_max)* nt_data

    # Load and prepare the patient image
    patient = nib.load(input_nii)
    pt_data = patient.get_data()[:,:,:]
    indx = np.where(pt_data < 0)
    # Set array values in the interval to val
    pt_data[indx] = 0
    pt_max = np.amax(pt_data)
    pt_data = (255/pt_max)* pt_data
    oldshape = pt_data.shape
    
    # We adjust first the low part of the histogram
    nt_data_low = nt_data
    indx = np.where(nt_data > (nbins/2))
    nt_data_low[indx] = 0
    pt_data_low = pt_data
    indx = np.where(pt_data > (nbins/2))
    pt_data_low[indx] = 0
    # We save data_low just for checking for the moment
    img = nib.Nifti1Image(pt_data_low, patient.affine, patient.header)
    nib.save(img, input_nii[0:-4] + '_low.nii')

    # Converts the data arrays to single dimension and histograms the first part of the data
    nt_data_array = nt_data_low.ravel()
    pt_data_array = pt_data_low.ravel()
    s_values, bin_idx, s_counts = np.unique(pt_data_array, return_inverse=True, return_counts=True)
    t_values, t_counts = np.unique(nt_data_array, return_counts=True)
    s_quantiles = np.cumsum(s_counts).astype(np.float64)
    s_quantiles /= s_quantiles[-1]
    t_quantiles = np.cumsum(t_counts).astype(np.float64)
    t_quantiles /= t_quantiles[-1]
    interp_t_values = np.interp(s_quantiles, t_quantiles, t_values)
    final_image_data_low = interp_t_values[bin_idx].reshape(oldshape)
    
    # We prepare now to adjust the second part of the histogram
    nt_data = template.get_data()[:,:,:]
    indx = np.where(nt_data < 0)
    nt_data[indx] = 0
    nt_max = np.amax(nt_data)
    nt_data = ((nbins-1)/nt_max)* nt_data

    pt_data = patient.get_data()[:,:,:]
    indx = np.where(pt_data < 0)
    pt_data[indx] = 0
    pt_max = np.amax(pt_data)
    pt_data = ((nbins-1)/pt_max)* pt_data
    oldshape = pt_data.shape

    nt_data_high = nt_data
    indx = np.where(nt_data <= (nbins/2))
    nt_data_high[indx] = 0
    indx = np.where(nt_data > (nbins/2))
    nt_data_high[indx] = nt_data_high[indx]-(nbins/2)
    pt_data_high = pt_data
    indx = np.where(pt_data <= (nbins/2))
    pt_data_high[indx] = 0
    indx_pth = np.where(pt_data_high > (nbins/2))
    pt_data_high[indx_pth] = pt_data_high[indx_pth]-(nbins/2)

    # Converts the data arrays to single dimension and histograms the second part of the data
    nt_data_array = nt_data_high.ravel()
    pt_data_array = pt_data_high.ravel()
    s_values, bin_idx, s_counts = np.unique(pt_data_array, return_inverse=True, return_counts=True)
    t_values, t_counts = np.unique(nt_data_array, return_counts=True)
    s_quantiles = np.cumsum(s_counts).astype(np.float64)
    s_quantiles /= s_quantiles[-1]
    t_quantiles = np.cumsum(t_counts).astype(np.float64)
    t_quantiles /= t_quantiles[-1]
    interp_t_values = np.interp(s_quantiles, t_quantiles, t_values)
    final_image_data_high = interp_t_values[bin_idx].reshape(oldshape)
    
    final_image_data_high[indx_pth] = final_image_data_high[indx_pth]+128
    final_image_data = final_image_data_high + final_image_data_low

    #Saves the output data
    img = nib.Nifti1Image(final_image_data, patient.affine, patient.header)
    nib.save(img, output_nii)


    return output_nii


def exact_histogram_matching(reference_nii, input_nii, output_nii, number_kernels=3, nbins=1024):

      
    template = nib.load(reference_nii)
    nt_data = template.get_data()[:,:,:]
    scaled_nt_data = np.round((nbins-1) * (nt_data / np.max(nt_data)))
    scaled_nt_data = np.asarray(scaled_nt_data, np.uint16)

    patient = nib.load(input_nii)
    pt_data = patient.get_data()[:,:,:]
    scaled_pt_data = np.round((nbins-1) * (pt_data / np.max(pt_data)))
    scaled_pt_data = np.asarray(scaled_pt_data, np.uint16)

    reference_histogram = ExactHistogramMatcher.get_histogram(scaled_nt_data, 10)
    matched_img = ExactHistogramMatcher.match_image_to_histogram(pt_data, reference_histogram,number_kernels)

    img = nib.Nifti1Image(matched_img, patient.affine, patient.header)
    nib.save(img, output_nii)

    return output_nii

