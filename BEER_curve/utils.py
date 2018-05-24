from scipy.signal import medfilt
import numpy as np
from statsmodels.robust.scale import mad

__all__ = ['bindata', 'median_boxcar_filter']

def median_boxcar_filter(data, window_length=None, endpoints='reflect'):
    """
    Creates median boxcar filter and deals with endpoints

    Parameters
    ----------
    data : numpy array 
       Data array
    window_length: int
        A scalar giving the size of the median filter window
    endpoints : str
        How to deal with endpoints. 
        Only option right now is 'reflect', which extends the data array
        on both ends by reflecting the data

    Returns
    -------
    filter : numpy array
        The filter array
    """

    filter_array = data
    # Create filter array
    if(endpoints == 'reflect'):
        last_index = len(data) - 1
        
        filter_array = np.concatenate((np.flip(data[0:window_length], 0), 
            data, 
            data[last_index - window_length:last_index]))

        # Make filter
        # Check that window_length is odd
        if(window_length % 2 == 0):
            window_length += 1
        filt = medfilt(filter_array, window_length)

        filt = filt[window_length:window_length + last_index + 1]

    return filt

def bindata(time, data, binsize, bin_calc='median', err_calc='mad'):
    """
    Bins data array

    Parameters
    ----------
    time : numpy array
        Array of times
    data : numpy array
        data array
    binsize : float
        Width of bins in same units at time
    bin_calc : str
        Method to use to calculate datum in each bin. 
        Can be either 'mean' or 'median'
    err_calc : str
        Method to use to calculate uncertainty in each bin. 
        Can be either 'std' or 'mad'
        Default ('mad') to using 1.4826 x median absolute deviation --
        https://en.wikipedia.org/wiki/Median_absolute_deviation

    Returns
    -------
    binned_time : numpy array
        Time binned
    binned_data : numpy array
        Data binned
    binned_err : numpy array
    """

    # 2018 May 23 - There are not always points in each time bin,
    #   so we will TRY to find points but will not always find them.
    times_to_try = np.arange(np.min(time) + 0.5*binsize, 
            np.max(time) - 0.5*binsize, binsize)

    binned_time = np.array([])
    binned_data = np.array([])
    binned_err = np.array([])

    if(bin_calc == 'median'):
        bin_calc_func = np.nanmedian
    elif(bin_calc == 'mean'):
        bin_calc_func = np.nanmean

    if(err_calc == 'mad'):
        err_calc_func = lambda x : mad(x)/np.sqrt(len(x))
    elif(err_calc == 'std'):
        err_calc_func = lambda x : np.nanstd(x)/np.sqrt(len(x))

    for i in range(len(times_to_try)):
        ind = np.argwhere(np.abs(times_to_try[i] - time) <= binsize)

        if(ind.size > 0): 
            # Remove nans, too
            cur_data = data[ind[~np.isnan(data[ind])]] 
            
            if(cur_data.size > 0):
                binned_time = np.append(binned_time, times_to_try[i])

                binned_data = np.append(binned_data, bin_calc_func(cur_data))

                # Check for bad error value
                try_error = err_calc_func(cur_data)
                if(try_error == 0.):
                    binned_err = np.append(binned_err, 1.)
                else:
                    binned_err = np.append(binned_err, err_calc_func(cur_data))

    return binned_time, binned_data, binned_err
