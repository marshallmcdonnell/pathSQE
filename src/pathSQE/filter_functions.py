import numpy as np
from mantid.simpleapi import *
import matplotlib.pyplot as plt


def evaluate_slice_quality(slice_name, filters, path_seg):
    # if not using filters, just say all slices are good
    if len(filters) == 0:
        return True
    
    filter_evals = []
    # filter based on fractional coverage of slice
    if 'coverage' in filters:
        slice_data = mtd[slice_name].getSignalArray()
        fractional_coverage = np.count_nonzero(~np.isnan(slice_data)) / slice_data.size
        print('fc: ', fractional_coverage)
        
        if fractional_coverage > 0.54:
            filter_evals.append(True)
        else:
            filter_evals.append(False)
    

    # UNDER DEV - filter looking at distribution of pixel intensities (contrast/sharpness) above elastic line (assumed bottom 3 rows)
    if 'contrast' in filters:
        slice_data = mtd[slice_name].getSignalArray()
        slice_data_withoutElasticLine = slice_data[:,0,0,3:]
        non_nan_values = slice_data_withoutElasticLine[~np.isnan(slice_data_withoutElasticLine)]
        #print(np.min(non_nan_values), np.max(non_nan_values), np.mean(non_nan_values), non_nan_values.shape)
        # Create histogram with log scale x-axis
        plt.figure()
        # Specify logarithmically spaced bin edges
        bins = np.logspace(-4.5, -1.5, 50)

        # Create histogram with log scale on both axes
        plt.hist(non_nan_values, bins=bins, log=True, edgecolor='black')
        plt.xlabel('Values (log scale)')
        plt.ylabel('Frequency (log scale)')
        plt.title('Histogram of Non-NaN Values (log-log scale)')
        plt.grid(True)
        plt.xscale('log')  # Set log scale on x-axis
        #plt.xlim([1e-4, np.max(non_nan_values)])

        # Save the plot as a PNG file
        plt.savefig(os.path.join('/SNS/ARCS/IPTS-5307/shared/Aiden/comprehensive/pathSQE_outputs/withFilters/',slice_name,'.png'))
        
        if 1 > 0.54:
            filter_evals.append(True)
        else:
            filter_evals.append(False)    


    # Filter to get rid of M point artifacts in 40meV FeSi dataset (IPTS-5307)
    if 'M_bragg_tails' in filters: 
        if path_seg == ('M', 'Gamma'):
            slice_data = mtd[slice_name].getSignalArray()
            slice_data = np.squeeze(slice_data)

            # Mask where bad signal is (approx region for M→Γ)
            shape = slice_data.shape
            feat_mask = np.zeros(shape)
            feat_mask[:13, 13:27] = 1  # E 3–7 meV and q 0–2/3 M→Γ

            try:
                # Region to apply percentile threshold
                E_ind_high = int(np.ceil(3 / 0.25))  
                slice_sub = slice_data[:, E_ind_high:]

                # Only consider non-NaN, non-zero values
                valid_mask = ~np.isnan(slice_sub) & (slice_sub != 0)
                if not np.any(valid_mask):
                    # If no valid data, don't reject the slice
                    print('no valid data - return True')
                    filter_evals.append(True)
                else:
                    segMax = np.percentile(slice_sub[valid_mask], 98)
                    perc_mask = slice_data >= segMax

                    # Apply both masks
                    if np.sum(feat_mask * perc_mask) > 2:
                        print('artifact found - return False')
                        filter_evals.append(False)
                    else:
                        print('no artifact found - return True')
                        filter_evals.append(True)
            except Exception as e:
                print(f"[Warning] M_bragg_tails filter skipped due to error: {e}")
                filter_evals.append(True)  # Assume safe if error
        else:
            filter_evals.append(True)



    # filter to get rid of M point artifacts in 70meV FeSi dataset (IPTS-21211)
    if 'M_bragg_tails_70meV' in filters: 
        if path_seg==('M','Gamma'):
            slice_data = mtd[slice_name].getSignalArray()
            slice_data = np.squeeze(slice_data)
            
            # mask where bad signal is
            shape = slice_data.shape
            feat_mask = np.zeros(shape)
            feat_mask[:10,10:17] = 1  # E 5.5-9 ish and q 0-0.2 on M to Gamma
            
            # mask data based on percentile
            E_ind_high = int(np.ceil(4/0.5))  
            slice = slice_data[:, E_ind_high:]
            mask = ~np.isnan(slice) & (slice != 0)
            segMax = np.percentile(slice[mask],95)
            perc_mask = slice_data >= segMax
                    
            if np.sum(feat_mask*perc_mask)>5:
                filter_evals.append(False)
            else:
                filter_evals.append(True)
        elif path_seg==('X','M'):
            slice_data = mtd[slice_name].getSignalArray()
            slice_data = np.squeeze(slice_data)
            
            # mask where bad signal is
            shape = slice_data.shape
            feat_mask = np.zeros(shape)
            feat_mask[12:,9:19] = 1  # E 5-10 ish and q 0.3-0.5 on X to M
            
            # mask data based on percentile
            E_ind_high = int(np.ceil(4/0.5))  
            slice = slice_data[:, E_ind_high:]
            mask = ~np.isnan(slice) & (slice != 0)
            segMax = np.percentile(slice[mask],95)
            perc_mask = slice_data >= segMax
                    
            if np.sum(feat_mask*perc_mask)>5:
                filter_evals.append(False)
            else:
                filter_evals.append(True)
        else:
            filter_evals.append(True)           

    


    # Filter to get rid of M point artifacts in 40meV FeSi dataset (IPTS-5307)
    if 'MnO_bragg_tails' in filters: 
        if path_seg == ('K', 'Gamma'):
            slice_data = mtd[slice_name].getSignalArray()
            slice_data = np.squeeze(slice_data)

            # Mask where bad signal is (approx region for M→Γ)
            shape = slice_data.shape
            feat_mask = np.zeros(shape)
            feat_mask[9:22, 7:18] = 1  # E 4–9 meV and q 1/3 – 3/4 K→Γ

            try:
                # Region to apply percentile threshold
                E_ind_high = int(np.ceil(4 / 0.5))  
                slice_sub = slice_data[:, E_ind_high:]

                # Only consider non-NaN, non-zero values
                valid_mask = ~np.isnan(slice_sub) & (slice_sub != 0)
                if not np.any(valid_mask):
                    # If no valid data, don't reject the slice
                    print('no valid data - return True')
                    filter_evals.append(True)
                else:
                    segMax = np.percentile(slice_sub[valid_mask], 98)
                    perc_mask = slice_data >= segMax

                    # Apply both masks
                    if np.sum(feat_mask * perc_mask) > 2:
                        print('artifact found - return False')
                        filter_evals.append(False)
                    else:
                        print('no artifact found - return True')
                        filter_evals.append(True)
            except Exception as e:
                print(f"[Warning] M_bragg_tails filter skipped due to error: {e}")
                filter_evals.append(True)  # Assume safe if error
        else:
            filter_evals.append(True)





    # more filters customized...
    
    
    # final check to see if slice passes all filters
    slice_good = np.all(filter_evals)
    return slice_good

