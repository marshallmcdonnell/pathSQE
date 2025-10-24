########################################################################################################
# Driver script for running pathSQE main program in conjuction with a pathSQE input file
# Author: Aiden Sable. Oct 2025.
########################################################################################################

# limit core usage to avoid extended Mantid consumption issues on shared clusters
from mantid.kernel import ConfigServiceImpl
config = ConfigServiceImpl.Instance()
config.setString('MultiThreaded.MaxCores', '8')
max_cores = config.get('MultiThreaded.MaxCores')
print("Updated Max Cores:", max_cores)

from mantid.simpleapi import *

def run_pathSQE(pathSQE_params, mde_data):
    import numpy as np
    import seekpath
    import os
    import sys
    sys.stdout = sys.__stdout__  # Reset stdout before the script exits
    import time
    import shutil


    import importlib as imp

    # Import refactored code from pathSQE_utils/ folder
    import pathSQE_utils.core as pathSQE_core
    imp.reload(pathSQE_core)

    import pathSQE_utils.plotting_and_reports as pathSQE_plotting_and_reports
    imp.reload(pathSQE_plotting_and_reports)

    import pathSQE_utils.helper as pathSQE_helper
    imp.reload(pathSQE_helper)

    import pathSQE_utils.simulations as pathSQE_simulations
    imp.reload(pathSQE_simulations)

    import pathSQE_utils.filter_functions as pathSQE_filter_functions
    imp.reload(pathSQE_filter_functions)

    import pathSQE_utils.reduce_data_to_MDE_07142023 as reduce_data_to_MDE_07142023
    imp.reload(reduce_data_to_MDE_07142023)

    import pathSQE_utils.slice_utils_07142023 as slice_utils_07142023
    imp.reload(slice_utils_07142023)

    #######################################################################################################
            ################## LOAD PARAMETERS, DATA, AND PREPARE FOR PROCESSING ###################
    #######################################################################################################

    # Find path, find spacegroup, load mde data, and create saveDir's if don't already exist, etc.
    start_time = time.time()

    # Normalize to an absolute path once and use it consistently
    pathSQE_params['output directory'] = os.path.abspath(pathSQE_params['output directory'])

    if not os.path.exists(pathSQE_params['output directory']):
        os.makedirs(pathSQE_params['output directory'])
    slices_dir = os.path.join(pathSQE_params['output directory'],'slices/')
    if not os.path.exists(slices_dir):
        os.makedirs(slices_dir)
    if not os.path.exists(os.path.join(slices_dir,'all_individual/')) and pathSQE_params['save individual slices']:
        os.makedirs(os.path.join(slices_dir,'all_individual/'))
    if not os.path.exists(os.path.join(slices_dir,'intraBZ_sims/')) and pathSQE_params['perform simulations']:
        os.makedirs(os.path.join(slices_dir,'intraBZ_sims/'))

    # Copy scripts being used to output folder
    current_folder = os.getcwd()  # Get the current working directory
    outputs_folder = pathSQE_params['output directory']  # Use the pre-defined output path
    all_individual_dir = os.path.join(slices_dir,'all_individual/')

    print(f"All files from {current_folder} being copied to {outputs_folder}")

    # Check if original pathSQE_input should be preserved
    input_file = os.path.join(outputs_folder, 'pathSQE_input.py')
    original_backup = os.path.join(outputs_folder, 'pathSQE_inputs_original.py')

    if os.path.exists(input_file) and not os.path.exists(original_backup):
        print("Preserving original pathSQE_input.py as pathSQE_inputs_original.py")
        shutil.move(input_file, original_backup)

    # Now copy all files from current_folder to outputs_folder
    for filename in os.listdir(current_folder):
        file_path = os.path.join(current_folder, filename)
        
        if os.path.isfile(file_path):
            shutil.copy(file_path, outputs_folder)

    # To write output to terminal and an output.txt at the same time, use this class
    class Tee(object):
        def __init__(self, *files):
            self.files = files

        def write(self, obj):
            for f in self.files:
                f.write(obj)
                f.flush()  # Flush to ensure the output is written immediately

        def flush(self):
            for f in self.files:
                f.flush()

    # Define the output file path
    output_file_path = os.path.join(pathSQE_params['output directory'], 'output.txt')

    # Open the file in write mode and Redirect stdout to both the terminal and the file
    output_file = open(output_file_path, 'a')
    sys.stdout = Tee(sys.stdout, output_file)


    if pathSQE_params['use seeKpath path']:
        poscar = pathSQE_helper.simple_read_poscar('POSCAR')
        user_defined_Qpoints = seekpath.get_path(structure=poscar)
    else:
        user_defined_Qpoints = pathSQE_params['user defined Qpoints']

    mtd_spacegroup = pathSQE_helper.find_matching_spacegroup(pathSQE_params)

    ###################################################
    print('Processing ', len(mde_data), ' datasets.')
    reduce_data_to_MDE_07142023.reduce_data_to_MDE(mde_data)


    if pathSQE_params['all symmetric 2d slices'] or pathSQE_params['all symmetric 1d cuts']:
        if not os.path.exists(os.path.join(pathSQE_params['output directory'],'Reports/')) and pathSQE_params['save reports']:
            os.makedirs(os.path.join(pathSQE_params['output directory'],'Reports/'))

        if isinstance(pathSQE_params['BZ to process'], list):
            BZ_list_init = [np.array(BZ) for BZ in pathSQE_params['BZ to process']] 
        else:
            soft_time = time.time()
            BZ_list_init, BZ_fractional_coverages, BZ_full_array = pathSQE_core.find_BZ_with_data(mde_data[0], pathSQE_params)
            print('BZ cov and thresh in ', time.time()-soft_time)

        BZ_list_init = np.array(BZ_list_init)
        print("Initial BZ list: ", len(BZ_list_init), BZ_list_init.tolist())

    if pathSQE_params['save individual slices']:
        sliceID = pathSQE_helper.get_next_sliceID(all_individual_dir)
    else:
        sliceID = 0


    #######################################################################################################
            ###################################### all symmetric 2d slices ###################################
    #######################################################################################################

    if pathSQE_params['all symmetric 2d slices']:
    ########## INTRA-BZ folding for each BZ specified ##########
        if not os.path.exists(os.path.join(slices_dir,'intraBZ_folded/')):
            os.makedirs(os.path.join(slices_dir,'intraBZ_folded/'))
        if not os.path.exists(os.path.join(slices_dir,'intraBZ_folded/DataAndNorms/')):
            os.makedirs(os.path.join(slices_dir,'intraBZ_folded/DataAndNorms/'))
        if not os.path.exists(os.path.join(pathSQE_params['output directory'],'folding_progess/')):
            os.makedirs(os.path.join(pathSQE_params['output directory'],'folding_progess/'))
        if not os.path.exists(os.path.join(slices_dir,'fully_folded/')):
            os.makedirs(os.path.join(slices_dir,'fully_folded/'))

        print("\nChecking for prior progress")
        BZ_list, removed_BZ = pathSQE_helper.resume_analysis(BZ_list_init, os.path.join(slices_dir,'intraBZ_folded/'), len(user_defined_Qpoints['path']), all_individual_dir)
        if len(BZ_list) == 0:
            print("\n \n ######## All BZ already processed - Skipping Reprocessing. ######## \n \n")
            if not os.path.exists(os.path.join(pathSQE_params['output directory'],'Reruns/')):
                os.makedirs(os.path.join(pathSQE_params['output directory'],'Reruns/'))
            #sys.exit(1)  # Exit with error status

        print('\nData being processed in {} BZs'.format(len(BZ_list)))
        print('\nBZs:', BZ_list.tolist())

        if pathSQE_params['perform simulations']:
            fully_folded_sims_data, fully_folded_sims_norm = pathSQE_simulations.load_previous_simulation_progress(os.path.join(slices_dir,'intraBZ_sims/'), user_defined_Qpoints, removed_BZ)
            print('loaded data and norm for sim lens ', len(fully_folded_sims_data), len(fully_folded_sims_norm))

        fold_timer = time.time()
        num_total_slices = pathSQE_helper.count_total_bzfold_slices(pathSQE_params, BZ_list, mtd_spacegroup)
        print(f"\n\n Total expected slices: {num_total_slices}\n\n")
        num_processed_slices = 0
        last_update_step = 0

        for B, BZ_offset in enumerate(BZ_list):
            last_BZ_fold_dsl = []
            dsl_fold = []
            all_sliceInfo_for_BZ = []
            all_sims_in_BZ = [[] for i in range(len(user_defined_Qpoints['path']))]
            for i in range(len(user_defined_Qpoints['path'])):
                path_seg = user_defined_Qpoints['path'][i]
                print('Processing BZ {}, path segments from {} to {}'.format(BZ_offset, path_seg[0], path_seg[1]))

                pt1_init = np.array(user_defined_Qpoints['point_coords'][path_seg[0]])
                pt2_init = np.array(user_defined_Qpoints['point_coords'][path_seg[1]])

                pt1_array, pt2_array = pathSQE_core.generate_unique_paths(mtd_spacegroup, pt1_init, pt2_init, pathSQE_params)
                print("\n{} symmetrically equivalent path segments".format(pt1_array.shape[0]))
                print(np.hstack((pt1_array,pt2_array)))

                all_sliceInfo_for_pathSeg = []
                good_slice_index = None
                for j in range(pt1_array.shape[0]):
                    q_dims_and_bins = pathSQE_core.choose_dims_and_bins(pathSQE_params, pt1_array[j], pt2_array[j],  pathSQE_params['u_vec'], pathSQE_params['v_vec'], BZ_offset=BZ_offset)
                    slice_desc = pathSQE_helper.make_slice_desc(pathSQE_params, q_dims_and_bins, pt1_array[j], pt2_array[j], path_seg)
                                    
                    # make slice and evaluate quality based on filters
                    slice_utils_07142023.make_slice(mde_data[0], slice_desc, ASCII_slice_folder='', MD_slice_folder='')
                    if pathSQE_params['save individual slices']:
                        SaveMD(mtd['_data'], Filename=os.path.join(all_individual_dir,'BZ{}_{}2{}_{}to{}_ID{}_data.nxs'.format(BZ_offset, path_seg[0], path_seg[1], pt1_array[j]+BZ_offset, pt2_array[j]+BZ_offset, sliceID)).replace("\\","/"), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
                        SaveMD(mtd['_norm'], Filename=os.path.join(all_individual_dir,'BZ{}_{}2{}_{}to{}_ID{}_norm.nxs'.format(BZ_offset, path_seg[0], path_seg[1], pt1_array[j]+BZ_offset, pt2_array[j]+BZ_offset, sliceID)).replace("\\","/"), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
                    slice_good = pathSQE_filter_functions.evaluate_slice_quality(slice_name=slice_desc['Name'], filters=pathSQE_params['slice filter functions'], path_seg=path_seg)
                    all_sliceInfo_for_pathSeg.append((slice_desc['Name'],sliceID,slice_good))

                    if not slice_good:
                        print('artifact found - returned False for ', slice_desc['Name'], sliceID)

                    # optionally do simulation workflow for analogous slice
                    if pathSQE_params['perform simulations']:
                        Qpoints = pathSQE_simulations.construct_sim_Qpts(pt1_array[j], pt2_array[j], q_dims_and_bins[0], pathSQE_params['qdim0 step size'], pathSQE_params['primitive to mantid transformation matrix'], BZ_offset)
                        #print(len(Qpoints))
                        sim_SQE_output = pathSQE_simulations.sim_SQE(pathSQE_params, Qpoints, Temperature=pathSQE_params['T and Ei conditions'][0][0])
                        #print(sim_SQE_output.shape)
                        if len(pathSQE_params['resolution blurring']) == 3:
                            BinnedSQE = pathSQE_simulations.SQE_to_2d_spectrum_advancedRes(pathSQE_params, sim_SQE_output)
                        else:
                            BinnedSQE = pathSQE_simulations.SQE_to_2d_spectrum(pathSQE_params, sim_SQE_output)
                        #print(BinnedSQE.shape)

                        if pathSQE_params['use experimental coverage mask']:
                            # Get the experimental signal array
                            slice_data = np.squeeze(mtd[slice_desc['Name']].getSignalArray())

                            # Ensure shapes match before masking
                            if slice_data.shape == BinnedSQE.shape:
                                BinnedSQE[np.isnan(slice_data)] = np.nan
                            else:
                                print(f"Warning: Shape mismatch between experimental data {slice_data.shape} and simulated data {BinnedSQE.shape}. Masking not applied.")

                        all_sims_in_BZ[i].append(BinnedSQE)
                        if pathSQE_params['save individual slices']:
                            np.save(os.path.join(all_individual_dir,'BZ{}_{}2{}_{}to{}_ID{}_sim.npy'.format(BZ_offset, path_seg[0], path_seg[1], pt1_array[j]+BZ_offset, pt2_array[j]+BZ_offset, sliceID)), BinnedSQE)
                        #print('num nans ', np.isnan(BinnedSQE).sum())  # Count NaN values

                    # only use good slices in the folding
                    if slice_good:
                        slice_desc['good_slice'] = True
                        if good_slice_index is None:
                            # First good slice, clone the data
                            data = mtd['_data'].clone()
                            norm = mtd['_norm'].clone()
                            dsl_fold.append(slice_desc.copy())
                            dsl_fold[i]['Name'] = '{}2{}_BZ{}_final'.format(path_seg[0], path_seg[1], BZ_offset)
                            good_slice_index = j  # Save the index of the first good slice
                            
                            # save slice desc of each segment for plotting if this is the recently processed BZ
                            last_BZ_fold_dsl.append(slice_desc.copy())
                            
                        else:
                            # Subsequent good slices, combine data
                            data, norm = pathSQE_core.combine_data_within_bz(data, norm, mtd['_data'], mtd['_norm'])

                            if not slice_good:
                                print('somehow folded artifact slice ', slice_desc['Name'], sliceID)

                    else:
                        slice_desc['good_slice'] = False

                    sliceID += 1

                    # Check if we should print an update
                    num_processed_slices += 1
                    should_print, last_update_step = pathSQE_helper.should_print_update(num_processed_slices, num_total_slices, last_update_step, update_percent=1)
                    if should_print:
                        pathSQE_helper.print_progress(num_processed_slices, num_total_slices, fold_timer, "BZ Folding")

                # Save the data for the last good slice (if any)
                if good_slice_index is not None:
                    saveDir = os.path.join(slices_dir,'intraBZ_folded/')
                    SaveMD(data, Filename=os.path.join(saveDir,'DataAndNorms/{}2{}_BZ{}_data.nxs'.format(path_seg[0], path_seg[1], BZ_offset)).replace("\\","/"), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
                    SaveMD(norm, Filename=os.path.join(saveDir,'DataAndNorms/{}2{}_BZ{}_norm.nxs'.format(path_seg[0], path_seg[1], BZ_offset)).replace("\\","/"), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)

                    output_ws = data / norm
                    name = '{}2{}_BZ{}_final.nxs'.format(path_seg[0], path_seg[1], BZ_offset)
                    SaveMD(output_ws, Filename=os.path.join(saveDir,name).replace("\\","/"), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
                    LoadMD(Filename=os.path.join(saveDir,name), OutputWorkspace=name, LoadHistory=False)
                    
                all_sliceInfo_for_BZ.append(all_sliceInfo_for_pathSeg)

            if pathSQE_params['perform simulations']:
                folded_sim_results = pathSQE_simulations.fold_symmetric_simulations(all_sims_in_BZ)

                intraBZ_sims = os.path.join(slices_dir, 'intraBZ_sims/')
                np.savez(os.path.join(intraBZ_sims, 'BZ{}_foldedSim_data.npz'.format(BZ_offset)), *[result[0] for result in folded_sim_results])
                np.savez(os.path.join(intraBZ_sims, 'BZ{}_foldedSim_norm.npz'.format(BZ_offset)), *[result[1] for result in folded_sim_results])
                np.savez(os.path.join(intraBZ_sims, 'BZ{}_foldedSim_sim.npz'.format(BZ_offset)), *[result[2] for result in folded_sim_results])
            
                if len(fully_folded_sims_data)==0:
                    fully_folded_sims_data = [results[0] for results in folded_sim_results]
                    fully_folded_sims_norm = [results[1] for results in folded_sim_results]
                else:
                    for tic in range(len(fully_folded_sims_data)):
                        #print('seg ', tic)
                        fully_folded_sims_data[tic] += folded_sim_results[tic][0]
                        fully_folded_sims_norm[tic] += folded_sim_results[tic][1]
            
            if pathSQE_params['save reports']:
                print('\nGenerating Report for BZ {}'.format(BZ_offset))
                if pathSQE_params['perform simulations']:
                    pathSQE_plotting_and_reports.generate_BZ_Report_with_Sims(pathSQE_params, BZ_offset, all_sliceInfo_for_BZ, dsl_fold, folded_sim_results, all_sims_in_BZ)
                else:
                    pathSQE_plotting_and_reports.generate_BZ_Report_morePages(pathSQE_params, BZ_offset, all_sliceInfo_for_BZ, dsl_fold)

            if not isinstance(pathSQE_params['BZ to process'], list):
                pathSQE_plotting_and_reports.update_BZ_folding_plot(BZ_full_array, BZ_offset, pathSQE_params)


            ########### INTER-BZ folding exp data over all BZ specified ###########
            if len(BZ_list_init)>1: 
                BZ_fold_list = np.vstack((removed_BZ, BZ_list[:B]))
                seg_names = pathSQE_core.fold_BZs(BZ_fold_list, BZ_offset, pathSQE_params, user_defined_Qpoints)

                # Adaptive colorbars, plotting, saving outputs
                pathMin, pathMax = pathSQE_helper.adaptive_colorbars(seg_names, pathSQE_params, percentile_min=10, percentile_max=95)
                folded_array = [np.squeeze(mtd[seg].getSignalArray()) for seg in seg_names]
                np.save(os.path.join(pathSQE_params['output directory'],'folding_progess/','folded_path_plot_{}.npy'.format(B+len(removed_BZ))), np.vstack(folded_array))
                fig = pathSQE_plotting_and_reports.plot_along_path_foldedBZ(foldedSegNames=seg_names, dsl_fold=last_BZ_fold_dsl, Ei=pathSQE_params['T and Ei conditions'][0][1], vmi=pathMin, vma=pathMax, cma='viridis')
                fig.savefig(os.path.join(pathSQE_params['output directory'],'folding_progess/','folded_path_plot_{}.png'.format(B+len(removed_BZ))))
                fig2 = pathSQE_plotting_and_reports.plot_SED_along_path_foldedBZ(foldedSegNames=seg_names, dsl_fold=last_BZ_fold_dsl, pathSQE_params=pathSQE_params)
                fig2.savefig(os.path.join(pathSQE_params['output directory'],'folding_progess/','folded_SED_path_plot_{}.png'.format(B+len(removed_BZ))))

                if pathSQE_params['perform simulations']: 
                    folded_sim_final = []
                    for ps in range(len(fully_folded_sims_data)):
                        final_norm = fully_folded_sims_norm[ps].copy()
                        final_norm[final_norm == 0] = np.nan
                        folded_sim_final.append( fully_folded_sims_data[ps] / final_norm )

                    simMax = 0
                    simMin = 1e4
                    # Find min/max percentiles across all slices
                    for path_sim in folded_sim_final:
                        mask = ~np.isnan(path_sim) & (path_sim > 0)
                        if np.any(mask):
                            segMin = np.percentile(path_sim[mask],60)
                            segMax = np.max(path_sim[mask])
                            simMax = max(pathMax, segMax)
                            simMin = min(pathMin, segMin)

                    fig3 = pathSQE_plotting_and_reports.plot_along_path_foldedBZ_sim(folded_sim_final,last_BZ_fold_dsl,pathSQE_params['T and Ei conditions'][0][1],simMin,simMax,'viridis')
                    fig3.savefig(os.path.join(pathSQE_params['output directory'],'folding_progess/','folded_sim_path_plot_{}.png'.format(B+len(removed_BZ))))
                    np.savez(os.path.join(pathSQE_params['output directory'], 'folding_progess/', 'folded_sim_path_plot_{}.npz'.format(B+len(removed_BZ))), *folded_sim_final)

                print('\nBZ {} done. Total time elapsed {}'.format(BZ_offset, pathSQE_helper.format_time(time.time()-start_time)))

        print('\nFolded {} BZ in {}'.format(len(BZ_list_init), pathSQE_helper.format_time(time.time()-start_time)))

        # if all data was already processed with mantid, skip the above and just regenerate outputs in Reruns folder
        if len(BZ_list) == 0:
            print("Folding pre-processed BZs: \n", removed_BZ)
            seg_names = pathSQE_core.fold_BZs(removed_BZ, removed_BZ[0], pathSQE_params, user_defined_Qpoints, rerun=True)

            # Adaptive colorbars, plotting, saving outputs
            pathMin, pathMax = pathSQE_helper.adaptive_colorbars(seg_names, pathSQE_params, percentile_min=10, percentile_max=95)
            folded_array = [np.squeeze(mtd[seg].getSignalArray()) for seg in seg_names]
            np.save(os.path.join(pathSQE_params['output directory'],'Reruns/','folded_path_plot.npy'), np.vstack(folded_array))
            fig = pathSQE_plotting_and_reports.plot_along_path_rerun(seg_names, user_defined_Qpoints, pathMin, pathMax)
            fig.savefig(os.path.join(pathSQE_params['output directory'],'Reruns/','folded_path_plot.png'))
            #fig2 = pathSQE_plotting_and_reports.plot_SED_along_path_foldedBZ(foldedSegNames=seg_names, dsl_fold=last_BZ_fold_dsl, pathSQE_params=pathSQE_params)
            #fig2.savefig(os.path.join(pathSQE_params['output directory'],'Reruns/','folded_SED_path_plot.png'))



    #######################################################################################################
            #################################### all symmetric 1d cuts ##################################
    #######################################################################################################

    if pathSQE_params['all symmetric 1d cuts']:
        sym_points_timer = time.time()
        num_total_slices = pathSQE_helper.count_total_symPt_slices(pathSQE_params, BZ_list_init, mtd_spacegroup, len(mde_data))
        print(f"Total expected slices: {num_total_slices}")
        num_processed_slices = 0
        last_update_step = 0

        for i in range(len(user_defined_Qpoints['1d_points'])):
            symPt_slice_names_allTemp = [[] for _ in range(len(mde_data))]
            symPt_slice_arrays = []
            
            # find unique sym points (arbitrary BZ)
            symPoint = user_defined_Qpoints['1d_points'][i]
            symPt_init = np.array(user_defined_Qpoints['point_coords'][symPoint])
            symPts_array = pathSQE_core.generate_unique_symPts(mtd_spacegroup, symPoint, symPt_init)
            print('\nPoint', symPoint)
            print('Equivalent points in a single BZ: ', symPts_array)

            # find all unique sym points in BZ coverage (avoid repeating same point in different BZ)
            allSymPts_array = pathSQE_core.generate_unique_symPts_inAllBZ(BZ_list_init, symPts_array)
            np.save(os.path.join(pathSQE_params['output directory'], 'probed_{}_pts.npy'.format(symPoint)), allSymPts_array)

            for j in range(allSymPts_array.shape[0]):
                for k in range(len(mde_data)):
                    if k == 0:
                        slice_desc = pathSQE_helper.make_slice_desc_1DSymPoints(pathSQE_params, symPoint, allSymPts_array[j])
                        slice_desc_temp2og = slice_desc.copy()
                        slice_desc['Name'] = slice_desc['Name']+'_{}K'.format(mde_data[k]['SampleLogVariables']['Temperature'])
                        slice_utils_07142023.make_slice(mde_data[k], slice_desc, ASCII_slice_folder='', MD_slice_folder='') 
                        
                        if np.count_nonzero((mtd[slice_desc['Name']].getSignalArray() != 0) & (~np.isnan(mtd[slice_desc['Name']].getSignalArray()))) > 5:
                            symPt_slice_arrays.append(allSymPts_array[j])
                            symPt_slice_names_allTemp[k].append(slice_desc['Name'])
                            if pathSQE_params['save individual slices']:
                                SaveMD(mtd[slice_desc['Name']], Filename=os.path.join(slices_dir,'all_individual/',slice_desc['Name']+'.nxs').replace("\\","/"), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
                
                    else:
                        # for other temps
                        slice_desc_temp2 = slice_desc_temp2og.copy()
                        slice_desc_temp2['Name'] = slice_desc_temp2['Name']+'_{}K'.format(mde_data[k]['SampleLogVariables']['Temperature'])
                        slice_utils_07142023.make_slice(mde_data[k], slice_desc_temp2, ASCII_slice_folder='', MD_slice_folder='') 
                    
                        if np.count_nonzero((mtd[slice_desc['Name']].getSignalArray() != 0) & (~np.isnan(mtd[slice_desc['Name']].getSignalArray()))) > 5:     
                            symPt_slice_names_allTemp[k].append(slice_desc_temp2['Name'])
                            if pathSQE_params['save individual slices']:
                                SaveMD(mtd[slice_desc_temp2['Name']], Filename=os.path.join(slices_dir,'all_individual/',slice_desc_temp2['Name']+'.nxs').replace("\\","/"), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
                            
                    # Check if we should print an update
                    num_processed_slices += 1
                    should_print, last_update_step = pathSQE_helper.should_print_update(num_processed_slices, num_total_slices, last_update_step, update_percent=1)
                    if should_print:
                        pathSQE_helper.print_progress(num_processed_slices, num_total_slices, sym_points_timer, "1D Sym Point")


            #print('sym point slice names: ', symPt_slice_names_allTemp)
            if pathSQE_params['save reports']:
                print('Generating Report for {} point'.format(symPoint))
                pathSQE_plotting_and_reports.generate_BZ_Report_1DSymPoints(pathSQE_params, symPoint, symPt_slice_arrays, symPt_slice_names_allTemp, pathSQE_params['u_vec'], pathSQE_params['v_vec'])


    #######################################################################################################
            ################################ SIMPLE 2D PATH ####################################
    #######################################################################################################

    if not pathSQE_params['all symmetric 2d slices'] and user_defined_Qpoints['path'] and len(user_defined_Qpoints['path']) > 0:

        dsl=[]
        seg_names = []
        all_sims = []
        for path_seg in user_defined_Qpoints['path']:
            pt1 = np.array(user_defined_Qpoints['point_coords'][path_seg[0]])
            pt2 = np.array(user_defined_Qpoints['point_coords'][path_seg[1]])
            q_dims_and_bins = pathSQE_core.choose_dims_and_bins(pathSQE_params, pt1, pt2, pathSQE_params['u_vec'], pathSQE_params['v_vec'])  
            slice_desc = pathSQE_helper.make_slice_desc(pathSQE_params, q_dims_and_bins, pt1, pt2, path_seg)
            dsl.append(slice_desc)
            seg_names.append(slice_desc['Name'])
            slice_utils_07142023.make_slice(mde_data[0], slice_desc, ASCII_slice_folder='', MD_slice_folder='')
            if pathSQE_params['save individual slices']:
                SaveMD(slice_desc['Name'], Filename=os.path.join(slices_dir,'all_individual/{}_to_{}.nxs'.format(path_seg[0], path_seg[1])).replace("\\","/"), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)

            # optionally do simulation workflow for analogous slice
            if pathSQE_params['perform simulations']:
                Qpoints = pathSQE_simulations.construct_sim_Qpts(pt1, pt2, q_dims_and_bins[0], pathSQE_params['qdim0 step size'], pathSQE_params['primitive to mantid transformation matrix'])
                sim_SQE_output = pathSQE_simulations.sim_SQE(pathSQE_params, Qpoints, Temperature=pathSQE_params['T and Ei conditions'][0][0])
                if len(pathSQE_params['resolution blurring']) == 3:
                    BinnedSQE = pathSQE_simulations.SQE_to_2d_spectrum_advancedRes(pathSQE_params, sim_SQE_output)
                else:
                    BinnedSQE = pathSQE_simulations.SQE_to_2d_spectrum(pathSQE_params, sim_SQE_output)

                if pathSQE_params['use experimental coverage mask']:
                    # Get the experimental signal array
                    slice_data = np.squeeze(mtd[slice_desc['Name']].getSignalArray())

                    # Ensure shapes match before masking
                    if slice_data.shape == BinnedSQE.shape:
                        BinnedSQE[np.isnan(slice_data)] = np.nan
                    else:
                        print(f"Warning: Shape mismatch between experimental data {slice_data.shape} and simulated data {BinnedSQE.shape}. Masking not applied.")

                all_sims.append(BinnedSQE)

        # Adaptive colorbars    
        pathMin, pathMax = pathSQE_helper.adaptive_colorbars(seg_names, pathSQE_params, percentile_min=10, percentile_max=95)
                
        fig = pathSQE_plotting_and_reports.plot_along_path_foldedBZ(foldedSegNames=seg_names, dsl_fold=dsl, Ei=pathSQE_params['T and Ei conditions'][0][1], vmi=pathMin, vma=pathMax, cma='viridis')
        fig.savefig(os.path.join(pathSQE_params['output directory'],'exp_path_plot.png'))
        combined_array = [np.squeeze(mtd[seg].getSignalArray()) for seg in seg_names] 
        np.savez(os.path.join(pathSQE_params['output directory'],'exp_path_plot.npz'), *combined_array)

        if pathSQE_params['perform simulations']:
            simMax = 0
            simMin = 1e4
            # Find min/max percentiles across all slices
            for path_sim in all_sims:
                mask = ~np.isnan(path_sim) & (path_sim > 0)
                if np.any(mask):
                    segMin = np.percentile(path_sim[mask],60)
                    segMax = np.max(path_sim[mask])
                    simMax = max(pathMax, segMax)
                    simMin = min(pathMin, segMin)

            fig3 = pathSQE_plotting_and_reports.plot_along_path_foldedBZ_sim(all_sims,dsl,pathSQE_params['T and Ei conditions'][0][1],simMin,simMax,'viridis')
            fig3.savefig(os.path.join(pathSQE_params['output directory'],'sim_path_plot.png'))
            np.savez(os.path.join(pathSQE_params['output directory'],'sim_path_plot.npz'), *all_sims)


    #######################################################################################################
            ################################ SIMPLE 1D POINTS ####################################
    #######################################################################################################

    if not pathSQE_params['all symmetric 1d cuts'] and user_defined_Qpoints['1d_points'] and len(user_defined_Qpoints['1d_points']) > 0:

        for i in range(len(user_defined_Qpoints['1d_points'])):
            slice_names = []
            
            # find unique sym points (arbitrary BZ)
            symPoint = user_defined_Qpoints['1d_points'][i]
            symPt_init = np.array(user_defined_Qpoints['point_coords'][symPoint])
            print('\nPoint', symPoint, symPt_init)

            for k in range(len(mde_data)):
                if k == 0:
                    slice_desc = pathSQE_helper.make_slice_desc_1DSymPoints(pathSQE_params, symPoint, symPt_init)
                    slice_desc_temp2og = slice_desc.copy()
                    slice_desc['Name'] = slice_desc['Name']+'_{}K'.format(mde_data[k]['SampleLogVariables']['Temperature'])
                    slice_names.append(slice_desc['Name'])
                    slice_utils_07142023.make_slice(mde_data[k], slice_desc, ASCII_slice_folder='', MD_slice_folder='') 
                    
                    if np.count_nonzero((mtd[slice_desc['Name']].getSignalArray() != 0) & (~np.isnan(mtd[slice_desc['Name']].getSignalArray()))) > 5 and pathSQE_params['save individual slices']:
                        SaveMD(slice_desc['Name'], Filename=os.path.join(slices_dir,'all_individual/',slice_desc['Name']+'.nxs').replace("\\","/"), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
            
                else:
                    # for other temps
                    slice_desc_temp2 = slice_desc_temp2og.copy()
                    slice_desc_temp2['Name'] = slice_desc_temp2['Name']+'_{}K'.format(mde_data[k]['SampleLogVariables']['Temperature'])
                    slice_names.append(slice_desc_temp2['Name'])
                    slice_utils_07142023.make_slice(mde_data[k], slice_desc_temp2, ASCII_slice_folder='', MD_slice_folder='') 
                
                    if np.count_nonzero((mtd[slice_desc['Name']].getSignalArray() != 0) & (~np.isnan(mtd[slice_desc['Name']].getSignalArray()))) > 5 and pathSQE_params['save individual slices']:     
                        SaveMD(slice_desc_temp2['Name'], Filename=os.path.join(slices_dir,'all_individual/',slice_desc_temp2['Name']+'.nxs').replace("\\","/"), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
            
            pathSQE_plotting_and_reports.plot_simple_1d_pt(pathSQE_params, symPoint, symPt_init, slice_names)




    # Remember to close the output file when done
    time.sleep(5)
    output_file.close()
    sys.stdout = sys.__stdout__  # Reset stdout before the script exits



if __name__ == "__main__":

    import importlib as imp
    import pathSQE_input
    imp.reload(pathSQE_input)
    import define_data
    imp.reload(define_data)

    pathSQE_params = pathSQE_input.define_pathSQE_params()
    #print(pathSQE_params)
    define_data_params = define_data.define_data_set(pathSQE_params['T and Ei conditions'])

    run_pathSQE(pathSQE_params, define_data_params)
