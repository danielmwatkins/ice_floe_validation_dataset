"""
This script links floes detected in both the Aqua and the Terra images. The goal is to have a set of
manually validated and corrected pairings between ice floes, such that we can examine the variance of the 
measures under expert floe edge detection. This serves as a baseline for the confidence in floe measurements.

The matching procedure includes the following
1. Identification of floes with a high degree of overlap (IOU > 0.5)
2. Overlaying imagery including potential matches for manual inspection
3. Manual addition of selected floes with low overlap
4. 

"""

import skimage as sk
import proplot as pplt
import numpy as np
import os
import pandas as pd
from skimage.color import rgba2rgb

modis_loc = "../data/modis/truecolor/"
labeled_loc = "../data/validation_images/labeled_floes/"
labeled_images = [f for f in os.listdir(labeled_loc) if "tiff" in f]

# check that both cases are there 
aqua_cases = list(filter(lambda x: 'aqua' in x, labeled_images))
cases = list(filter(lambda x: x.replace('aqua', 'terra') in labeled_images, aqua_cases))
cases.sort()

# get info on number of floes per image
case_info = pd.DataFrame(columns=['n_aqua', 'n_terra', 'mean_area_aqua', 'mean_area_terra', 'area_overlap'],
                         index=[c.split('-')[0] for c in cases])
for row, case in zip(case_info.index, cases):
    fpath_lb_aqua = labeled_loc + case
    fpath_lb_terra = labeled_loc + case.replace('aqua', 'terra')

    # Convert float to int for label images
    # Also drop the unused channel layer
    lb_aqua = sk.io.imread(fpath_lb_aqua).astype(int)
    lb_terra = sk.io.imread(fpath_lb_terra).astype(int)
    if len(lb_aqua.shape) != 2:
        lb_aqua = lb_aqua[:,:,0]
    if len(lb_terra.shape) != 2:
        lb_terra = lb_terra[:,:,0]
    area_aqua = sk.measure.regionprops_table(lb_aqua, properties=['area'])['area']
    area_terra = sk.measure.regionprops_table(lb_terra, properties=['area'])['area']

    area_overlap = np.sum((lb_aqua > 0) & (lb_terra > 0))
    
    case_info.loc[row, :] = [len(area_aqua), len(area_terra), np.mean(area_aqua), np.mean(area_terra), area_overlap]
    if lb_aqua.shape[0] != 400:
        print(row, 'aqua', lb_aqua.shape)
    if lb_terra.shape[0] != 400:
        print(row, 'terra', lb_terra.shape)
    
case_info = case_info.fillna(0) 
case_info['pair_floes'] = True
case_info.loc[case_info[['n_aqua', 'n_terra']].min(axis=1) == 0, 'pair_floes'] = False
print('Number of candidate cases: ', case_info['pair_floes'].sum(), '/', len(cases))

######## Function to load images for creating the tables and plotting ############
def load_case(case, modis_loc="../data/modis/truecolor/", labeled_loc="../data/validation_images/labeled_floes/"):
    """
    Load the labeled and truecolor images for the case, grab the region props tables, and make a first guess of matches.
    """
    cn, rg, dt, sat, imtype = case.split('-')
    fpath_tc_aqua = modis_loc + '-'.join([cn, rg, '100km', dt]) + '.'.join(['.aqua', 'truecolor', '250m', 'tiff'])
    fpath_tc_terra = modis_loc + '-'.join([cn, rg, '100km', dt]) + '.'.join(['.terra', 'truecolor', '250m', 'tiff'])
    fpath_lb_aqua = labeled_loc + case
    fpath_lb_terra = labeled_loc + case.replace('aqua', 'terra')

    # Convert to RGB -- we don't use the opacity layer
    tc_aqua = rgba2rgb(sk.io.imread(fpath_tc_aqua))
    tc_terra = rgba2rgb(sk.io.imread(fpath_tc_terra))
    
    # Convert float to int for label images
    # Also drop the unused channel layer
    lb_aqua = sk.io.imread(fpath_lb_aqua).astype(int)
    lb_terra = sk.io.imread(fpath_lb_terra).astype(int)
    
    # A couple cases have odd dimensions. Most are 400 by 400 by 3 but some lack the 3rd dim.
    if len(lb_aqua.shape) != 2:
        lb_aqua = lb_aqua[:,:,0]
    if len(lb_terra.shape) != 2:
        lb_terra = lb_terra[:,:,0]
    
    # Load land masks and landfast ice masks
    # Also drop the unused channel layer
    lf_aqua = sk.io.imread(fpath_lb_aqua.replace('labeled_floes', 'binary_landfast').replace('tiff', 'png')).astype(int)
    lf_terra = sk.io.imread(fpath_lb_terra.replace('labeled_floes', 'binary_landfast').replace('tiff', 'png')).astype(int)
    if len(lf_aqua.shape) != 2:
        lf_aqua = lf_aqua[:,:,0]
    if len(lf_terra.shape) != 2:
        lf_terra = lf_terra[:,:,0]

    lm_aqua = sk.io.imread(fpath_lb_aqua.replace('labeled_floes', 'binary_landfast').replace('tiff', 'png')).astype(int)
    lm_terra = sk.io.imread(fpath_lb_terra.replace('labeled_floes', 'binary_landfast').replace('tiff', 'png')).astype(int)
    
    if len(lm_aqua.shape) != 2:
        lm_aqua = lm_aqua[:,:,0]
    if len(lm_terra.shape) != 2:
        lm_terra = lm_terra[:,:,0]

    ####### Find regions of overlap ########
    regions_aqua = pd.DataFrame(sk.measure.regionprops_table(
        lb_aqua, properties=['label', 'area', 'convex_area', 'centroid', 'perimeter', 'axis_major_length', 'axis_minor_length'])).set_index('label')
    regions_terra = pd.DataFrame(sk.measure.regionprops_table(
        lb_terra, properties=['label', 'area', 'convex_area', 'centroid', 'perimeter', 'axis_major_length', 'axis_minor_length'])).set_index('label')

    # Labels where there is some overlap
    aqua_labels = np.ravel(lb_aqua)
    terra_labels = np.ravel(lb_terra)
    labels = np.unique(aqua_labels[(aqua_labels > 0) & (terra_labels > 0)])

    matches = []
    for label in labels:
        idx_aqua = aqua_labels == label
        candidates = np.unique(terra_labels[idx_aqua & (terra_labels > 0)])
        # print(candidates)
        for c in candidates:
            idx_terra = terra_labels == c
            iou = np.sum(idx_terra & idx_aqua) / np.sum(idx_aqua | idx_terra) 
            aqua_area = np.sum(idx_aqua)
            terra_area = np.sum(idx_terra)
            joint_area = np.sum(idx_aqua & idx_terra)
            matches.append([label, c, aqua_area, terra_area, joint_area, iou])
            
    matches = pd.DataFrame(matches, columns=['aqua_label', 'terra_label', 'aqua_area',
                                             'terra_area', 'joint_area', 'iou'])
    
    if len(matches) > 0:
        for row, data in matches.iterrows():
            r_aqua, c_aqua = regions_aqua.loc[data.aqua_label.astype(int), ['centroid-0', 'centroid-1']]
            r_terra, c_terra = regions_terra.loc[data.terra_label.astype(int), ['centroid-0', 'centroid-1']]
            matches.loc[row, 'r_aqua'] = r_aqua
            matches.loc[row, 'c_aqua'] = c_aqua
            matches.loc[row, 'r_terra'] = r_terra
            matches.loc[row, 'c_terra'] = c_terra
        
        matches['drows'] = matches['r_aqua'] - matches['r_terra']
        matches['dcols'] = matches['c_aqua'] - matches['c_terra']

    fwd_filter = matches.sort_values(['aqua_label', 'iou'])
    fwd_filter = fwd_filter.loc[~fwd_filter['aqua_label'].duplicated(keep='last')]
    bwd_filter = fwd_filter.sort_values(['terra_label', 'iou'])
    bwd_filter = bwd_filter.loc[~bwd_filter['terra_label'].duplicated(keep='last')]
    matches = bwd_filter.copy()

    return regions_aqua, regions_terra, matches, lb_aqua, lb_terra, tc_aqua, tc_terra

######## Create dictionary with the regionprops tables #######
matches_dict = {}
aqua_dict = {}
terra_dict = {}

for case_idx in range(len(cases)):
    regions_aqua, regions_terra, matches, lb_aqua, lb_terra, tc_aqua, tc_terra = load_case(cases[case_idx])
    case_number = cases[case_idx].split('-')[0]
    matches['case_number'] = case_number
    regions_aqua['case_number'] = case_number
    regions_terra['case_number'] = case_number
    
    matches_dict[case_number] = matches
    aqua_dict[case_number] = regions_aqua
    terra_dict[case_number] = regions_terra


    
    regions_aqua.to_csv('../data/floe_property_tables/aqua/' + cases[case_idx].replace('labeled_floes.tiff', 'floe_properties.csv'))
    regions_terra.to_csv('../data/floe_property_tables/terra/' + cases[case_idx].replace(
        'labeled_floes.tiff', 'floe_properties.csv').replace('aqua', 'terra'))
                         
######### Lookup dictionary for manual selection of overlapping floes #########
# Format: case number, (all, some, none), empty if all/none, list of aqua labels if label = some
# all: keep every overlap
# all_iou: keep all over 0.5, drop the tentative
# some: keep ones in list
# none: no good matches

 
filter_dict = {
    '001': ['all', []], 
    '004': ['all', []],
    '005': ['none', []], # Complex case, some potential for matches but its ambiguous
    '006': ['all', []], # Some cases esp. in Aqua image that should be refined
    '007': ['some', [18, 20, 41, 48]], # Additional matches possible with velocity estimates
    '008': ['all', []], # Refine Terra image -- additional large floe possible to discern
    '009': ['some', [1, 5, 8]], # Additional floes possible in Aqua
    '010': ['all', []],
    '011': ['all', []],
    '012': ['all', []], # Fix the Terra image -- we are not including floes that intersect the edge
    '013': ['all', []], # Some in top right overlap due to drift, not matching
    '015': ['all', []], # In bottom right, can link after accounting for drift
    '016': ['some', [3, 6, 9, 12, 13, 18, 27, 32, 42, 94, 118]], # After correcting floe borders, additional matches are possible
    '019': ['some', [5,  7, 13, 14, 18, 19, 24, 25, 31, 32, 33, 48, 65, 67, 70]], # One on the right hand side should be skipped
    '021': ['all', []],
    '022': ['some', [3, 5, 17, 29, 36, 39, 40, 44, 45, 46, 47, 66]],
    '023': ['all', []],
    '029': ['all', []],   # Some to add with larger drift
    '033': ['all', []],
    '043': ['all', []],
    '044': ['all', []],
    '046': ['all', []],
    '047': ['all', []],
    '048': ['all', []],
    '051': ['none', []], # missing floe in Terra image
    '053': ['all', []],
    '056': ['all', []], 
    '062': ['all', []],
    '065': ['none', []], # looks like large drift between images
    '067': ['none', []], # large drift -- potential for linking manually
    '068': ['none', []], # large drift -- potential for manual links
    '071': ['none', []], # large drift -- potential for manual links
    '075': ['all', []],
    '081': ['some', [8]],
    '093': ['all', []], # floe in Aqua image looks like it's incorrectly split into two
    '095': ['all', []],
    '100': ['all', []],
    '107': ['all', []],
    '108': ['all', []],
    '109': ['some', [10, 17]], # Worth checking the three small floes - seems the automatic gets mixed up.
    '110': ['all', []],
    '111': ['all', []], # Additional matches possible with drift correction
    '112': ['some', [4,  5,  7, 10, 12, 14, 16, 17, 20, 22, 26, 27,
                     28, 29, 30, 33, 36, 38, 41, 43, 51, 52, 54, 55, 57, 58, 62]], # dropped 8 with low IOU 
    '115': ['some', [2, 4, 9, 10, 18, 21, 23]],
    '116': ['some', [13, 15, 26]], # potential for additional floes to be identified
    '118': ['all', []],
    '119': ['all', []],
    '121': ['some', [1, 2, 12, 47, 49, 50, 53, 56, 60, 61, 62, 63,
                     64, 65, 67, 69, 75]], # Potential for more matches
    '128': ['all', []],
    '129': ['all', []],
    '138': ['all', []], # Might need to look closer -- seems pretty good though.
    '141': ['all', []],
    '144': ['all', []],
    '148': ['some', [8, 9, 10, 11, 17]], # Some can be fixed with drift correction
    '150': ['none', []], # Looks like the terra images are misaligned, double check
    '152': ['all', []],
    '155': ['some', [29]], # Multiple matcehs possible with drift correction
    '156': ['all_iou', []], # unclear other matches
    '157': ['all', []],
    '158': ['all', []],
    '160': ['all', []],
    '161': ['all', []],
    '164': ['some', [12, 13, 14]], # Drift correction can add matches
    '166': ['none', []], # This one has tons, but drift correction needed
    '171': ['none', []], # Overlap leads to wrong match
    }

#### Create lookup dict for pairs to add #####
add_dict = {}
# tbd

def adjust_matches(matches_dict, filter_dict, add_dict = {}):
    """
    matches_dict: {case_number, match_df}
    filter_dict: of the flagged matches, which should be kept?
    add_dict: of the unflagged floes, which pairs should be added?
    Returns updated matches dataframe.
    TBD: add column for match type. high_iou, low_iou, added
    TBD: removed -- any auto matches that shouldn't be there?
    """
    updated_matches = {}
    for case in matches_dict:
        init_idx = list(matches_dict[case].loc[matches_dict[case]['iou'] >= 0.5].index)
        if case in filter_dict:
            if filter_dict[case][0] == 'all':
                updated_matches[case] = matches_dict[case].copy()
                updated_matches[case]['method'] = 'high_iou'
            elif filter_dict[case][0] == 'some':
                manu_idx = [x for x in matches_dict[case].index if matches_dict[case].loc[x, 'aqua_label'] in filter_dict[case][1]]
                keep_idx = init_idx + manu_idx
                keep_idx.sort()
                updated_matches[case] = matches_dict[case].loc[keep_idx, :].copy()
                updated_matches[case].loc[init_idx, 'method'] = 'high_iou'
                updated_matches[case].loc[manu_idx, 'method'] = 'low_iou_manual'
            else:
                print(case)
                updated_matches[case] = pd.DataFrame(data=np.nan, columns=matches_dict[case].columns, index=[0])
        
    return updated_matches

updated_matches = adjust_matches(matches_dict, filter_dict)
for case_number in updated_matches:
    for case in cases:
        if case_number == case.split('-')[0]:
            fname = case.replace('labeled_floes.tiff', 'floe_properties.csv').replace('aqua', 'matched')
            updated_matches[case_number].to_csv('../data/floe_property_tables/matched/' + fname)
        # else:
        #     fname = case_number + '-matched-floe_properties.csv'
            
    
    


                                            
##### Overlay floes and matches ######
# TBD: add step to load landfast and landmask images
# TBD: add check that matches df exists, only plot if it does
def plot_match_images(regions_aqua, regions_terra, updated_matches, lb_aqua, lb_terra, tc_aqua, tc_terra, case_number):
    """
    Overlay images with candidate matches
    """
    lbmask = lambda im: np.ma.masked_array(im, mask=im==0)
    d = sk.morphology.disk(3)
    
    matches = updated_matches[case_number]

    fig, axs = pplt.subplots(ncols=3, nrows=2)
    ### Format ###
    # | TC_aqua         | TC_terra         | Floes Alone | 
    # | TC_aqua + masks | TC_terra + masks | Floes Matched |
    
    axs[0, 0].imshow(tc_aqua)
    axs[0, 1].imshow(tc_terra)
    axs[1, 0].imshow(tc_aqua)
    axs[1, 1].imshow(tc_terra)
    
    # Overlay outlines
    axs[1, 0].imshow(lbmask(lb_aqua - sk.morphology.erosion(lb_aqua, footprint=d)), color='r')
    axs[1, 1].imshow(lbmask(lb_terra - sk.morphology.erosion(lb_terra, footprint=d)), color='b')

    # # Overlay masks
    # axs[1, 0].imshow(lbmask(lf_aqua), color='gold')
    # axs[1, 1].imshow(lbmask(lf_terra), color='gold')
    # axs[1, 0].imshow(lbmask(lm_aqua), color='tab:red')
    # axs[1, 1].imshow(lbmask(lm_terra), color='tab:red')

    #### Floe shapes ####
    axs[0, 2].imshow(lbmask(lb_aqua), vmin=0, vmax=1, color='r', alpha=0.25)
    axs[0, 2].imshow(lbmask(lb_terra), vmin=0, vmax=1, color='b', alpha=0.25)
    axs[1, 2].imshow(lbmask(lb_aqua), vmin=0, vmax=1, color='r', alpha=0.25)
    axs[1, 2].imshow(lbmask(lb_terra), vmin=0, vmax=1, color='b', alpha=0.25)
    for row, data in matches.iterrows():
        # Overlay outlines of matched objects
        lb_aqua_floe = (lb_aqua == data.aqua_label).astype(int)
        lb_aqua_floe = lb_aqua_floe - sk.morphology.erosion(lb_aqua_floe, footprint=d)
        lb_terra_floe = (lb_terra == data.terra_label).astype(int)
        lb_terra_floe = lb_terra_floe - sk.morphology.erosion(lb_terra_floe, footprint=d)
        axs[1, 2].imshow(lbmask(lb_aqua_floe), color='r')
        axs[1, 2].imshow(lbmask(lb_terra_floe), color='b')
        if data.iou < 0.5:
            axs[1, 2].plot(data.c_aqua, data.r_aqua, c='k', m='+')
            axs[1, 2].text(data.c_aqua, data.r_aqua,
                           int(data.aqua_label), c='r')
    # Things for legend:
    # landfast ice, landmask, aqua floes, terra floes, matched.
    for c, lw, m, ms, lb in zip(['r', 'b', 'k', 'tab:red', 'gold'], # c
                                [1, 1, 0, 0, 0], # lw
                                ['', '', '+', 's', 's'], # m
                                [0, 0, 5, 5, 5], # ms
                                ['Aqua Floes', 'Terra Floes', 'IOU < 0.5', 'Land', 'Landfast Ice']):
        axs[1, 2].plot([], [], color=c, lw=lw, ms=ms, m=m, label=lb)
    fig.legend(loc='b', ncols=5)
    axs[1, 2].format(yreverse=True, xreverse=False)
    case = ''
    for idx in range(len(cases)):
        if case_number == cases[idx].split('-')[0]:
            case = cases[idx]
    for ax, title in zip(axs, ['Aqua truecolor', 'Terra truecolor', 'Labeled Floes',
                               'Aqua labeled', 'Terra labeled',  'Paired Floes']):
        ax.format(title=title)
        cn, rg, dt, sat, imtype = case.split('-')
    axs.format(yreverse=True, xreverse=False,
               suptitle='Case {cn}: {rg}, {dt}'.format(cn=cn,
               rg=rg.replace('_', ' ').title(),
               dt=pd.to_datetime(dt).strftime('%Y-%m-%d')))
    
    fig.save('../data/validation_images/matching_test_images/' + case_number + '_test_image.png', dpi=300)
    pplt.close(fig)

for case_idx in range(len(cases)):
    regions_aqua, regions_terra, matches, lb_aqua, lb_terra, tc_aqua, tc_terra = load_case(cases[case_idx])
    case_number = cases[case_idx].split('-')[0]
    # matches_dict[case_number] = matches
    plot_match_images(regions_aqua, regions_terra, updated_matches, lb_aqua, lb_terra, tc_aqua, tc_terra, case_number)