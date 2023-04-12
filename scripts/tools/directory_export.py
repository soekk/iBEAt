import os
import dbdicom as db

source_dir = 'C://Users//UOS//Desktop//Data_iBEAT//Leeds_Sample'
resultsPath = 'C://Users//UOS//Desktop//masks_iBEAt'

file_list = os.listdir(source_dir)

out_desc        = 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'
left_mask_desc  = 'LK'
right_mask_desc = 'RK'


for file_name in file_list:
    
    studie = source_dir + "//" + file_name

    folder = db.database(path=studie)

    out_ph = folder.series(SeriesDescription=out_desc)
    mask_l = folder.series(SeriesDescription=left_mask_desc)
    mask_r = folder.series(SeriesDescription=right_mask_desc)

    exportToFolder = out_ph + mask_l + mask_r

    for series in exportToFolder:
        print(series.SeriesDescription)    
        series.export_as_dicom(resultsPath)


