from os.path import join
from glob import glob
import pandas as pd

configfile: 'config_3T.yml'

#load participants.tsv file (list of HCP subject IDs),
df = pd.read_table(config['participants_tsv'])
subjects = df.participant_id.to_list() 

run = config['run']
vox_res = config['voxel_size']
vtx_res = config['vertex_nr'] 

wildcard_constraints:
    subject="[0-9]+"

rule all:
    input: expand('funcparc/clustering/seed-{s}_method-spectralcosine_k-{k}_cluslabels.nii.gz',s=config['seed'], k=range(2,config['max_k']+1),allow_missing=True)

# rule unzip_packages:
#     input: 
#         packages = expand(join(config['hcp1200_zip_dir'],'{subject}_{package}.zip'), package=config['hcp_package_dict'].values(), allow_missing=True)
#     params:
#         out_dir = 'hcp1200',
#         files_in_pkg = expand('{filename}',filename=config['hcp_package_dict'].keys())
#     output: 
#         files = expand('hcp1200/{filename}',filename=config['hcp_package_dict'].keys())
#     group: preprocessing
#     run: 
#         for pkg,file in zip(input.packages, params.files_in_pkg):
#             shell('unzip {pkg} {file} -d {params.out_dir}')

# rule here to make custom subcort atlas
rule merge_roi:
    input:
        atlas = config['subcort_atlas'],
        roi = config['seed_roi']
    output: expand('funcparc/atlas/BigBrain{vs}.nii.gz',vs=vox_res)
    group: 'atlas'
    run:    
        import numpy as np
        import nibabel as nib

        atlas_in = nib.load(input.atlas)
        roi_in = nib.load(input.roi)
        roi = roi_in.get_fdata()

        atlas_new = atlas_in.get_fdata()
        atlas_new[roi==1] = 16

        img = nib.Nifti1Image(atlas_new, affine=atlas_in.affine, header=atlas_in.header)
        nib.save(img,output[0])

rule create_atlas:
    input:
        atlas = rules.merge_roi.output,
        labels = config['subcort_labels'],
        lh = config['lh_mmp'],
        rh = config['rh_mmp']
    output:
        nifti = expand('funcparc/atlas/HCP-MMP_BigBrain{vs}.nii.gz',vs=vox_res),
        cifti = expand('funcparc/atlas/HCP-MMP_BigBrain{vs}.{nv}_fs_LR.dlabel.nii',vs=vox_res,nv=vtx_res) 
    group: 'atlas'
    singularity: config['singularity_connectomewb']
    shell:
        'wb_command -volume-label-import {input.atlas} {input.labels} {output.nifti} &&'
        'wb_command -cifti-create-label {output.cifti} -volume {output.nifti} {output.nifti} -left-label {input.lh} -right-label {input.rh}'

# Prepare subcortical rs-fMRI data as done by HCP
rule prepare_subcort:
    input:
        vol = expand(join(config['hcp_dir'],'{{subject}}/MNINonLinear/Results/{r}/{r}_hp2000_clean.nii.gz'),r=run),
        rois = rules.create_atlas.output.nifti
    output:
        out = expand('funcparc/{{subject}}/input/{r}_AtlasSubcortical.nii.gz',r=run)
    params:
        sigma = 1.6,
        temp = 'funcparc/{subject}/tmp'
    group: 'preprocessing'
    singularity: config['singularity_connectomewb']
    log: 'logs/prepare_subcort/{subject}.log'
    threads: 8
    resources:
        mem_mb = 32000
    shell: 'scripts/prep_subcortical.sh {input.vol} {input.rois} {params.temp} {params.sigma} {output.out} &> {log}'

# Extract cortical timeseries data
rule cifti_separate:
    input: 
        dtseries = expand(join(config['hcp_dir'],'{{subject}}/MNINonLinear/Results/{r}/{r}_Atlas_MSMAll_hp2000_clean.dtseries.nii'),r=run)
    output: 
        lh = expand('funcparc/{{subject}}/input/{r}.L.{nv}_fs_LR.func.gii',nv=vtx_res,r=run),
        rh = expand('funcparc/{{subject}}/input/{r}.R.{nv}_fs_LR.func.gii',nv=vtx_res,r=run)
    params:
        nr_vertices = config['vertex_nr']
    group: 'preprocessing'
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000
    shell:
        'wb_command -cifti-separate {input} COLUMN -metric CORTEX_LEFT {output.lh} -metric CORTEX_RIGHT {output.rh}'

# Create CIFTI timeseries data
rule create_dtseries:
    input: 
        vol = rules.prepare_subcort.output.out,
        rois = rules.create_atlas.output.nifti,
        lh = rules.cifti_separate.output.lh,
        rh = rules.cifti_separate.output.rh
    output: expand('funcparc/{{subject}}/input/{r}.{nv}_fs_LR.dtseries.nii',nv=vtx_res,r=run) 
    group: 'preprocessing'
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000
    shell:
        'wb_command -cifti-create-dense-timeseries {output} -volume {input.vol} {input.rois} -left-metric {input.lh} -right-metric {input.rh}'

# # Extract confounds for cleaning rs-fMRI data
# rule extract_confounds:
#     input:
#         vol = rules.prepare_subcort.input.vol,
#         rois = join(config['hcp_dir'],'{subject}/MNINonLinear/ROIs/Atlas_wmparc.1.60.nii.gz'),
#         movreg = join(config['hcp_dir'],'{subject}/MNINonLinear/Results/rfMRI_REST2_7T_AP/Movement_Regressors_dt.txt')
#     output: 'funcparc/{subject}/input/confounds.tsv'
#     group: 'preprocessing'
#     log: 'logs/extract_confounds/{subject}.log'
#     resources:
#         mem_mb = 32000
#     script: 'scripts/extract_confounds.py'

# # Clean rs-fMRI data
# rule clean_tseries:
#     input:
#         dtseries = rules.create_dtseries.output,
#         #confounds = rules.extract_confounds.output
#     output: expand('funcparc/{{subject}}/input_cleaned/{r}.{nv}_fs_LR.dtseries.nii',nv=vtx_res,r=run) 
#     group: 'preprocessing'
#     singularity: config['singularity_ciftify']
#     log: 'logs/clean_dtseries/{subject}.log'
#     threads: 8
#     resources:
#         mem_mb = 32000
#     shell:
#         'ciftify_clean_img --output-file={output} --detrend --standardize --tr=0.720 --verbose {input.dtseries} &> {log}'

rule parcellate_tseries:
    input:
        dtseries = rules.create_dtseries.output,
        rois = rules.create_atlas.output.cifti
    output: expand('funcparc/{{subject}}/input_parcellated/{r}.{nv}_fs_LR.ptseries.nii',nv=vtx_res,r=run)
    group: 'preprocessing'
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000
    shell:
        'wb_command -cifti-parcellate {input.dtseries} {input.rois} COLUMN {output}'

rule compute_correlation:
    input:
        ptseries = rules.parcellate_tseries.output,
        vol = rules.create_dtseries.output
    params:
        seed = config['seed']
    output: 'funcparc/{subject}/output/correlation_matrix.npz'
    group: 'preprocessing'
    script: 'scripts/compute_correlation.py'

rule combine_correlation:
    input: expand('funcparc/{subject}/output/correlation_matrix.npz',subject=subjects,allow_missing=True)
    output: 'funcparc/clustering/correlation_matrix_group.npz'
    group: 'clustering'
    run:
        import numpy as np

        data = np.load(input[0])
        nsubjects = len(input)
        combined = np.zeros([nsubjects,data['corr'].shape[0],data['corr'].shape[1]])

        for i,npz in enumerate(input):
            data = np.load(npz)
            combined[i,:,:] = data['corr']

        np.savez(output[0], corr_group=combined,indices=data['indices'])

rule spectral_clustering:
    input:
        correlation = rules.combine_correlation.output,
        rois = rules.create_atlas.output.nifti[0]
    output:
        niftis = expand('funcparc/clustering/seed-{s}_method-spectralcosine_k-{k}_cluslabels.nii.gz',s=config['seed'], k=range(2,config['max_k']+1),allow_missing=True),
        labels = 'funcparc/clustering/cluster_labels.csv'
    params:
        max_k = config['max_k']
    group: 'clustering'
    script: 'scripts/spectral_clustering.py'
