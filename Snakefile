from os.path import join
from glob import glob
import pandas as pd

configfile: 'config.yml'

#load participants.tsv file (list of HCP subject IDs),
df = pd.read_table(config['participants_tsv'])
subjects = df.participant_id.to_list() 

wildcard_constraints:
    subject="[0-9]+"

rule all:
    input: expand('funcparc/clustering/seed-ZIR_method-spectralcosine_k-{k}_cluslabels.nii.gz',k=range(2,config['max_k']+1),allow_missing=True)

rule unzip_packages:
    input: 
        packages = expand(join(config['hcp1200_zip_dir'],'{subject}_{package}.zip'), package=config['hcp_package_dict'].values(), allow_missing=True)
    params:
        out_dir = 'hcp1200',
        files_in_pkg = expand('{filename}',filename=config['hcp_package_dict'].keys())
    output: 
        files = expand('hcp1200/{filename}',filename=config['hcp_package_dict'].keys())
    group: preprocessing
    run: 
        for pkg,file in zip(input.packages, params.files_in_pkg):
            shell('unzip {pkg} {file} -d {params.out_dir}')

rule cifti_separate:
    input: 
        dtseries = 'hcp1200/{subject}/MNINonLinear/Results/rfMRI_REST2_7T_AP/rfMRI_REST2_7T_AP_Atlas_1.6mm.dtseries.nii'
    output: 
        lh = 'funcparc/{subject}/input/rfMRI_REST2_7T_AP.L.59k_fs_LR.func.gii',
        rh = 'funcparc/{subject}/input/rfMRI_REST2_7T_AP.R.59k_fs_LR.func.gii'
    group: preprocessing
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000
    shell:
        'wb_command -cifti-separate {input} COLUMN -metric CORTEX_LEFT {output.lh} -metric CORTEX_RIGHT {output.rh}'

rule prepare_subcort:
    input:
        vol = 'hcp1200/{subject}/MNINonLinear/Results/rfMRI_REST2_7T_AP/rfMRI_REST2_7T_AP.nii.gz',
        rois = config['subcort_atlas']
    output:
        out = 'funcparc/{subject}/input/rfMRI_REST2_7T_AP_AtlasSubcortical.nii.gz'
    params:
        sigma = 1.6,
        temp = 'funcparc/{subject}/temp'
    group: preprocessing
    singularity: config['singularity_connectomewb']
    log: 'logs/prepare_subcort/{subject}.log'
    threads: 8
    resources:
        mem_mb = 32000
    shell: 'scripts/prep_subcortical.sh {input.vol} {input.rois} {params.temp} {params.sigma} {output.out} &> {log}'

rule create_dtseries:
    input: 
        vol = rules.prepare_subcort.output.out,
        rois = config['subcort_atlas'],
        lh = rules.cifti_separate.output.lh,
        rh = rules.cifti_separate.output.rh
    output: 'funcparc/{subject}/input/rfMRI_REST2_7T_AP.59k_fs_LR.dtseries.nii'
    group: preprocessing
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000
    shell:
        'wb_command -cifti-create-dense-timeseries {output} -volume {input.vol} {input.rois} -left-metric {input.lh} -right-metric {input.rh}'

rule extract_confounds:
    input:
        vol = 'hcp1200/{subject}/MNINonLinear/Results/rfMRI_REST2_7T_AP/rfMRI_REST2_7T_AP.nii.gz',
        rois = 'hcp1200/{subject}/MNINonLinear/ROIs/Atlas_wmparc.1.60.nii.gz',
        movreg = 'hcp1200/{subject}/MNINonLinear/Results/rfMRI_REST2_7T_AP/Movement_Regressors_dt.txt'
    output: 'funcparc/{subject}/input/confounds.tsv'
    group: preprocessing
    log: 'logs/extract_confounds/{subject}.log'
    resources:
        mem_mb = 32000
    script: 'scripts/extract_confounds.py'

rule clean_tseries:
    input:
        dtseries = rules.create_dtseries.output,
        confounds = rules.extract_confounds.output
    output: 'funcparc/{subject}/input_cleaned/rfMRI_REST2_7T_AP.59k_fs_LR.dtseries.nii'
    group: preprocessing
    singularity: config['singularity_ciftify']
    log: 'logs/clean_dtseries/{subject}.log'
    threads: 8
    resources:
        mem_mb = 32000
    shell:
        'ciftify_clean_img --output-file={output} --detrend --standardize --confounds-tsv={input.confounds} --low-pass=0.08 --high-pass=0.009 --tr=1 --verbose {input.dtseries} &> {log}'

rule parcellate_tseries:
    input:
        dtseries = rules.clean_tseries.output,
        rois = config['hcpmmp_atlas']
    output: 'funcparc/{subject}/input_parcellated/rfMRI_REST2_7T_AP.59k_fs_LR.ptseries.nii'
    group: preprocessing
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000
    shell:
        'wb_command -cifti-parcellate {input.dtseries} {input.rois} COLUMN {output}'

rule compute_correlation:
    input:
        ptseries = rules.parcellate_tseries.output,
        vol = rules.clean_tseries.output
    output: 'funcparc/{subject}/output/correlation_matrix.npz'
    group: preprocessing
    script: 'scripts/compute_correlation.py'

rule combine_correlation:
    input: expand('funcparc/{subject}/output/correlation_matrix.npz',subject=subjects,allow_missing=True)
    output: 'funcparc/clustering/correlation_matrix_group.npz'
    group: clustering
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
        rois = config['subcort_atlas'],
    output:
        niftis = expand('funcparc/clustering/seed-ZIR_method-spectralcosine_k-{k}_cluslabels.nii.gz',k=range(2,config['max_k']+1),allow_missing=True),
        labels = 'funcparc/clustering/cluster_labels.csv'
    params:
        max_k = config['max_k']
    group: clustering        
    script: 'scripts/spectral_clustering.py'
