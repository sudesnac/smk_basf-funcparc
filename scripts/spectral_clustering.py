import numpy as np
import pandas as pd
import nibabel as nib
from sklearn.cluster import SpectralClustering

# Mute warning in case seed has same number of voxels as target ROIs
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

# Function to save niftis
def save_label_nii (labelimg,affine,header,out_nifti):
    img = nib.Nifti1Image(labelimg,affine=affine,header=header)
    nib.save(img,out_nifti)

# Load data
data = np.load(snakemake.input.correlation[0])
correlation = data['corr_group']
indices = data['indices']

afile = snakemake.input.rois
atlas = nib.load(afile)
atlas_data = atlas.get_fdata()

# Reshape and concatenate subjects
corr = np.moveaxis(correlation,0,2)
corr_concat = corr.reshape([corr.shape[0],corr.shape[1]*corr.shape[2]])
corr_concat += 1 # Spectral clustering doesn't like negative input apparantly

# Output
out_nii_list = snakemake.output.niftis
cluster_range = range(2,snakemake.params.max_k+1)
labels = np.zeros((corr_concat.shape[0],len(cluster_range)))

# Run spectral clustering and save results to nifti
for i,k in enumerate(cluster_range):
    clustering = SpectralClustering(n_clusters=k, assign_labels="discretize",random_state=0,affinity='cosine').fit(corr_concat)
    labels[:,i] = clustering.labels_
    
    labelimg = np.zeros(atlas_data.shape)
    for j in range(0,len(atlas_data[atlas_data==16])):
        labelimg[indices[j][0],indices[j][1],indices[j][2]] = labels[j,i]+1
    print(f'i={i}, k={k},saving {out_nii_list[i]}')
    save_label_nii(labelimg,atlas.affine,atlas.header,out_nii_list[i])

# Save results to CSV file
df = pd.DataFrame(labels,columns=cluster_range)
df.to_csv(snakemake.output.labels)