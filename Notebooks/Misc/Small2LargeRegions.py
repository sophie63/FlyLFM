import numpy as np
from scipy import stats
from scipy import io
import nibabel as nb

filename='/home/sophie/RegionList'
with open(filename) as f:
    content = f.readlines()
Names=[Line.split('\t') for Line in content]
RegionName=[Names[i][0] for i in range(75)]
Num=[int(Names[i][2]) for i in range(75)]

LargerRegionDic={'':'','AME_R':'OL','LO_R':'OL','NO':'CX','BU_R':'CX','PB':'CX','LH_R':'LH','LAL_R':'LX','SAD':'PENP',
                'CAN_R':'PENP','AMMC_R':'PENP','ICL_R':'INP','VES_R':'VMNP','IB_R':'INP','ATL_R':'INP','CRE_R':'INP',
                'MB_PED_R':'MB','MB_VL_R':'MB','MB_ML_R':'MB','FLA_R':'PENP','LOP_R':'OL','EB':'CX','AL_R':'AL',
                'ME_R':'OL','FB':'CX','SLP_R':'SNP','SIP_R':'SNP','SMP_R':'SNP','AVLP_R':'VLNP','PVLP_R':'VLNP',
                'IVLP_R':'VLNP','PLP_R':'VLNP','AOTU_R':'VLNP','GOR_R':'VMNP','MB_CA_R':'MB','SPS_R':'VMNP',
                'IPS_R':'VMNP','SCL_R':'INP','EPA_R':'VMNP','GNG':'GNG','PRW':'PENP','GA_R':'LX','AME_L':'OL',
                'LO_L':'OL','BU_L':'CX','LH_L':'LH','LAL_L':'LX','CAN_L':'PENP','AMMC_L':'PENP','ICL_L':'INP',
                'VES_L':'VMNP','IB_L':'INP','ATL_L':'INP','CRE_L':'INP','MB_PED_L':'MB','MB_VL_L':'MB',
                'MB_ML_L':'MB','FLA_L':'PENP','LOP_L':'OL','AL_L':'AL','ME_L':'OL','SLP_L':'SNP','SIP_L':'SNP',
                'SMP_L':'SNP','AVLP_L':'VLNP','PVLP_L':'VLNP','IVLP_L':'VLNP','PLP_L':'VLNP','AOTU_L':'VLNP',
                'GOR_L':'VMNP','MB_CA_L':'MB','SPS_L':'VMNP','IPS_L':'VMNP','SCL_L':'INP','EPA_L':'VMNP','GA_L':'LX'}

LargerRegionInd={'OL':1,'VLNP':2,'VMNP':3,'AL':4,'MB':5,'LH':6,'SNP':7,'CX':8,'LX':9,'INP':10,'PENP':11,'GNG':12,'':13}
SmallRegionLargeId={k: LargerRegionInd[v] for (k,v) in LargerRegionDic.items()}

with open('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/paper_dataset_name_lists.txt') as f:
	filenameM_list = [line.strip() for line in f.readlines()]

for filenameM in filenameM_list:
	print filenameM
	img1 = nb.load(filenameM)
	Masks = img1.get_data()
	Sm=Masks.shape
	Masks=np.array(Masks)

	MaskLargeReg=np.zeros([Sm[0],Sm[1],Sm[2],12])
	for k in range(13):
	    for j in range(74):
	        if SmallRegionLargeId[RegionName[j+1]]==k:
	            MaskLargeReg[:,:,:,k-1]=MaskLargeReg[:,:,:,k-1]+Masks[:,:,:,Num[j+1]-1]
	        
	NewName='/'.join(filenameM.split('/')[:-1])+'/JFRCTransformedLarge.nii'
	nim=nb.Nifti1Image(MaskLargeReg,np.eye(4))
	nb.save(nim,NewName)






