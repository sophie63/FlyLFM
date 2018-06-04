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

LargerRegionDic={'':'','AME_R':'OL_R','LO_R':'OL_R','NO':'CX','BU_R':'CX','PB':'CX','LH_R':'LH_R','LAL_R':'LX_R','SAD':'PENP',
                'CAN_R':'PENP','AMMC_R':'PENP','ICL_R':'INP_R','VES_R':'VMNP_R','IB_R':'INP_R','ATL_R':'INP_R','CRE_R':'INP_R',
                'MB_PED_R':'MB_R','MB_VL_R':'MB_R','MB_ML_R':'MB_R','FLA_R':'PENP','LOP_R':'OL_R','EB':'CX','AL_R':'AL_R',
                'ME_R':'OL_R','FB':'CX','SLP_R':'SNP_R','SIP_R':'SNP_R','SMP_R':'SNP_R','AVLP_R':'VLNP_R','PVLP_R':'VLNP_R',
                'IVLP_R':'VLNP_R','PLP_R':'VLNP_R','AOTU_R':'VLNP_R','GOR_R':'VMNP_R','MB_CA_R':'MB_R','SPS_R':'VMNP_R',
                'IPS_R':'VMNP_R','SCL_R':'INP_R','EPA_R':'VMNP_R','GNG':'GNG','PRW':'PENP','GA_R':'LX_R','AME_L':'OL_L',
                'LO_L':'OL_L','BU_L':'CX','LH_L':'LH_L','LAL_L':'LX_L','CAN_L':'PENP','AMMC_L':'PENP','ICL_L':'INP_L',
                'VES_L':'VMNP_L','IB_L':'INP_L','ATL_L':'INP_L','CRE_L':'INP_L','MB_PED_L':'MB_L','MB_VL_L':'MB_L',
                'MB_ML_L':'MB_L','FLA_L':'PENP','LOP_L':'OL_L','AL_L':'AL_L','ME_L':'OL_L','SLP_L':'SNP_L','SIP_L':'SNP_L',
                'SMP_L':'SNP_L','AVLP_L':'VLNP_L','PVLP_L':'VLNP_L','IVLP_L':'VLNP_L','PLP_L':'VLNP_L','AOTU_L':'VLNP_L',
                'GOR_L':'VMNP_L','MB_CA_L':'MB_L','SPS_L':'VMNP_L','IPS_L':'VMNP_L','SCL_L':'INP_L','EPA_L':'VMNP_L','GA_L':'LX_L'}

LargerRegionInd={'OL_R':1,'OL_L':2,'VLNP_R':3,'VLNP_L':4,'VMNP_R':5,'VMNP_L':6,'AL_R':7,'AL_L':8,'MB_R':9,'MB_L':10,'LH_R':11,'LH_L':12,'SNP_R':13,'SNP_L':14,'CX':15,'LX_R':16,'LX_L':17,'INP_R':18,'INP_L':19,'PENP':20,'GNG':21,'':22}
SmallRegionLargeId={k: LargerRegionInd[v] for (k,v) in LargerRegionDic.items()}

with open('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/paper_dataset_mask_lists.txt') as f:
	filenameM_list = [line.strip() for line in f.readlines()]


for filenameM in filenameM_list:
	print filenameM
	img1 = nb.load(filenameM)
	Masks = img1.get_data()
	Sm=Masks.shape
	Masks=np.array(Masks)

	MaskLargeReg=np.zeros([Sm[0],Sm[1],Sm[2],21])
	for k in range(22):
	    for j in range(74):
	        if SmallRegionLargeId[RegionName[j+1]]==k:
	            MaskLargeReg[:,:,:,k-1]=MaskLargeReg[:,:,:,k-1]+Masks[:,:,:,Num[j+1]-1]
	        
	NewName='/'.join(filenameM.split('/')[:-1])+'/JFRCTransformedLargeLR.nii'
	nim=nb.Nifti1Image(MaskLargeReg,np.eye(4))
	nb.save(nim,NewName)






