#!/bin/bash

# Calibration and Reconstruction Automation
# Dominick Nguyen
#
# Note: Requires prexisting g2.8x large instance
#       Also, edit all paths to pem files for ssh, scp commands
#
# args:
# 1.tif_file - AVG Tif file path
# 2.instance_path - Path to any text file to write instances + pems
# 3.directory_path - Path to directory containing image sequences
# 4.num_slices - Number of slices
# 5.region - Region for reconstruction step. See note in reconstruction script for more
# 6.instance_id - ec2 instance id for reconstruction
# 7.reconstr_path - path to store files after reconstruction step
# 8.save_path - path to save for tiff to nii
# usage ./calibrate_script.sh <tif_file> <instance_file> <directory_file> <num_slices> <region>
#         <instance_id> <reconstr_path> <save_path>

# Assumes name of file is formatted Data100XXX
FILE_NAME=$(basename $1)
TIF_NUM=$(echo $FILE_NAME | cut -c9-14)
TIF_FILE_NAME="${TIF_NUM}_ss1_f1250_6microns.lfc"

echo $FILE_NAME
echo $TIF_NUM
echo $TIF_FILE_NAME

# start up instance for calibration(Note: change to correct region depending on passed instance)
aws --region us-west-1 ec2 start-instances --instance-ids $6

sleep 1m

# get public dns of instance 
public_dns="$(aws --region us-west-1 ec2 describe-instances --instance-ids $6 --query "Reservations[].Instances[].PublicDnsName" --output=text)"

echo $public_dns

# Note:replace with own pem location for 
# scp tif file to ec2 instance
scp -o StrictHostKeyChecking=no -i ~/Downloads/Feb2018.pem $1 "ubuntu@${public_dns}:~/" 

# run command for calibration
ssh -o StrictHostKeyChecking=no -i ~/Downloads/Feb2018.pem "ubuntu@${public_dns}" "python ~/stanford_lfanalyze_v0.4/lfcalibrate.py ${FILE_NAME} --pixel-size 6.5 --pitch 125 --focal-length 3125 --magnification 40 --na 0.8 --tubelens-focal-length 180.0 --wavelength 510 --medium-index 1.33 --num-slices ${4} --um-per-slice 6.0 --num-threads 15 --supersample 1 -o ${TIF_FILE_NAME}"

# if the calibration was successful
# edit scp command
RESULT=$?
if [ $RESULT -eq 0 ]; then
	echo "Complete. Downloading calibration file back to local as ~/${TIF_FILE_NAME}"
	scp -o StrictHostKeyChecking=no -i ~/Downloads/Feb2018.pem "ubuntu@${public_dns}:~/${TIF_FILE_NAME}" .
	aws --region us-west-1 ec2 stop-instances --instance-ids $6
	echo "Download done. Stopping instance and Running reconstruction script."
	./run_matlab.sh $2 $3 $5 $TIF_FILE_NAME $7 $8
else
	echo "Nonzero Exit Code. Failed to generate calibration file. Stopping instance now."
	aws --region us-west-1 ec2 stop-instances --instance-ids $6
fi
