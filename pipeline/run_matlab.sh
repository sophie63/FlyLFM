#!/bin/bash

echo "Creating new reconstruction instances"

if [[ "$3" = "us-west-1" ]]; then
    pem="Feb2018"
    instance_num=5
elif [[ "$3" = "us-west-2" ]]; then
    pem="February2018V2"
    instance_num=5
elif [[ "$3" = "us-east-1" ]]; then
    pem="February2018V3"
    instance_num=10
fi

aws --region $3 ec2 run-instances --launch-template LaunchTemplateName=reconstruction_template --count $instance_num

echo "writing public dns names to file ~/InstancesListc"

aws --region $3 ec2 describe-instances --query "Reservations[*].Instances[*].[PublicDnsName]" --output=text | xargs -i echo "{},$pem" > $1

sleep 1m

echo "Beginning upload"

matlab -nodesktop -nodisplay -r "tif_scp('$1','$2');quit" &
(sleep 30 && echo "Beginning deconvolve" && matlab -nodesktop -nodisplay -r "run_deconvolve('$1','$2','$4');quit") &

wait
echo "Uploading and deconvolving complete"
echo "downloading completed files"

matlab -nodesktop -nodisplay -r "download_ss('$1', '$2', '$5');quit"

wait
echo "done"
echo "terminating all instances"

# terminate all instances
aws --region $3 ec2 terminate-instances --instance-ids $(aws --region $3 ec2 describe-instances --filters "Name=instance-state-name,Values=pending,running,stopped,stopping" --query "Reservations[].Instances[].[InstanceId]" --output text | tr '\n' ' ')

# Autotiff2matlab.py
python Auto_tiff2matlab.py $5 $6
