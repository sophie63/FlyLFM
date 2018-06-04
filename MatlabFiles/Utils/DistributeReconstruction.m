fileID = fopen('~/InstancesListc','r');
formatSpec='%s';
D = textscan(fileID,formatSpec,'Delimiter',',')
fclose(fileID);
S=size(D{1})/2;

folder_name = uigetdir;
a=dir([folder_name '/*.tif']);
out=size(a,1);
folder_name2=strrep(folder_name,' ','\\ ');

% Copy data to server
for i=1:S(1)
    cmdStr=['for i in {',num2str((i-1)*ceil(out/S(1))),'..',num2str(i*ceil(out/S(1))),'}; do printf -v j %%04d $i; scp -o "StrictHostKeyChecking no" -r -i ~/Downloads/',D{1}{2*(i-1)+2},'.pem ',folder_name2,'/',folder_name2(end-9:end),'_$j.tif ubuntu@',D{1}{2*(i-1)+1},':~/.; done;\n'];
    fprintf(cmdStr);
end
%%% ##########Copy 1######

for i=1:S(1)
    cmd=['scp -o "StrictHostKeyChecking no" -r -i ~/Downloads/',D{1}{2*(i-1)+2},'.pem 100624_ss1_f1250_6microns.lfc ubuntu@',D{1}{2*(i-1)+1},':~/.\n'];
    fprintf(cmd);
end
%%% ########## copy2 #####
% % Connect to instances 
% for i=1:S(1)
%     cmd=['sudo ssh -o "StrictHostKeyChecking no" -i ~/Downloads/',D{1}{2*(i-1)+2},'.pem ubuntu@',D{1}{2*(i-1)+1},'\n'];
%     fprintf(cmd);
% end

% Code per instance, each goes to different terminal
for i=1:S(1)
%    cmd0=['gnome-terminal "'];
%    fprintf(cmd0);
    cmd=['ssh -o "StrictHostKeyChecking no" -i ~/Downloads/',D{1}{2*(i-1)+2},'.pem ubuntu@',D{1}{2*(i-1)+1},' \n'];
    cmdStr=['for i in {',num2str((i-1)*ceil(out/S(1))),'..',num2str(i*ceil(out/S(1))),'}; do printf -v j \"%%04d\" $i; python2.7 ~/stanford_lfanalyze_v0.4/lfdeconvolve.py ~/',folder_name2(end-9:end),'_$j.tif -c 100624_ss1_f1250_6microns.lfc -o ~/',folder_name2(end-5:end),'ss1/',folder_name2(end-5:end),'ss1-$j.tif --max-iter=30; rm ~/',folder_name2(end-9:end),'_$j.tif; done; \n'];
    
    fprintf(cmd);
    fprintf(cmdStr);

end
% ############ Copy3 ######

% Can be separated to terminals to save time copying
for i=1:S(1)
    cmd=['scp -o "StrictHostKeyChecking no" -r -i ~/Downloads/',D{1}{2*(i-1)+2},'.pem ubuntu@',D{1}{2*(i-1)+1},':~/*ss1 /media/sophie/New\\ Volume/.\n'];    
    fprintf(cmd);
end
  
% 
% for i=1:S
%     
%     cmd=['sudo ssh -o "StrictHostKeyChecking no" -i ~/Downloads/',D{1}{2*(i-1)+2},'.pem ubuntu@',D{1}{2*(i-1)+1},' ''rm *tif''&\n'];
%     fprintf(cmd);
% end
% 
% ubuntu@ec2-54-215-205-131.us-west-1.compute.amazonaws.com:~/*ss1 /media/sophie/New\ Volume/.
% 
% sudo ssh -o "StrictHostKeyChecking no" -i ~/Downloads/Jan2018.pem ubuntu@ec2-13-57-42-137.us-west-1.compute.amazonaws.com 'rm *tif'
% 
% 
% 
% gnome-terminal --window-with-profile=NOCLOSEPROFILE -e "ssh -i ~/Downloads/Jan2018.pem ubuntu@ec2-13-57-189-223.us-west-1.compute.amazonaws.com 'for i in {466..699}; do printf -v j "%04d" $i; python2.7 ~/stanford_lfanalyze_v0.4/lfdeconvolve.py ~/Data100660_$j.tif -c 100659_ss1_f3125_6microns.lfc -o ~/100660ss1/100660ss1-$j.tif --max-iter=30; done;'"
% 
% gnome-terminal "ssh -i ~/Downloads/Jan2018.pem ubuntu@ec2-13-57-189-223.us-west-1.compute.amazonaws.com 'for i in {466..699}; do printf -v j "%04d" $i; python2.7 ~/stanford_lfanalyze_v0.4/lfdeconvolve.py ~/Data100660_$j.tif -c 100659_ss1_f3125_6microns.lfc -o ~/100660ss1/100660ss1-$j.tif --max-iter=30; done;'"
% 
% gnome-terminal --window-with-profile=NOCLOSEPROFILE -e "ssh -i ~/Downloads/Jan2018.pem ubuntu@ec2-13-57-189-223.us-west-1.compute.amazonaws.com"
% 
% 
% 
% 
%    gnome-terminal -x sh -c "!!; bash" 
% !!; exec bash\""