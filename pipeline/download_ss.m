function download_ss(instance_path, directory_path, reconstr_path)
    fileID = fopen(instance_path,'r');
    formatSpec='%s';
    D = textscan(fileID,formatSpec,'Delimiter',',');
    fclose(fileID);
    S=size(D{1})/2

    
    folder_name = directory_path;
    a=dir([folder_name '/*.tif']);
    out=size(a,1)
    folder_name2=strrep(folder_name,' ','\\ ');
    folder_name3=reconstr_path;

    %parfor i=1:S(1)
    parfor i=1:S(1)    
        username = 'ubuntu'
        hostname = D{1}{2*(i-1)+1}
        private_key_path = strcat('/home/sophie/Downloads/', D{1}{2*(i-1)+2},'.pem')
        ssh2_conn = ssh2_config_publickey(hostname, username, private_key_path, '')
        for j=(1+ceil((i-1)*out/S(1)):ceil(i*out/S(1)))
            
            file_name = strcat(folder_name2(end-5:end),'ss1','-',num2str(j, '%05d') ,'.tif ')
            %file_name = strcat(folder_name2(end-5:end),'ss1','-',num2str(j) ,'.tif ')
            ssh2_conn = scp_get(ssh2_conn, file_name, folder_name3, strcat('/home/ubuntu/', folder_name2(end-5:end),'ss1/'))        
       
        end
        strcat('Set ',num2str(i),' finished')
        ssh2_conn = ssh2_close(ssh2_conn)
    end
    'All sets finished'
end
