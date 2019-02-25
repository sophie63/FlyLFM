function tif_scp(instance_path, directory_path)
    % Upload .tif files to AWS via scp
    fileID = fopen(instance_path,'r');
    formatSpec='%s';
    D = textscan(fileID,formatSpec,'Delimiter',',');
    fclose(fileID);
    S=size(D{1})/2

    folder_name =directory_path;
    a=dir([folder_name '/*.tif']);
    out=size(a,1);
    folder_name2=strrep(folder_name,' ','\\ ')

    % Copy data to server
    parfor i=1:S(1)
        % Create ssh connection to EC2 instance
        % login information for EC2 instance
        username = 'ubuntu';
        hostname = D{1}{2*(i-1)+1}
        
        private_key_path = strcat('/home/sophie/Downloads/', D{1}{2*(i-1)+2},'.pem');
        ssh2_conn = ssh2_config_publickey(hostname, username, private_key_path, '');
        for j=(1+ceil((i-1)*out/S(1)):ceil(i*out/S(1)))
            % upload file using scp
            file_to_upload = strcat(folder_name2(end-9:end), '_', num2str(j, '%05d'), '.tif')
            %file_to_upload = strcat(folder_name2(end-9:end), '_', num2str(j), '.tif')            
            path_to_upload = '/home/ubuntu';
            ssh2_conn = scp_put(ssh2_conn, file_to_upload, path_to_upload, folder_name);
        end
        % close connection to instance
        ssh2_conn = ssh2_close(ssh2_conn)
    end
end
