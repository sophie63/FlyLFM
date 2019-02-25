function run_deconvolve(instance_path, directory_path, lfc_file)
    % Upload calibration file / Run lfdeconvolve
    fileID = fopen(instance_path,'r');
    formatSpec='%s';
    D = textscan(fileID,formatSpec,'Delimiter',',');
    fclose(fileID);
    S=size(D{1})/2;

    folder_name = directory_path;
    a=dir([folder_name '/*.tif']);
    out=size(a,1);
    folder_name2=strrep(folder_name,' ','\\ ')
    
    fprintf('start uploading calibration file')
    parfor i=1:S(1)
        username = 'ubuntu'
        i   
        hostname = D{1}{2*(i-1)+1}
        private_key_path = strcat('/home/sophie/Downloads/', D{1}{2*(i-1)+2},'.pem')
        ssh2_conn = ssh2_config_publickey(hostname, username, private_key_path, '')
        ssh2_conn = scp_put(ssh2_conn, lfc_file, '/home/ubuntu/', '.')
        ssh2_conn = ssh2_close(ssh2_conn)
    end
    fprintf('finished uploading calibration file')
    
    parfor i=1:S(1)
        username = 'ubuntu'
        hostname = D{1}{2*(i-1)+1}
        private_key_path = strcat('/home/sophie/Downloads/', D{1}{2*(i-1)+2},'.pem')
        ssh2_conn = ssh2_config_publickey(hostname, username, private_key_path, '')
        for j=(1+ceil((i-1)*out/S(1)):ceil(i*out/S(1)))
            file_num = num2str(j, '%05d')
            %file_num=num2str(j);
            cmd_str = ['python2.7 ~/stanford_lfanalyze_v0.4/lfdeconvolve.py /home/ubuntu/' folder_name2(end-9:end) '_' file_num '.tif' ' ' '-c ' lfc_file ' -o /home/ubuntu/' folder_name2(end-5:end) 'ss1/' folder_name2(end-5:end) 'ss1-' file_num '.tif --max-iter=30'];
            cmd_str2 = strcat('rm ~/',folder_name2(end-9:end), '_', file_num, '.tif')
            ssh2_conn = ssh2_command(ssh2_conn, cmd_str);
            ssh2_conn = ssh2_command(ssh2_conn, cmd_str2)
        end
        strcat('Set ',num2str(i),' finished')
        ssh2_conn = ssh2_close(ssh2_conn)
    end
    'All sets finished'
end
