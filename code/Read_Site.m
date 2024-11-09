function [site_all,site_num]=Read_Site(file_site)
%% Read Station List
% args:
%     file_site: sitelist file
% made by Caiya Yue @ CUMTB and CASM
% ----
%%
site_num=0;
site_all.name = zeros(300,1);
fid=fopen(file_site,'r');
while(1)
    strTemp=fgets(fid);    
    if(strTemp == -1)
        break;
    end
    if(strcmp(strTemp(1:2),'X_')); continue; end
    site_num=site_num+1;
    site_all(site_num).name=strTemp;
end
fclose(fid);