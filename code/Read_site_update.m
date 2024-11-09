function [site_all_ud,site_num_ud]=Read_site_update(site_all,site_num,cfg)
%% Update by reading the calculated station file
% args:
%     site_all: sitelist
%     site_num: site number
%     cfg     : configuration information
% made by Caiya Yue @ CUMTB and CASM
% ----
%% interpreting control information
dir=cfg.dir;
yyyy=num2str(cfg.yyyy);
doy=num2str(cfg.doy);
if(cfg.doy<10); doy=strcat('00',doy); end
if(cfg.doy<100 && cfg.doy>9); doy=strcat('0',doy); end
%% update sitelist
site_num_ud=0;
site_all_ud.name = zeros(300,1);
fout=fopen(strcat(cfg.out_dir,'/site_name_update.txt'),'w');
for s=1:site_num
    ln=length(site_all(s).name);
    if(ln>=30); file_tmp=strcat(site_all(s).name(1:12),yyyy,doy,site_all(s).name(20:38),'.amb1'); end
    if(ln<=30); file_tmp=strcat(lower(site_all(s).name(1:4)),doy,'0.',yyyy(3:4),'o','.amb1'); end
    file_tmp_t=strcat(dir,file_tmp);
    fid=fopen(file_tmp_t,'r');
    if(fid>0)
        site_num_ud=site_num_ud+1;
        site_all_ud(site_num_ud).name=file_tmp;
    else
        continue;
    end
    fclose(fid);
    % output new sitelist
    if(ln>=30); file_tmp=strcat(site_all(s).name(1:12),yyyy,doy,site_all(s).name(20:38)); end
    if(ln<=30); file_tmp=strcat(lower(site_all(s).name(1:4)),doy,'0.',yyyy(3:4),'o'); end
    fprintf(fout,'%s\n',file_tmp);
end
fclose(fout);