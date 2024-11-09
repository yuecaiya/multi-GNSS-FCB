function [site_amb,site_num_amb]=Read_gamp_result_cslip(yyyy,doy,dir,site_all,site_num,flag,cfg) 
%% Excluding stations with frequent cycle jumps
% args:
%     dir: cycle file catalogue
%     site_all: sitelist
%     site_num: site number
%     cfg     : configuration information
% made by Caiya Yue @ CUMTB and CASM
% ----
%% interpreting control information
if(cfg.NavSystem==1); gps_ln=357; end % minimum number of characters in a row for input file
if(cfg.NavSystem==2); gps_ln=437; end % 
if(cfg.NavSystem==3); gps_ln=437; end % 
if(cfg.NavSystem==4); gps_ln=437; end % 
nprn=cfg.NPRN;
nprnc=cfg.nprn_C2;
inter=cfg.inter;
nline=86400/inter;
%% body
del_sta=strcat(cfg.out_dir,'/noexist_',flag,'_sta.txt');
fid_del_sta=fopen(del_sta,'w');
site_num_amb=site_num; % count the remaining number of stations
i=0;
for ii=1:site_num
    file_nlamb=strcat(dir,site_all(ii).name,'.',flag); 
    exe_status=strcat('Reading F-amb file: ',file_nlamb,' Num. ',num2str(ii),'/',num2str(site_num));
    exe_s=strcat(site_all(ii).name,'.',flag);
    fprintf('    %s\n',exe_s);
    fid=fopen(file_nlamb,'r');    
    % output non-existent station
    if(fid <= 0)
        del_site=strcat(site_all(ii).name,' noexist');
        fprintf(fid_del_sta,'%1c',del_site);
        fprintf(fid_del_sta,'\n',del_site);
        site_num_amb=site_num_amb-1;
        fclose(fid);
        continue;
    end    
    % read data   
    kt=0;
    tmp_cs_tol=[]; 
    while(1)
        strNlamb=fgets(fid);
        if(strNlamb == -1)
            break;
        end               
        if(length(strNlamb)<=gps_ln)
            num_str=str2num(strNlamb);
            lnn=length(num_str);
            kt=kt+1;
            tmp_cs_tol(kt,1:lnn-8)=num_str(9:lnn); 
        end
    end
    if(sum(sum(tmp_cs_tol))==0); continue; end
    i=i+1;
    site_amb(i).cslip=sum(sum(tmp_cs_tol));    
    site_amb(i).name = site_all(ii).name;    
    fclose(fid);
end
fclose(fid_del_sta);
%% sort and output sitelist based on the number of cycle slips
tmp_cs_1=zeros(i,1);
for k=1:i
    tmp_cs_1(k)=site_amb(k).cslip;
end
[a,b]=sort(tmp_cs_1);
yz_cs=fix(i*cfg.cslip_thr);
yyyy=num2str(cfg.yyyy);
doy=num2str(cfg.doy);
if(cfg.doy<10); doy=strcat('00',doy); end
if(cfg.doy<100 && cfg.doy>9); doy=strcat('0',doy); end
% output sitelist
site_num_ud=0;
site_all_ud.name = zeros(300,1);
fout=fopen(strcat(cfg.out_dir,'/site_name_update_cslip.txt'),'w');
for kk=1:yz_cs
    lna=length(site_amb(b(kk)).name);
    fprintf(fout,'%s\n',site_amb(b(kk)).name(1:lna-1));
end
fclose(fout);
end
