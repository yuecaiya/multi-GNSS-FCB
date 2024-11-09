function [site_amb,site_num_amb]=Read_gamp_result(yyyy,doy,dir,site_all,site_num,flag,cfg) 
%% Read ambiguity, cycle slip and satellite azimth
% args:
%     dir: file catalogue of ambiguity, cycle slip and satellite azimth
%     site_all: sitelist
%     site_num: site number
%     cfg     : configuration information
% made by Caiya Yue @ CUMTB and CASM
% ----
%% interpreting control information
if(cfg.NavSystem==1); gps_ln=357; end % gps 
if(cfg.NavSystem==2); gps_ln=437; end % Galileo 
if(cfg.NavSystem==3); gps_ln=637; end % BDS2 
if(cfg.NavSystem==4); gps_ln=637; end % BDS2 
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
        del_site=strcat(site_all(ii).name,' noexist')
        fprintf(fid_del_sta,'%1c',del_site);
        fprintf(fid_del_sta,'\n',del_site);
        site_num_amb=site_num_amb-1;
        fclose(fid);
        continue;
    end    
    % read data  
    i=i+1;
    site_amb(i).name = site_all(ii).name;
    site_amb(i).gamb = zeros(nline,nprn+1); % the first row and the first column are satellites and time, respectively
    j=0;
    while(1)
        strNlamb=fgets(fid);
        if(strNlamb == -1)
            break;
        end
        ln=length(strNlamb);
        if(length(strNlamb)<=gps_ln)
            if(rem((str2num(strNlamb(12:13))*3600+str2num(strNlamb(15:16))*60+str2num(strNlamb(18:19))),30)~=0); continue; end
            j=(str2num(strNlamb(12:13))*3600+str2num(strNlamb(15:16))*60+str2num(strNlamb(18:19)))/30;
            j=j+1;        
            site_amb(i).gamb(j,1)=(str2num(strNlamb(12:13))*3600+str2num(strNlamb(15:16))*60+str2num(strNlamb(18:19)))/30; % (second)
            site_amb(i).gamb(j,2:nprn+1)=str2num(strNlamb(36:ln))';
            if(cfg.NavSystem==3); site_amb(i).gamb(j,nprnc+2:nprn+1)=99999.000; end % BDS-2
            if(cfg.NavSystem==4); site_amb(i).gamb(j,2:nprnc+1)=99999.000; end % BDS-3
        end            
    end
    fclose(fid);
end
fclose(fid_del_sta);

if(strcmp(flag,'L3_ifb'))
    site_fcb=site_amb;
    site_num_fcb=site_num_amb;
    sf=strcat(cfg.out_dir,'\Read_',flag,'_',yyyy,doy);
    save(sf,'site_fcb','site_num_fcb');
elseif(strcmp(flag,'elev'))
    site_elev=site_amb;
    site_num_elev=site_num_amb;
    sf=strcat(cfg.out_dir,'\Read_',flag,'_',yyyy,doy);
    save(sf,'site_elev','site_num_elev');
elseif(strcmp(flag,'amb1'))
    site_fcb1=site_amb;
    site_num1_fcb=site_num_amb;
    sf=strcat(cfg.out_dir,'\Read_',flag,'_',yyyy,doy);
    save(sf,'site_fcb1','site_num1_fcb');
elseif(strcmp(flag,'amb2'))
    site_fcb2=site_amb;
    site_num2_fcb=site_num_amb;
    sf=strcat(cfg.out_dir,'\Read_',flag,'_',yyyy,doy);
    save(sf,'site_fcb2','site_num2_fcb');
elseif(strcmp(flag,'amb3'))
    site_fcb3=site_amb;
    site_num3_fcb=site_num_amb;
    sf=strcat(cfg.out_dir,'\Read_',flag,'_',yyyy,doy);
    save(sf,'site_fcb3','site_num3_fcb');
elseif(strcmp(flag,'cslip'))
    site_slip=site_amb;
    site_num_slip=site_num_amb;
    sf=strcat(cfg.out_dir,'\Read_',flag,'_',yyyy,doy);
    save(sf,'site_slip','site_num_slip');
end
