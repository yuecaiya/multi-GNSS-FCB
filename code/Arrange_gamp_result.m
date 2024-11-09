function [gps_Famb]=Arrange_gamp_result(yyyy,doy,site_amb,site_num_amb,flag,cfg)
%% Data format preprocessing for ambiguity, satellite azimuth and cycle slip
% args:
%     site_amb    : ambiguity, satellite azimuth and cycle slip data
%     site_num_amb: site number
%     cfg         : configuration information
%     flag        : indicator of ambiguity, satellite azimuth and cycle slip data
% made by Caiya Yue @ CUMTB and CASM
% ----
%% Reorganize floating ambiguity  with a sampling interval of 30 seconds
gps_Famb.name=zeros(site_num_amb,1);
nprn=cfg.NPRN;
rm_prn=str2num(cfg.rem_prn);
lnrm=length(rm_prn);
inter=cfg.inter;
nline=86400/inter;
i=0;
while(1)        
        i=i+1;
        if(i>site_num_amb)
            i=i-1;
            break;  
        end
        gps_Famb(i).name=site_amb(i).name;
        gps_Famb(i).amb=zeros(nline+1,nprn+1); % initialize preprocessed array matrix
        % Data is stored starting from the second row
        for ln=2:nprn+1
           gps_Famb(i).amb(1,ln)=ln-1+100; % first "1" in "101" represents GPS, and "01" represents satellite PRN
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop reading data
        dj=0; 
        k=0;
        for ep=2:nline+1
           k=k+1;
           if(site_amb(i).gamb(k,1)==dj)
               gps_Famb(i).amb(ep,1)=dj;
               for p=2:nprn+1
                  if(site_amb(i).gamb(k,p)==99999.000); site_amb(i).gamb(k,p)=0; end
                  gps_Famb(i).amb(ep,p)=site_amb(i).gamb(k,p);
                  % not estimating BDS-2 GEO satellites FCB
                  if(cfg.NavSystem==3&&p<=6&&cfg.geo_flag==0) 
                      gps_Famb(i).amb(ep,p)=0; 
                  end 
                  % excluded satellites
                  for rmn=1:lnrm
                     if(rm_prn(rmn)==p-1); gps_Famb(i).amb(ep,p)=0; end
                  end
               end               
           else
               gps_Famb(i).amb(ep,1)=dj;
               k=k-1;
           end
           dj=dj+1;
        end        
end
%%
if(strcmp(flag,'L3_ifb'))
    site_fcb=gps_Famb;
    sf=strcat(cfg.out_dir,'\Arra_',flag,'_',yyyy,doy);
    save(sf,'site_fcb');
elseif(strcmp(flag,'elev'))
    site_elev=gps_Famb;
    sf=strcat(cfg.out_dir,'\Arra_',flag,'_',yyyy,doy);
    save(sf,'site_elev');
elseif(strcmp(flag,'amb1'))
    site_fcb1=gps_Famb;
    sf=strcat(cfg.out_dir,'\Arra_',flag,'_',yyyy,doy);
    save(sf,'site_fcb1');
elseif(strcmp(flag,'amb2'))
    site_fcb2=gps_Famb;
    sf=strcat(cfg.out_dir,'\Arra_',flag,'_',yyyy,doy);
    save(sf,'site_fcb2');
elseif(strcmp(flag,'amb3'))
    site_fcb3=gps_Famb;
    sf=strcat(cfg.out_dir,'\Arra_',flag,'_',yyyy,doy);
    save(sf,'site_fcb3');
elseif(strcmp(flag,'cslip'))
    site_slip=gps_Famb;
    sf=strcat(cfg.out_dir,'\Arra_',flag,'_',yyyy,doy);
    save(sf,'site_slip');
end