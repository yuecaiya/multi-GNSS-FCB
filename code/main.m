clc;
clear;
format long
%% FCB/UPD estimation softer
% FCB/UPD estimation on both raw frequency and combined frequency
% Default system: GPS(G), Galileo(E), BDS-2(C2) and BDS-3(C3)
% Default frequency: G(L1 L2 L3), E(E1 E5a E5b), C2/3(B1I B3I B2I/B2a)
% made by Caiya Yue @ China University of Mining & Technology,Beijing (CUMTB)
%                   @ Chinese Academy of Surveying and Mapping (CASM)
%                   @ Liaocheng University (LCU)
% ----
%
%% Common parameter settings
data_dir='..\data\2023038_E\'; 
yyyy=2023;
doy=038;
NavSystem=2;   % 1:GPS   2:Galileo (GAL)   3:BDS-2 (BD2)   4:BDS-3 (BD3)
fprintf('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n');
if(NavSystem==1); fprintf('>>>> GPS UPD/FCB is being estimated >>>>\n'); end
if(NavSystem==2); fprintf('>>>> GAL UPD/FCB is being estimated >>>>\n'); end
if(NavSystem==3); fprintf('>>>> BD2 UPD/FCB is being estimated >>>>\n'); end
if(NavSystem==4); fprintf('>>>> BD3 UPD/FCB is being estimated >>>>\n'); end
fprintf('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n');
fprintf('\n');
%% Read configuration file
fprintf('** Read configuration file');
[cfg]=UD_read_configuration(data_dir,yyyy,doy,NavSystem);
fprintf(' ---Normal end\n');
fprintf('\n');
%% Create solution directory and open output file
mkdir(cfg.out_dir);
file_log=strcat('./',cfg.out_dir,'/sol_log.txt');
fid=fopen(file_log,'w');
fclose(fid);
diary(file_log);
diary on;
%% Read Station List - All Stations
fprintf('** Reading sitelist');
yyyy=num2str(cfg.yyyy);
doy=num2str(cfg.doy);
dir=cfg.dir;
file_site=strcat('../data/','sitename_all.txt');
site_num=0;
[site_all,site_num]=Read_Site(file_site);
fprintf(' ---Normal end\n');
fprintf('\n');
%% Update by reading the calculated station file
fprintf('** Updating sitelist according to observation');
if(cfg.upd_sitelist==1)
     [site_all_ud,site_num_ud]=Read_site_update(site_all,site_num,cfg);
end
[site_all,site_num]=Read_Site(strcat(cfg.out_dir,'/site_name_update.txt'));
fprintf(' ---Normal end\n');
fprintf('\n');
%% Excluding stations with frequent cycle jumps
fprintf('** Updating sitelist according to cycle slip\n');
sdoy=num2str(cfg.doy);
[site_slip0,site_num_slip0]=Read_gamp_result_cslip(yyyy,sdoy,dir,site_all,site_num,'cslip',cfg);
[site_all,site_num]=Read_Site(strcat(cfg.out_dir,'/site_name_update_cslip.txt'));
fprintf(' ---Normal end\n');
fprintf('\n');
%% Read the original frequency floating-point ambiguity, in meters
fprintf('** Reading float ambiguity\n');
amb1=cfg.amb(1:4);
amb2=cfg.amb(6:9);
amb3=cfg.amb(11:14);
sdoy=num2str(cfg.doy);
if(cfg.doy<10); sdoy=strcat('00',sdoy); end
if(cfg.doy<100 && cfg.doy>9); sdoy=strcat('0',sdoy); end
ra1=strcat(cfg.out_dir,'\Read_amb1_',yyyy,sdoy,'.mat');
ra2=strcat(cfg.out_dir,'\Read_amb2_',yyyy,sdoy,'.mat');
ra3=strcat(cfg.out_dir,'\Read_amb3_',yyyy,sdoy,'.mat');
% amb1
if ~exist(ra1,'file')
    [site_fcb1,site_num1_fcb]=Read_gamp_result(yyyy,sdoy,dir,site_all,site_num,amb1,cfg);
else
    load(ra1);
end
% amb2
if ~exist(ra2,'file')
    [site_fcb2,site_num2_fcb]=Read_gamp_result(yyyy,sdoy,dir,site_all,site_num,amb2,cfg);
else
    load(ra2);
end
% amb3
if ~exist(ra3,'file')
    [site_fcb3,site_num3_fcb]=Read_gamp_result(yyyy,sdoy,dir,site_all,site_num,amb3,cfg);
else
    load(ra3);
end
if(site_num1_fcb~=site_num2_fcb||site_num1_fcb~=site_num3_fcb); stop; end
site_num_fcb=site_num1_fcb;
fprintf(' ---Normal end\n');
fprintf('\n');
%% Data format preprocessing for float ambiguity
fprintf('** Arranging float ambiguity\n');
ra1=strcat(cfg.out_dir,'\Arra_amb1_',yyyy,sdoy,'.mat');
ra2=strcat(cfg.out_dir,'\Arra_amb2_',yyyy,sdoy,'.mat');
ra3=strcat(cfg.out_dir,'\Arra_amb3_',yyyy,sdoy,'.mat');
% amb1
if ~exist(ra1,'file')
    [gps_fcb1]=Arrange_gamp_result(yyyy,sdoy,site_fcb1,site_num_fcb,amb1,cfg); % --单位：m
else
    load(ra1);
    gps_fcb1=site_fcb1;
end
% amb2
if ~exist(ra2,'file')
    [gps_fcb2]=Arrange_gamp_result(yyyy,sdoy,site_fcb2,site_num_fcb,amb2,cfg);
else
    load(ra2);
    gps_fcb2=site_fcb2;
end
% amb3
if ~exist(ra3,'file')
    [gps_fcb3]=Arrange_gamp_result(yyyy,sdoy,site_fcb3,site_num_fcb,amb3,cfg);
else
    load(ra3);
    gps_fcb3=site_fcb3;
end
fprintf(' ---Normal end\n');
fprintf('\n');
%% Read satellite altitude angle in degrees
fprintf('** Reading elevation\n');
rae=strcat(cfg.out_dir,'\Read_elev_',yyyy,sdoy,'.mat');
if ~exist(rae,'file')
    [site_elev,site_num_elev]=Read_gamp_result(yyyy,sdoy,dir,site_all,site_num,'elev',cfg);
else
    load(rae);
end
fprintf(' ---Normal end\n');
fprintf('\n');
%% Data format preprocessing
fprintf('** Arranging elevation\n');
rae=strcat(cfg.out_dir,'\Arra_elev_',yyyy,sdoy,'.mat');
if ~exist(rae,'file')
    [gps_elev]=Arrange_gamp_result(yyyy,sdoy,site_elev,site_num_elev,'elev',cfg);
else
    load(rae);
    gps_elev=site_elev;
end
fprintf(' ---Normal end\n');
fprintf('\n');
%% Read cycle slip information
fprintf('** Reading cslip\n');
rac=strcat(cfg.out_dir,'\Read_cslip_',yyyy,sdoy,'.mat');
if ~exist(rac,'file')
    [site_slip,site_num_slip]=Read_gamp_result(yyyy,sdoy,dir,site_all,site_num,'cslip',cfg);
else
    load(rac);    
end
fprintf(' ---Normal end\n');
fprintf('\n');
%% Data format preprocessing
fprintf('** Arranging cslip\n');
rac=strcat(cfg.out_dir,'\Arra_cslip_',yyyy,sdoy,'.mat');
if ~exist(rac,'file')
    [gps_slip]=Arrange_gamp_result(yyyy,sdoy,site_slip,site_num_slip,'cslip',cfg);
else
    load(rac);
    gps_slip=site_slip;
end
fprintf(' ---Normal end\n');
fprintf('\n');
%% BDS-3 B1C and B2a frequency band settings
gps_fcbt=0;
Freq=0;
if(cfg.NavSystem==4 && cfg.new_ban==1) % BDS-3 B1C
    gps_fcbt=gps_fcb1;
    gps_fcb1=gps_fcb2;
    gps_fcb2=gps_fcbt;
    Freq=cfg.Freq(8);
    cfg.Freq(8)=cfg.Freq(10);
    cfg.Freq(10)=Freq;
end
if(cfg.NavSystem==4 && cfg.new_ban==2) % BDS-3 B2a
    cfg.Freq(11)=cfg.Freq(12);
end
%% Linear combination {ultra wide lane (UWL), wide lane (WL) and narrow lane (NL)} -- (cycle)
fprintf('** Executing UD_N1N2N3_amb');
[Namb23,Namb12,Namb43,Pelev23,Pelev12,Pelev43]=UD_N1N2N3_amb(gps_fcb1,gps_fcb2,gps_fcb3,gps_elev,site_num_fcb,cfg);
fprintf(' ---Normal end\n');
fprintf('\n');

%% Merge within arc segments based on variance and cycle slip information
fprintf('** Executing UD_combine_amb_arc\n');
[Tamb23,Tamb12,Tamb43,Aelev23,Aelev12,Aelev43]=UD_combine_amb_arc(Namb23,Namb12,Namb43,Pelev23,Pelev12,Pelev43,gps_slip,site_num_fcb,cfg);
fprintf('   ----Normal end\n');
fprintf('\n');

%% Extract the merged floating ambiguity fractional part
fprintf('** Executing UD_combine_Fra_amb');
[FTamb23,FTamb12,FTamb43]=UD_combine_Fra_amb(Tamb23,Tamb12,Tamb43,site_num_fcb,cfg);
fprintf(' ---Normal end\n');
fprintf('\n');

%% Reference stations and satellite adaptive judgment
fprintf('** Executing UD_total_sat_num');
% if(cfg.ref_prn==0)
    [num_sats]=UD_total_sat_num(gps_fcb3,site_num_fcb);
    [data,ref_prn]=max(num_sats);
% end
% if(cfg.ref_prn==1)
    [num_recs]=UD_total_rec_num(gps_fcb3,site_num_fcb);
    [data,ref_rec]=max(num_recs);
% end
% if(cfg.ref_rec>0); ref_rec=cfg.ref_prn; end
if(cfg.ref_prn==0); ref_prn_sat=ref_prn; end
if(cfg.ref_prn==1); ref_prn_sat=ref_rec; end
fprintf(' ---Normal end\n');fprintf('\n');
%% two calculation methods for UPD/FCB
fprintf('** Executing UD_estimate_upd\n');
save(strcat(cfg.out_dir,'\FTamb23'),'FTamb23');
save(strcat(cfg.out_dir,'\FTamb12'),'FTamb12');
save(strcat(cfg.out_dir,'\FTamb43'),'FTamb43');
save(strcat(cfg.out_dir,'\Aelev23'),'Aelev23');
save(strcat(cfg.out_dir,'\Aelev12'),'Aelev12');
save(strcat(cfg.out_dir,'\Aelev43'),'Aelev43');
load(strcat(cfg.out_dir,'\FTamb23.mat'));
load(strcat(cfg.out_dir,'\FTamb12.mat'));
load(strcat(cfg.out_dir,'\FTamb43.mat'));
load(strcat(cfg.out_dir,'\Aelev23.mat'));
load(strcat(cfg.out_dir,'\Aelev12.mat'));
load(strcat(cfg.out_dir,'\Aelev43.mat'));
%% Estimation on raw frequency FCB/UPD
if(cfg.FCB_MOD==1)
    % FCB/UPD estimation
    [final_upd,QSupd]=UD_estimate_upd(yyyy,doy,FTamb23,FTamb12,FTamb43,Aelev23,Aelev12,Aelev43,site_num_fcb,ref_prn_sat,cfg,ref_prn); %---- 原始方法
    fprintf('       ----Normal end\n');
    fprintf('\n');
    % output and plot FCB/UPD 
    fprintf('** Outputing the raw and LC UPD\n');
    [stat]=out_upd_file_nl(yyyy,doy,final_upd,cfg,QSupd,ref_prn);
    fprintf(' ----Normal end\n');
    fprintf('\n');
end
%% UWL-WL-NL FCB/UPD estimation
if(cfg.FCB_MOD==2)
    ref_prn_sat=ref_prn;
    % FCB/UPD estimation
    [final_upd_u,QSupd_fix_amb_u]=UD_estimate_upd_uwn(yyyy,sdoy,Tamb23,Tamb12,Tamb43,FTamb23,FTamb12,FTamb43,Aelev23,Aelev12,Aelev43,site_num_fcb,ref_prn_sat,cfg,ref_prn,0,'u'); %---- UWL
    [final_upd_w,QSupd_fix_amb_w]=UD_estimate_upd_uwn(yyyy,sdoy,Tamb23,Tamb12,Tamb43,FTamb23,FTamb12,FTamb43,Aelev23,Aelev12,Aelev43,site_num_fcb,ref_prn_sat,cfg,ref_prn,0,'w'); %---- WL
    % Single day UPD/FCB synthesis
%     [final_upd_w_day]=UD_compose_day_FCB(QSupd_fix_amb_w,cfg);
%     QSupd_fix_amb_w1=final_upd_w_day;
    [final_upd_n,QSupd_fix_amb_n]=UD_estimate_upd_uwn(yyyy,sdoy,Tamb23,Tamb12,Tamb43,FTamb23,FTamb12,FTamb43,Aelev23,Aelev12,Aelev43,site_num_fcb,ref_prn_sat,cfg,ref_prn,QSupd_fix_amb_w,'n'); %---- NL
    fprintf('       ----Normal end\n');
    fprintf('\n');
    % FCB/UPD output and plot
    [a,ln]=size(final_upd_u);
    lc=length(final_upd_u(1).upd(:,1));
    for e=1:ln
        final_upd(e).upd(:,1)=final_upd_u(e).upd(:,1); % UWL FCB
        final_upd(e).upd(:,2)=final_upd_w(e).upd(:,1); %  WL FCB
        final_upd(e).upd(:,3)=final_upd_n(e).upd(:,1); %  NL FCB
    end    
    fprintf('** Outputing the raw and LC UPD\n');
    QSupd=0;
    [stat]=out_upd_file_nl_uwn(yyyy,sdoy,final_upd,cfg,QSupd,ref_prn);
    fprintf(' ----Normal end\n');
    fprintf('\n');
end
%% Output the content of the command window into a file
diary off;
