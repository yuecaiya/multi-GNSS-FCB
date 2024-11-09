function [cfg]=UD_read_configuration(data_dir,yyyy,doy,NavSystem)
%% Read configuration file
%
% ---
%% Constant section
CLIGHT    = 299792458.0;           % speed of light (m/s)
% GPS
FREQ1_G     = 1.57542E9;           % L1/E1  frequency (Hz)
FREQ2_G     = 1.22760E9;           % L2     frequency (Hz)
FREQ3_G     = 1.17645E9;           % L5/E5a frequency (Hz)
% GLONASS
FREQ1_R = 1.60200E9;               % GLONASS G1 base frequency (Hz)
DFRQ1_R = 0.56250E6;               % GLONASS G1 bias frequency (Hz/n)
FREQ2_R = 1.24600E9;               % GLONASS G2 base frequency (Hz)
DFRQ2_R = 0.43750E6;               % GLONASS G2 bias frequency (Hz/n)
FREQ3_R = 1.202025E9;              % GLONASS G3 frequency (Hz)
% Galileo
FREQ1_E     = 1.57542E9;           % E1  frequency (Hz)
FREQ2_E     = 1.17645E9;           % E5a frequency (Hz)
FREQ3_E     = 1.20714E9;           % E5b frequency (Hz)
FREQ4_E     = 2.492028E9;          % S   frequency (Hz)
% BDS
FREQ1_CMP = 1.561098E9;            % BeiDou B1 frequency (Hz)
FREQ2_CMP = 1.20714E9;             % BeiDou B2 frequency (Hz)
FREQ3_CMP = 1.26852E9;             % BeiDou B3 frequency (Hz)
FREQ4_CMP = 1.57542E9;             % BeiDou B1C frequency (Hz)
FREQ5_CMP = 1.17645E9;             % BeiDou B2a frequency (Hz)
FREQ6_CMP = 1.20714E9;             % BeiDou B2b frequency (Hz)
cfg.Freq=[FREQ1_G FREQ2_G FREQ3_G FREQ1_E FREQ2_E FREQ3_E FREQ4_E FREQ1_CMP FREQ2_CMP FREQ3_CMP FREQ4_CMP FREQ5_CMP FREQ6_CMP];
%% Configuration option
cfg.elev=15;                   % satellite cut-off altitude angle£¨¡ã£©
cfg.inter=30;                  % raw data sampling interval£¨s£©
cfg.span_t=30;                 % UPD/FCB output sampling interval£¨min£©
cfg.com3='4 -3 0';             % NL combination coefficients for estimating raw frequency FCB/UPD
cfg.amb='amb1 amb2 amb3';      % floating ambiguity
cfg.plot_flag=1;               % 0:not ploting;  1: plot FCB/UPD
cfg.upd_sitelist=1;            % 0:Do not update station list;  1:update station list
cfg.amb_thh='0.15 0.2 0.25';   % median exclusion threshold for UWL, WL and NL
cfg.sigma='3.0 3.0 2.0';       % threshod of n times sigma for  UWL, WL and NL
cfg.fix_amb='0.25 0.25 0.9';   % adjustment weight reduction threshold (cycle), fixed ambiguity threshold (cycle), adjustment weight reduction value
cfg.weight=1;                  % 0:equal-weighted;  1:STD weight calculation
cfg.original_upd=0;            % 0: first epoch£»  1£ºeach epoch
cfg.cslip_thr=0.90;            % Sort the stations based on the total number of cycle slips and retain stations within this rate
cfg.dir=data_dir; 
cfg.yyyy=yyyy;
cfg.doy=doy;
cfg.FCB_MOD=2;                   % 1:raw FCB   2:UWL-WL-NL FCB                      
cfg.NavSystem=NavSystem;
    if(cfg.NavSystem==1); cfg.NPRN=32; end  % number of satellites                                  %
    if(cfg.NavSystem==2); cfg.NPRN=40; end                                                          %
    if(cfg.NavSystem==3); cfg.NPRN=60; end                                                          %
    if(cfg.NavSystem==4); cfg.NPRN=60; end                                                          %
    if(cfg.NavSystem==4)                                                                            %
        cfg.new_ban=2;           % 1:B1C  2:B2a                                                     %
    end                                                                                             %
    %---------------------------------Arc segment merging control strategy---------------------------
    if(cfg.NavSystem==1); cfg.comb_time=90; end  % minimum arc segment for merging £¨min£©
    if(cfg.NavSystem==2); cfg.comb_time=120; end                                                     %
    if(cfg.NavSystem==3); cfg.comb_time=120; end                                                     %
    if(cfg.NavSystem==4); cfg.comb_time=120; end                                                     %
    %---------------------------------Excluded satellites---------------------------------------------
    if(cfg.NavSystem==1); cfg.rem_prn='00 00'; end   
    if(cfg.NavSystem==2); cfg.rem_prn='00 00'; end                                                      
    if(cfg.NavSystem==3); cfg.rem_prn='00 00'; end                                                      
    if(cfg.NavSystem==4); cfg.rem_prn='00 36 41 42'; end 
cfg.comb_arc=30; % merge the data from the last 30 minutes as the ambiguity parameter
cfg.nprn_C2=16;  % number of BDS-2 satellites                                                
cfg.geo_flag=0;  % 0£ºdo not estimate GEO satellite;  1:estimate GEO satellite                              
cfg.ref_prn=0;   % 0£ºmax satellite  1£ºmax site                                            
%*------------------------------Output section settings----------------------------------------------%
if(cfg.doy<100 && cfg.doy>9); s_doy=strcat('0',num2str(cfg.doy)); end
if(cfg.doy<10); s_doy=strcat('00',num2str(cfg.doy)); end
if(cfg.doy>99); s_doy=num2str(cfg.doy); end
if(cfg.NavSystem==1); out_dir=strcat('../result/result_G_',num2str(cfg.yyyy),s_doy); end                                         
if(cfg.NavSystem==2); out_dir=strcat('../result/result_E_',num2str(cfg.yyyy),s_doy); end                                                      
if(cfg.NavSystem==3); out_dir=strcat('../result/result_C2_',num2str(cfg.yyyy),s_doy); end                                                      
if(cfg.NavSystem==4); out_dir=strcat('../result/result_C3_',num2str(cfg.yyyy),s_doy); end                                                      
cfg.out_dir=out_dir;
end