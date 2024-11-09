function [Namb23,Namb12,Namb43,Pelev23,Pelev12,Pelev43]=UD_N1N2N3_amb(gps_fcb1,gps_fcb2,gps_fcb3,gps_elev,site_num_fcb,cfg)
%% Linear combination {ultra wide lane (UWL), wide lane (WL) and narrow lane (NL)} -- (cycle)
% args:
%     gps_fcb1    : ambiguity on L1 frequency
%     gps_fcb2    : ambiguity on L2 frequency
%     gps_fcb3    : ambiguity on L3 frequency
%     gps_elev    : satellite altitude angle
%     site_num_fcb: site number
%     cfg         : configuration information
% return:
%     Namb23 : UWL combination ambiguity
%     Namb12 :  WL combination ambiguity
%     Namb43 :  NL combination ambiguity
%     Pelev23: weight for UWL
%     Pelev12: weight for  WL
%     Pelev43: weight for  NL
% made by Caiya Yue @ CUMTB and CASM
% ----
%% interpreting control files
xs=str2num(cfg.com3);
x1=xs(1);
x2=xs(2);
x3=xs(3);
elev=cfg.elev;
%%
CLIGHT    = 299792458.0;         % speed of light (m/s) */
% GPS
if(cfg.NavSystem==1)
    lam1=CLIGHT/cfg.Freq(1);
    lam2=CLIGHT/cfg.Freq(2);
    lam3=CLIGHT/cfg.Freq(3);
    lc1=-lam2^2/(lam1^2-lam2^2);
    lc2= lam1^2/(lam1^2-lam2^2);
end
% Galileo
if(cfg.NavSystem==2)
    lam1=CLIGHT/cfg.Freq(4);
    lam2=CLIGHT/cfg.Freq(5);
    lam3=CLIGHT/cfg.Freq(6);
    lc1=-lam2^2/(lam1^2-lam2^2);
    lc2= lam1^2/(lam1^2-lam2^2);
end
% BDS2
if(cfg.NavSystem==3)
    lam1=CLIGHT/cfg.Freq(8);
    lam2=CLIGHT/cfg.Freq(10);
    lam3=CLIGHT/cfg.Freq(9); % b2i
    lc1=-lam2^2/(lam1^2-lam2^2);
    lc2= lam1^2/(lam1^2-lam2^2);
end
% BDS3
if(cfg.NavSystem==4)
    lam1=CLIGHT/cfg.Freq(8);
    lam2=CLIGHT/cfg.Freq(10);
    lam3=CLIGHT/cfg.Freq(11); % b1c or b2a
    lc1=-lam2^2/(lam1^2-lam2^2);
    lc2= lam1^2/(lam1^2-lam2^2);
end
%% 
for site=1:site_num_fcb
    ln=length(gps_fcb1(site).amb(:,1))-1;
    lc=length(gps_fcb1(site).amb(1,:))-1;
    Namb23(site).name=gps_fcb1(site).name;
    Namb12(site).name=gps_fcb2(site).name;
    Namb43(site).name=gps_fcb3(site).name;
    Namb23(site).amb=zeros(ln,lc);
    Namb12(site).amb=zeros(ln,lc);
    Namb43(site).amb=zeros(ln,lc);
    Pelev23(site).elev=zeros(ln,lc); %£¨0 1 -1£©
    Pelev12(site).elev=zeros(ln,lc); %£¨1 -1 0£©
    Pelev43(site).elev=zeros(ln,lc); %£¨4 -3 0£©
     %% constructing a linear combination of ambiguity
    for ep=1:ln  
       for prn=1:lc
          if(gps_elev(site).amb(ep+1,prn+1)>30)
              sigma_e=1;
          else
              el_arc=gps_elev(site).amb(ep+1,prn+1)*pi/180;
              sigma_e=2*sin(el_arc);
          end 
          % UWL
          if((gps_fcb2(site).amb(ep+1,prn+1)*gps_fcb3(site).amb(ep+1,prn+1)~=0) && gps_elev(site).amb(ep+1,prn+1)>elev) % (0 1 -1)
               Namb23(site).amb(ep,prn)=gps_fcb2(site).amb(ep+1,prn+1)/lam2-gps_fcb3(site).amb(ep+1,prn+1)/lam3;               
               Pelev23(site).elev(ep,prn)=((1/lam2)^2+(1/lam3)^2)*(1/sigma_e); % variance or weight
          end 
          % WL and free ionospheric combination ambiguity
          if((gps_fcb1(site).amb(ep+1,prn+1)*gps_fcb2(site).amb(ep+1,prn+1)~=0) && gps_elev(site).amb(ep+1,prn+1)>elev) %(1 -1 0)
              Namb12(site).amb(ep,prn)=gps_fcb1(site).amb(ep+1,prn+1)/lam1-gps_fcb2(site).amb(ep+1,prn+1)/lam2; % variance or weight
              Pelev12(site).elev(ep,prn)=((1/lam1)^2+(1/lam2)^2)*(1/sigma_e);
              % Raw FCB
              if(cfg.FCB_MOD==1)  % (4 -3 0)
                  Namb43(site).amb(ep,prn)=x1*gps_fcb1(site).amb(ep+1,prn+1)/lam1+x2*gps_fcb2(site).amb(ep+1,prn+1)/lam2+x3*gps_fcb3(site).amb(ep+1,prn+1)/lam3;
                  Pelev43(site).elev(ep,prn)=x1^2*(1/sigma_e)/lam1^2+x2^2*(1/sigma_e)/lam2^2+x3^2*(1/sigma_e)/lam3^2; % variance or weight
              end
              % UWL-WL-NL
              if(cfg.FCB_MOD==2) %free ionospheric combination
                    Namb43(site).amb(ep,prn)=lc1*gps_fcb1(site).amb(ep+1,prn+1)+lc2*gps_fcb2(site).amb(ep+1,prn+1);
                    Pelev43(site).elev(ep,prn)=lc1^2*(1/sigma_e)+lc2^2*(1/sigma_e); % variance or weight
              end
          end
       end       
    end 
end