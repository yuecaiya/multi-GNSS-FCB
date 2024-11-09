function [v,B,v0,B0,P,prn_flag,rec_prn_flag,rm_site_n,rm_site]=make_coe_matrix_uwn1(FTamb23,FTamb12,FTamb43,Aelev23,Aelev12,Aelev43,site_num_fcb,ref_prn,cfg,ep,upd_wm,uwn_flag,gps_upd)
%% Constructing coefficient matrix
% args:
%     FTamb23 : UWL combination ambiguity
%     FTamb12 :  WL combination ambiguity
%     FTamb43 :  NL combination ambiguity
%     Aelev23 : weight for UWL
%     Aelev12 : weight for  WL
%     Aelev43 : weight for  NL
%     ref_prn : reference satellite or receiver
%     site_num_fcb: site number
%     cfg         : configuration information
% return:
%     v : residual vector 
%     B : coefficient matrix
%     P : weight matrix
% made by Caiya Yue @ CUMTB and CASM
% ----
%% interpreting control information
CLIGHT    = 299792458.0;         % speed of light (m/s) */
% GPS
if(cfg.NavSystem==1)
    f1=cfg.Freq(1);
    f2=cfg.Freq(2);
    coef=f2/(f1-f2);
    lam_nl=CLIGHT/(f1+f2);
end
% Galileo
if(cfg.NavSystem==2)
    f1=cfg.Freq(4);
    f2=cfg.Freq(5);
    coef=f2/(f1-f2);
    lam_nl=CLIGHT/(f1+f2);
end
% BDS2
if(cfg.NavSystem==3)
    f1=cfg.Freq(8);
    f2=cfg.Freq(10);
    coef=f2/(f1-f2);
    lam_nl=CLIGHT/(f1+f2);
end
% BDS3
if(cfg.NavSystem==4)
    f1=cfg.Freq(8);
    f2=cfg.Freq(10);
    coef=f2/(f1-f2);
    lam_nl=CLIGHT/(f1+f2);
end
sigma0=0.001;
lc=length(FTamb23(1).amb(1,:));
lam_wl=CLIGHT/(f1-f2);
%% perform the first coefficient design for all stations
% BT¡¾column: UWL WL NL¡¿
%   ¡¾raw   :site satellite¡¿
bt=zeros(site_num_fcb*lc+1,site_num_fcb+lc);
vt=zeros(site_num_fcb*lc+1,1);
pt=zeros(site_num_fcb*lc+1,site_num_fcb*lc+1);

TMPamb=0;
if(strcmp(uwn_flag,'u'))
    TMPamb=FTamb23;
    TMPele=Aelev23;
end
if(strcmp(uwn_flag,'w'))
    TMPamb=FTamb12;
    TMPele=Aelev12;
end
if(strcmp(uwn_flag,'n'))
    TMPamb=FTamb43;
    TMPele=Aelev43;
end
    
rm_site=zeros(site_num_fcb,1);
rm_site_n=0;
for site=1:site_num_fcb
    flag_tm1=0;
    for prn=1:lc      
        if(TMPamb(site).amb(ep,prn)~=0)
            if(strcmp(uwn_flag,'u') || strcmp(uwn_flag,'w')); vt((site-1)*lc+(prn-1)+1,1)=TMPamb(site).amb(ep,prn); end
            if(strcmp(uwn_flag,'n'))
                [upd_w_m]=UD_estimate_upd_tq_fix_amb(TMPamb,upd_wm,ep,prn,site,cfg);
                if(upd_w_m==9999.999); continue; end
                 vt((site-1)*lc+(prn-1)+1,1)=TMPamb(site).amb(ep,prn)/lam_nl-coef*upd_w_m;
%                  vt((site-1)*lc+(prn-1)+1,1)=(TMPamb(site).amb(ep,prn)-coef*upd_w_m*lam_wl)/lam_nl;
                 vt((site-1)*lc+(prn-1)+1,1)=vt((site-1)*lc+(prn-1)+1,1)-fix(vt((site-1)*lc+(prn-1)+1,1));
            end
            bt((site-1)*lc+(prn-1)+1,(site-1)+1)=1;
            bt((site-1)*lc+(prn-1)+1,site_num_fcb+(prn-1)+1)=-1;
            pt((site-1)*lc+(prn-1)+1,(site-1)*lc+(prn-1)+1)=sigma0/TMPele(site).elev(ep,prn);
            flag_tm1=1;
        end
    end
    if(flag_tm1==0 )
        rm_site_n=rm_site_n+1;
        rm_site(rm_site_n)=site;
    end
end
%% update FTamb23 FTamb12,FTamb43,Aelev23,Aelev12,Aelev43 based on the excluded stations
nsite=0;
for site=1:site_num_fcb
    site_tc_flag=0;
    for ns=1:rm_site_n
       if(site==rm_site(ns)); site_tc_flag=1; end
    end
    if(site_tc_flag==1); continue; end
    nsite=nsite+1;
    UDTMPamb(nsite)=TMPamb(site);
    UDTMPele(nsite)=TMPele(site);
end
%% implement a second coefficient design based on the excluded stations
bt=zeros(nsite*lc+1,nsite+lc);
vt=zeros(nsite*lc+1,1);
pt=zeros(nsite*lc+1,nsite*lc+1);
k=0;
rec_prn_flag=[];
for site=1:nsite
    for prn=1:lc  
        if(UDTMPamb(site).amb(ep,prn)~=0)    
            if(strcmp(uwn_flag,'u') || strcmp(uwn_flag,'w')); vt((site-1)*lc+(prn-1)+1,1)=UDTMPamb(site).amb(ep,prn); end
            if(strcmp(uwn_flag,'n'))
                [upd_w_m]=UD_estimate_upd_tq_fix_amb(UDTMPamb,upd_wm,ep,prn,site,cfg);
                if(upd_w_m==9999.999); continue; end
                 vt((site-1)*lc+(prn-1)+1,1)=UDTMPamb(site).amb(ep,prn)/lam_nl-coef*upd_w_m;
%                  vt((site-1)*lc+(prn-1)+1,1)=(TMPamb(site).amb(ep,prn)-coef*upd_w_m*lam_wl)/lam_nl;
                 vt((site-1)*lc+(prn-1)+1,1)=vt((site-1)*lc+(prn-1)+1,1)-fix(vt((site-1)*lc+(prn-1)+1,1));
            end
            bt((site-1)*lc+(prn-1)+1,(site-1)+1)=1;
            bt((site-1)*lc+(prn-1)+1,nsite+(prn-1)+1)=-1;
            pt((site-1)*lc+(prn-1)+1,(site-1)*lc+(prn-1)+1)=sigma0/UDTMPele(site).elev(ep,prn);
            k=k+1;
            rec_prn_flag(k,1)=site;
            rec_prn_flag(k,2)=prn;
        end       
    end
end
% if(cfg.ref_prn==0) % % reference satellite
    vt(nsite*lc+1,1)=0;
    pt(nsite*lc+1,nsite*lc+1)=100000;
    bt(nsite*lc+1,nsite+(ref_prn-1)+1)=-1;
% end
% if(cfg.ref_prn==1) % reference receiver
%     vt(nsite*lc+1,1)=0;
%     pt(nsite*lc+1,nsite*lc+1)=100000;
%     bt(nsite*lc+1,ref_prn)=1;
% end
B0=bt;
v0=vt;
%% remove rows with all zeros
all_ln=nsite*lc+1;
btt=[];
vtt=[];
ptt=[];
k=0;
for n=1:all_ln
    if(sum(abs(bt(n,:)))~=0)
       k=k+1;
       btt(k,:)=bt(n,:);
       vtt(k,1)=vt(n,1);
       ptt(k,k)=pt(n,n);
    end
end
%% remove columns with all zeros
all_lc=nsite+lc;
k=0;
tc_lc=zeros(all_lc,1);
for c=1:all_lc
    if(sum(abs(btt(:,c)))==0)
       k=k+1;
       tc_lc(k)=c;
    end
end
% UPD/FCB estimation flag
tq=0;
prn_flag=zeros(lc,1);
for c=1:k
    if(tc_lc(c)<=nsite)
        site=tc_lc(c);
        fprintf('Warning: The station of %s have no data!',UDFTamb23(site).name);
%         stop
    else
        prn=tc_lc(c)-nsite; % prn
        tc=tc_lc(c)-tq;
        btt(:,tc)=[];
        prn_flag(prn)=999.99;
        tq=tq+1;
    end
end
B=btt;
v=vtt;
P=ptt;
end