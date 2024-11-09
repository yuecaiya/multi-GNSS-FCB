function [v,B,P,prn_flag,rec_flag,rm_site_n,rm_site]=make_coe_matrix(FTamb23,FTamb12,FTamb43,Aelev23,Aelev12,Aelev43,site_num_fcb,ref_prn,cfg,ep)
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
xs=str2num(cfg.com3); 
x1=xs(1);
x2=xs(2);
x3=xs(3);
sigma0=0.001;
lc=length(FTamb23(1).amb(1,:));
%% perform the first coefficient design for all stations
% BT¡¾column: UWL WL NL¡¿
%   ¡¾raw   :site satellite¡¿
bt=zeros(site_num_fcb*3*lc+3,site_num_fcb*3+lc*3);
vt=zeros(site_num_fcb*3*lc+3,1);
pt=zeros(site_num_fcb*3*lc+3,site_num_fcb*3*lc+3);
prn_flag=zeros(lc,3); % UPD/FCB estimation flag: L1 L2 L3
for prn=1:lc
   prn_flag(prn,1)=prn;
   prn_flag(prn,2)=prn;
   prn_flag(prn,3)=prn;
end
rm_site=zeros(site_num_fcb,1);
rm_site_n=0;
for site=1:site_num_fcb
    flag_tm1=0;
    flag_tm2=0;
    flag_tm3=0;
    for prn=1:lc  
        % ---UWL£¨0 1 -1£©        
        if(FTamb23(site).amb(ep,prn)~=0)
            vt((site-1)*lc*3+(prn-1)*3+1,1)=FTamb23(site).amb(ep,prn);
            bt((site-1)*lc*3+(prn-1)*3+1,(site-1)*3+1)=0;
            bt((site-1)*lc*3+(prn-1)*3+1,(site-1)*3+2)=1;
            bt((site-1)*lc*3+(prn-1)*3+1,(site-1)*3+3)=-1;
            bt((site-1)*lc*3+(prn-1)*3+1,site_num_fcb*3+(prn-1)*3+1)=0;
            bt((site-1)*lc*3+(prn-1)*3+1,site_num_fcb*3+(prn-1)*3+2)=-1;
            bt((site-1)*lc*3+(prn-1)*3+1,site_num_fcb*3+(prn-1)*3+3)=1;
            pt((site-1)*lc*3+(prn-1)*3+1,(site-1)*lc*3+(prn-1)*3+1)=sigma0/Aelev23(site).elev(ep,prn);
            flag_tm1=1;
        end
        % ---WL£¨1 -1 0£©
        if(FTamb12(site).amb(ep,prn)~=0)
            vt((site-1)*lc*3+(prn-1)*3+2,1)=FTamb12(site).amb(ep,prn);
            bt((site-1)*lc*3+(prn-1)*3+2,(site-1)*3+1)=1;
            bt((site-1)*lc*3+(prn-1)*3+2,(site-1)*3+2)=-1;
            bt((site-1)*lc*3+(prn-1)*3+2,(site-1)*3+3)=0;
            bt((site-1)*lc*3+(prn-1)*3+2,site_num_fcb*3+(prn-1)*3+1)=-1;
            bt((site-1)*lc*3+(prn-1)*3+2,site_num_fcb*3+(prn-1)*3+2)=1;
            bt((site-1)*lc*3+(prn-1)*3+2,site_num_fcb*3+(prn-1)*3+3)=0;
            pt((site-1)*lc*3+(prn-1)*3+2,(site-1)*lc*3+(prn-1)*3+2)=sigma0/Aelev12(site).elev(ep,prn);
            flag_tm2=2;
        end
        % ---NL£¨4 -3 0£©
        if(FTamb43(site).amb(ep,prn)~=0)
            vt((site-1)*lc*3+(prn-1)*3+3,1)=FTamb43(site).amb(ep,prn);
            bt((site-1)*lc*3+(prn-1)*3+3,(site-1)*3+1)=x1;
            bt((site-1)*lc*3+(prn-1)*3+3,(site-1)*3+2)=x2;
            bt((site-1)*lc*3+(prn-1)*3+3,(site-1)*3+3)=x3;
            bt((site-1)*lc*3+(prn-1)*3+3,site_num_fcb*3+(prn-1)*3+1)=-x1;
            bt((site-1)*lc*3+(prn-1)*3+3,site_num_fcb*3+(prn-1)*3+2)=-x2;
            bt((site-1)*lc*3+(prn-1)*3+3,site_num_fcb*3+(prn-1)*3+3)=-x3;
            pt((site-1)*lc*3+(prn-1)*3+3,(site-1)*lc*3+(prn-1)*3+3)=sigma0/Aelev43(site).elev(ep,prn);
            flag_tm3=3;
        end
    end
    if(flag_tm1*flag_tm2*flag_tm3==0 )
        rm_site_n=rm_site_n+1;
        rm_site(rm_site_n)=site;
    end
end
%% update FTamb23 FTamb12,FTamb43,Aelev23,Aelev12,Aelev43 based on the excluded stations
nsite=0;
for site=1:site_num_fcb
    % excluded stations
    site_tc_flag=0;
    for ns=1:rm_site_n
       if(site==rm_site(ns)); site_tc_flag=1; end
    end
    if(site_tc_flag==1); continue; end
    nsite=nsite+1;
    UDFTamb23(nsite)=FTamb23(site);
    UDFTamb12(nsite)=FTamb12(site);
    UDFTamb43(nsite)=FTamb43(site);
    UDAelev23(nsite)=Aelev23(site);
    UDAelev12(nsite)=Aelev12(site);
    UDAelev43(nsite)=Aelev43(site);
end
%% implement a second coefficient design based on the excluded stations
bt=zeros(nsite*3*lc+3,nsite*3+lc*3);
vt=zeros(nsite*3*lc+3,1);
pt=zeros(nsite*3*lc+3,nsite*3*lc+3);
k=0;
rec_flag=[];
for site=1:nsite
    for prn=1:lc  
        % ---UWL£¨0 1 -1£©
        if(UDFTamb23(site).amb(ep,prn)~=0)
            vt((site-1)*lc*3+(prn-1)*3+1,1)=UDFTamb23(site).amb(ep,prn);
            bt((site-1)*lc*3+(prn-1)*3+1,(site-1)*3+1)=0;
            bt((site-1)*lc*3+(prn-1)*3+1,(site-1)*3+2)=1;
            bt((site-1)*lc*3+(prn-1)*3+1,(site-1)*3+3)=-1;
            bt((site-1)*lc*3+(prn-1)*3+1,nsite*3+(prn-1)*3+1)=0;
            bt((site-1)*lc*3+(prn-1)*3+1,nsite*3+(prn-1)*3+2)=-1;
            bt((site-1)*lc*3+(prn-1)*3+1,nsite*3+(prn-1)*3+3)=1;
            pt((site-1)*lc*3+(prn-1)*3+1,(site-1)*lc*3+(prn-1)*3+1)=sigma0/UDAelev23(site).elev(ep,prn);
            k=k+1;
            rec_flag(k,1)=site;
            rec_flag(k,2)=prn;
            rec_flag(k,3)=1;
        end
        % ---WL£¨1 -1 0£©
        if(UDFTamb12(site).amb(ep,prn)~=0)
            vt((site-1)*lc*3+(prn-1)*3+2,1)=UDFTamb12(site).amb(ep,prn);
            bt((site-1)*lc*3+(prn-1)*3+2,(site-1)*3+1)=1;
            bt((site-1)*lc*3+(prn-1)*3+2,(site-1)*3+2)=-1;
            bt((site-1)*lc*3+(prn-1)*3+2,(site-1)*3+3)=0;
            bt((site-1)*lc*3+(prn-1)*3+2,nsite*3+(prn-1)*3+1)=-1;
            bt((site-1)*lc*3+(prn-1)*3+2,nsite*3+(prn-1)*3+2)=1;
            bt((site-1)*lc*3+(prn-1)*3+2,nsite*3+(prn-1)*3+3)=0;
            pt((site-1)*lc*3+(prn-1)*3+2,(site-1)*lc*3+(prn-1)*3+2)=sigma0/UDAelev12(site).elev(ep,prn);
            k=k+1;
            rec_flag(k,1)=site;
            rec_flag(k,2)=prn;
            rec_flag(k,3)=2;
        end
        % ---NL£¨4 -3 0£©
        if(UDFTamb43(site).amb(ep,prn)~=0)
            vt((site-1)*lc*3+(prn-1)*3+3,1)=UDFTamb43(site).amb(ep,prn);
            bt((site-1)*lc*3+(prn-1)*3+3,(site-1)*3+1)=x1;
            bt((site-1)*lc*3+(prn-1)*3+3,(site-1)*3+2)=x2;
            bt((site-1)*lc*3+(prn-1)*3+3,(site-1)*3+3)=x3;
            bt((site-1)*lc*3+(prn-1)*3+3,nsite*3+(prn-1)*3+1)=-x1;
            bt((site-1)*lc*3+(prn-1)*3+3,nsite*3+(prn-1)*3+2)=-x2;
            bt((site-1)*lc*3+(prn-1)*3+3,nsite*3+(prn-1)*3+3)=-x3;
            pt((site-1)*lc*3+(prn-1)*3+3,(site-1)*lc*3+(prn-1)*3+3)=sigma0/UDAelev43(site).elev(ep,prn);
            k=k+1;
            rec_flag(k,1)=site;
            rec_flag(k,2)=prn;
            rec_flag(k,3)=3;
        end
    end
end
if(cfg.ref_prn==0) % reference satellite
    vt(nsite*3*lc+1,1)=0;
    vt(nsite*3*lc+2,1)=0;
    vt(nsite*3*lc+3,1)=0;
    pt(nsite*3*lc+1,nsite*3*lc+1)=10000;
    pt(nsite*3*lc+2,nsite*3*lc+2)=10000;
    pt(nsite*3*lc+3,nsite*3*lc+3)=10000;
    bt(nsite*3*lc+1,nsite*3+(ref_prn-1)*3+1)=-1;
    bt(nsite*3*lc+2,nsite*3+(ref_prn-1)*3+2)=-1;
    bt(nsite*3*lc+3,nsite*3+(ref_prn-1)*3+3)=-1;
end
if(cfg.ref_prn==1) % reference receiver
    vt(nsite*3*lc+1,1)=0;
    vt(nsite*3*lc+2,1)=0;
    vt(nsite*3*lc+3,1)=0;
    pt(nsite*3*lc+1,nsite*3*lc+1)=10000;
    pt(nsite*3*lc+2,nsite*3*lc+2)=10000;
    pt(nsite*3*lc+3,nsite*3*lc+3)=10000;
    bt(nsite*3*lc+1,(ref_prn-1)*3+1)=1;
    bt(nsite*3*lc+2,(ref_prn-1)*3+2)=1;
    bt(nsite*3*lc+3,(ref_prn-1)*3+3)=1;
end

%% remove rows with all zeros
all_ln=nsite*3*lc+3;
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
all_lc=nsite*3+lc*3;
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
for c=1:k
    if(tc_lc(c)<=nsite*3)
        site=fix(tc_lc(c)/3)+1;
        fprintf('Warning: the station of %s have no data!',UDFTamb23(site).name);
%         stop
    else
        if(rem(tc_lc(c)-nsite*3,3)>=1); prn=fix((tc_lc(c)-nsite*3)/3)+1; end % PRN
        if(rem(tc_lc(c)-nsite*3,3)==0); prn=fix((tc_lc(c)-nsite*3)/3); end % PRN
        amb=tc_lc(c)-nsite*3-(prn-1)*3;  % L1 L2 L3
        tc=tc_lc(c)-tq;
        btt(:,tc)=[];
        prn_flag(prn,amb)=999.99;
        tq=tq+1;
    end
end
B=btt;
v=vtt;
P=ptt;
end
