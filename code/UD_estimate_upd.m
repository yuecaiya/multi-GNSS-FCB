function [fina_upd,QSupd]=UD_estimate_upd(yyyy,doy,FTamb23,FTamb12,FTamb43,Aelev23,Aelev12,Aelev43,site_num_fcb,ref_prn_rec,cfg,ref_prn)
%% Estimation on raw frequency FCB/UPD
% args:
%     FTamb23 : fractional part of UWL combination ambiguity
%     FTamb12 : fractional part of WL combination ambiguity
%     FTamb43 : fractional part of NL combination ambiguity
%     Aelev23 : weight for UWL
%     Aelev12 : weight for  WL
%     Aelev43 : weight for  NL
%     ref_prn : reference satellite or receiver
%     site_num_fcb: site number
%     cfg         : configuration information
% return:
%     fina_upd : FCB/UPD values
%     QSupd    : interface to be assigned
% made by Caiya Yue @ CUMTB and CASM
% ----
%% FCB/UPD estimation epoch by epoch
ln=length(FTamb23(1).amb(:,1));
lc=length(FTamb23(1).amb(1,:));
out_debug=strcat(cfg.out_dir,'/upd_debug',yyyy,doy,'.txt');
fid=fopen(out_debug,'w');
for ep=1:ln  
    % constructing coefficient matrix B=[LC23 LC12 LC43]'; v=[famb23 famb12 famb43]'; Q->P
   [v,B,P,prn_flag,rec_flag,rm_site_n,rm_site]=make_coe_matrix(FTamb23,FTamb12,FTamb43,Aelev23,Aelev12,Aelev43,site_num_fcb,ref_prn_rec,cfg,ep);
    % UPD/FCB initial value determination
    v_f=v;
    x0=pinv((B'*P*B))*B'*P*v;
    % iterative calculation of UPD/FCB
   [gps_upd,QSupd1]=UD_estimate_upd_sub_l(x0,FTamb23,FTamb12,FTamb43,rm_site_n,rm_site,v,B,P,v_f,prn_flag,rec_flag,site_num_fcb,cfg,ep,fid);
    % store UPD/FCB
    QSupd(ep,:)=QSupd1;
    fina_upd(ep).upd=gps_upd;
end
fclose(fid);
