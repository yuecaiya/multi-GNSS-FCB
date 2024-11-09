function [final_upd_uwn,QSupd_fix_amb]=UD_estimate_upd_uwn(yyyy,doy,Tamb23,Tamb12,Tamb43,FTamb23,FTamb12,FTamb43,Aelev23,Aelev12,Aelev43,site_num_fcb,ref_prn,cfg,ref_prn1,QSupd_fix_amb_w,uwn_flag)
%% UWL-WL-NL FCB/UPD estimation
% args:
%     FTamb23 : fractional part of UWL combination ambiguity 
%     FTamb12 : fractional part of WL combination ambiguity
%     FTamb43 : fractional part of NL combination ambiguity
%     Tamb23  : UWL combination ambiguity 
%     Tamb12  : WL combination ambiguity
%     Tamb43  : NL combination ambiguity
%     Aelev23 : weight for UWL
%     Aelev12 : weight for  WL
%     Aelev43 : weight for  NL
%     ref_prn : reference satellite or receiver
%     ref_prn1: reference satellite
%     site_num_fcb: site number
%     cfg         : configuration information
% return:
%     fina_upd : FCB/UPD values
%     QSupd_fix_amb: interface to be assigned
% made by Caiya Yue @ CUMTB and CASM
% ----
%% FCB/UPD estimation epoch by epoch
ln=length(FTamb23(1).amb(:,1));
lc=length(FTamb23(1).amb(1,:));
if(strcmp(uwn_flag,'u')); out_debug=strcat(cfg.out_dir,'/upd_debug_u_',yyyy,doy,'.txt'); end
if(strcmp(uwn_flag,'w')); out_debug=strcat(cfg.out_dir,'/upd_debug_w_',yyyy,doy,'.txt'); end
if(strcmp(uwn_flag,'n')); out_debug=strcat(cfg.out_dir,'/upd_debug_n_',yyyy,doy,'.txt'); end
if(strcmp(uwn_flag,'u')); fprintf('    **************** Solving UWL FCB ****************\n'); end
if(strcmp(uwn_flag,'w')); fprintf('    **************** Solving WL FCB *****************\n'); end
if(strcmp(uwn_flag,'n')); fprintf('    **************** Solving NL FCB *****************\n'); end
fid=fopen(out_debug,'w');
QSupd_fix_amb.amb=zeros(ln,1);
for ep=1:ln 
   if(ep==1); gps_upd=0; end
   % constructing coefficient matrix B=[LC23 LC12 LC43]'; v=[famb23 famb12 famb43]'; Q->P
   if(strcmp(uwn_flag,'u') || strcmp(uwn_flag,'w'))
       [v,B,v0,B0,P,prn_flag,rec_flag,rm_site_n,rm_site]=make_coe_matrix_uwn(FTamb23,FTamb12,FTamb43,Aelev23,Aelev12,Aelev43,site_num_fcb,ref_prn,cfg,ep,uwn_flag);
   end
   if(strcmp(uwn_flag,'n'))
       upd_wm=QSupd_fix_amb_w;
       [v,B,v0,B0,P,prn_flag,rec_flag,rm_site_n,rm_site]=make_coe_matrix_uwn1(FTamb23,FTamb12,FTamb43,Aelev23,Aelev12,Aelev43,site_num_fcb,ref_prn,cfg,ep,upd_wm,uwn_flag,gps_upd);
   end
    % UPD/FCB initial value determination
    v_f=v;
%     x0=pinv((B'*P*B))*B'*P*v;
    nsite=site_num_fcb-rm_site_n;
    if(cfg.ref_prn==0)
        [x0]=UD_estimate_UPD_x0(v0,B0,P,prn_flag,rec_flag,ref_prn,nsite,FTamb23,cfg);
    end
    lnx0=length(x0);
    lnxb=length(B(1,:));
    if(lnx0~=lnxb)
        fprintf('Warning:***UPD / FCB solution error at receiver and satellite end\n');
        x0=pinv((B'*P*B))*B'*P*v;
    end
    % iterative calculation of UPD/FCB
   [gps_upd,QSupd1]=UD_estimate_upd_sub_l_uwn(x0,FTamb23,FTamb12,FTamb43,Tamb23,Tamb12,Tamb43,rm_site_n,rm_site,v,B,P,v_f,prn_flag,rec_flag,site_num_fcb,cfg,ep,fid,uwn_flag,ref_prn);
    % store UPD/FCB
    final_upd_uwn(ep).upd(:,1)=gps_upd;
    QSupd_fix_amb(ep).amb=QSupd1;   
end
fclose(fid);
%% plot UPD/FCB
% UWL
if(strcmp(uwn_flag,'u'))    
    figure(1)
    for s=1:lc
        tmp_fcb=[];
        k=0;
        for e=1:ln
            if(final_upd_uwn(e).upd(s,1)<100 && final_upd_uwn(e).upd(ref_prn1,1)<100)
                k=k+1;
                tmp_fcb(k,1)=e;
                tmp_fcb(k,2)=final_upd_uwn(e).upd(s,1)-final_upd_uwn(e).upd(ref_prn1,1);
            end
        end
        if(k==0); continue; end
        plot(tmp_fcb(:,1),tmp_fcb(:,2),'.');
        hold on
    end 
    title('UWL FCB');
    ylim([-1.2 1.2]);
end
% WL
if(strcmp(uwn_flag,'w'))    
    figure(2)
    for s=1:lc
        tmp_fcb=[];
        k=0;
        for e=1:ln
            if(final_upd_uwn(e).upd(s,1)<100 && final_upd_uwn(e).upd(ref_prn1,1)<100)
                k=k+1;
                tmp_fcb(k,1)=e;
                tmp_fcb(k,2)=final_upd_uwn(e).upd(s,1)-final_upd_uwn(e).upd(ref_prn1,1);
            end
        end
        if(k==0); continue; end
        plot(tmp_fcb(:,1),tmp_fcb(:,2),'.');
        hold on
    end   
    title('WL FCB');
    ylim([-1.2 1.2]);
end
% NL
if(strcmp(uwn_flag,'n'))   
    figure(3)
    for s=1:lc
        tmp_fcb=[];
        k=0;
        for e=1:ln
            if(final_upd_uwn(e).upd(s,1)<100)
                k=k+1;
                tmp_fcb(k,1)=e;
                tmp_fcb(k,2)=final_upd_uwn(e).upd(s,1);
            end
        end
        if(k==0); continue; end
        plot(tmp_fcb(:,1),tmp_fcb(:,2),'.');
        hold on
       % fprintf('%d %5.3f\n',s,std(tmp_fcb(:,2)));
    end  
    title('NL FCB');
    ylim([-1.2 1.2]);
end
end
