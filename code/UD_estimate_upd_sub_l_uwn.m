function [gps_upd,QSupd1]=UD_estimate_upd_sub_l_uwn(x0,FTamb23,FTamb12,FTamb43,Tamb23,Tamb12,Tamb43,rm_site_n,rm_site,v,B,P,v_f,prn_flag,rec_flag,site_num_fcb,cfg,ep,fid,uwn_flag,ref_prn)
%% iterative calculation of UPD/FCB
% args:
%          x0 : initial value of UPD/FCB
%     FTamb23 : fractional part of UWL combination ambiguity 
%     FTamb12 : fractional part of WL combination ambiguity
%     FTamb43 : fractional part of NL combination ambiguity
%     Tamb23  : UWL combination ambiguity 
%     Tamb12  : WL combination ambiguity
%     Tamb43  : NL combination ambiguity
%     Aelev23 : weight for UWL
%     Aelev12 : weight for  WL
%     Aelev43 : weight for  NL
%     v       : residual vector 
%     B       : coefficient matrix
%     P       : weight matrix
%     ref_prn : reference satellite or receiver
%     site_num_fcb: site number
%     cfg         : configuration information
% return:
%     gps_upd : FCB/UPD values
%     QSupd    : interface to be assigned
% made by Caiya Yue @ CUMTB and CASM
% ----
%% interpreting control information
lan=length(B(1,:));
lnc=length(prn_flag(:,1));
fix_yz=str2num(cfg.fix_amb);
tmp_fix_n1=0; 
QSupd=zeros(1,lnc);
%% update FTamb23,FTamb12,FTamb43,Aelev23,Aelev12,Aelev43 based on the excluded stations
TMPamb=0;
if(strcmp(uwn_flag,'u')); FTMPamb=FTamb23; end
if(strcmp(uwn_flag,'w')); FTMPamb=FTamb12; end
if(strcmp(uwn_flag,'n')); FTMPamb=FTamb43; end
if(strcmp(uwn_flag,'u')); TMPamb=Tamb23; end
if(strcmp(uwn_flag,'w')); TMPamb=Tamb12; end
if(strcmp(uwn_flag,'n')); TMPamb=Tamb43; end
ys_nsite=site_num_fcb;
nsite=0;
for site=1:site_num_fcb
    % excluded stations
    site_tc_flag=0;
    for ns=1:rm_site_n
       if(site==rm_site(ns)); site_tc_flag=1; end
    end
    if(site_tc_flag==1); continue; end
    nsite=nsite+1;
    F_TMPamb(nsite)=FTMPamb(site);
    Y_TMPamb(nsite)=TMPamb(site);
end
site_num_fcb=nsite;
%% iterative estimation
% store the satellite FCB in an array ---------------------------------(1)
x=x0;
pk=0;
m=1;
gps_upd=prn_flag;
for i=site_num_fcb+1:length(x) 
    for j=m:lnc+1
        if(gps_upd(j)==999.99); continue; end
        gps_upd(j)=x(i);
        m=j+1;
        break;
    end
end
% initial update ---------------------------------(2)
for i=1:length(v_f)
    if(i>length(v_f)-1)
       pk=pk+1;
       v(pk,1)=0;
       continue;
    end
    site_i=rec_flag(i,1);
    prn_i=rec_flag(i,2);
%     if(uwn_i==1)
        namb=v_f(i)-(x(site_i)-gps_upd(prn_i));
        pk=pk+1;
        v(pk,1)=namb-round(namb);
        if(namb>fix_yz(1)); P(pk,pk)=fix_yz(3)*P(pk,pk); end
%     end
end
% continuous iterative estimation ---------------------------------(3)
iteration_i=0; % total number of iterations
fprintf('    Numbers of iterations for %2dst epoch: ',ep);
outf=strcat(cfg.out_dir,'/ResFcb_ep',num2str(ep));
fout=fopen(outf,'w');
if(strcmp(uwn_flag,'u')); fix_yz_uwn=fix_yz(2); end
if(strcmp(uwn_flag,'w')); fix_yz_uwn=fix_yz(2); end
if(strcmp(uwn_flag,'n')); fix_yz_uwn=fix_yz(2); end
ka2=0; 
while(1)
    iteration_i=iteration_i+1;   
    fprintf(fout,'iteration %3d\n',iteration_i);
    % least squares solution
    dx=pinv((B'*P*B))*B'*P*v; 
    x=x+dx;
    % interpret satellite UPD/FCB
    m=1;
    gps_upd=prn_flag;
    for i=site_num_fcb+1:length(x)    
        for j=m:lnc+1
            if(gps_upd(j)==999.99); continue; end
            gps_upd(j)=x(i);
            m=j+1;
            break;
        end
    end
    % fixed ambiguity
    ka1=0; % total number of station-satellite
    k1=0;  % fixed number of station-satellite (<0.25c)
    k11=0; % fixed number of station-satellite (<0.15c)    
    pk=0;  
    for i=1:length(v_f)
       if(i>length(v_f)-1)
           pk=pk+1;
%            v(pk,1)=0;
           v(pk,1)=gps_upd(ref_prn);
           continue;
       end
        site_i=rec_flag(i,1);
        prn_i=rec_flag(i,2);
        prn=prn_i;
        res_uwn=zeros(lnc,1);
%         if(uwn_i==1)
            namb=v_f(i)-(x(site_i)-gps_upd(prn_i));
            ka1=ka1+1;
            pk=pk+1;
            v(pk,1)=namb-round(namb);
            res_uwn(prn,1)=namb-round(namb);
            namb=abs(round(namb)-namb);
            if(namb<fix_yz_uwn) %fix_yz(2)
                k1=k1+1;
                P(pk,pk)=1*P(pk,pk);
                if(namb<0.15); k11=k11+1; end
            elseif(namb+1<fix_yz_uwn) %fix_yz(2)
                k1=k1+1;
                P(pk,pk)=1*P(pk,pk);
                v_f(i)=v_f(i)+1;
                if(namb+1<0.15); k11=k11+1; end
            elseif(abs(namb-1)<fix_yz_uwn) %fix_yz(2)
                k1=k1+1;
                P(pk,pk)=1*P(pk,pk);
                v_f(i)=v_f(i)-1;
                if(namb+1<0.15); k11=k11+1; end
            elseif(namb>fix_yz_uwn) %fix_yz(2)
                 P(pk,pk)=fix_yz(3)*P(pk,pk);
            end 
%         end        
        % residual output for each iteration
        if(strcmp(uwn_flag,'u'))
            fprintf(fout,'UW_%s ',F_TMPamb(nsite).name(1:4)); %UWL
            for prn=1:lnc
                fprintf(fout,'%6.3f ',res_uwn(prn,1));
            end
            fprintf(fout,'\n');
        end
        %
        if(strcmp(uwn_flag,'w'))
            fprintf(fout,'WL_%s ',F_TMPamb(nsite).name(1:4)); %WL
            for prn=1:lnc
                fprintf(fout,'%6.3f ',res_uwn(prn,1));
            end
            fprintf(fout,'\n');
        end
        %
        if(strcmp(uwn_flag,'n'))
            fprintf(fout,'NL_%s ',F_TMPamb(nsite).name(1:4)); %NL
            for prn=1:lnc
                fprintf(fout,'%6.3f ',res_uwn(prn,1));
            end
            fprintf(fout,'\n');
        end
    end
    %% fixed ambiguity rate  
    if(k1==tmp_fix_n1 || iteration_i>25)
        if(strcmp(uwn_flag,'u')); fprintf(fid,'Numbers of iterations and fix rate of UWL FCB %3d %6.3f %6.3f\n',iteration_i,k1/ka1,k11/ka1); end
        if(strcmp(uwn_flag,'w')); fprintf(fid,'Numbers of iterations and fix rate of WL FCB %3d %6.3f %6.3f\n',iteration_i,k1/ka1,k11/ka1); end
        if(strcmp(uwn_flag,'n')); fprintf(fid,'Numbers of iterations and fix rate of NL FCB %3d %6.3f %6.3f\n',iteration_i,k1/ka1,k11/ka1); end
        break;
    else
        tmp_fix_n1=k1;
        fprintf('%3d ',iteration_i);
    end  
end
fclose(fout);
%% assign initial value to the next epoch
x0=x;
fprintf('\n');
%% fixed ambiguity extraction
QSupd1=0;
% interpretation of satellite FCB
m=1;
gps_upd=prn_flag;
for i=site_num_fcb+1:length(x)    
    for j=m:lnc+1
        if(gps_upd(j)==999.99); continue; end
        gps_upd(j)=x(i);
        m=j+1;
        break;
    end
end
% ambiguity extraction
k1=0;
ka=0;
for s=1:lnc
    k=0;
    tmpupd=[];
    for i=1:length(v_f)-1
        site_i=rec_flag(i,1);
        prn_i=rec_flag(i,2);
        if(prn_i==s)
            k=k+1;
            ka=ka+1;
            tmpupd(k,1)=Y_TMPamb(site_i).amb(ep,prn_i)-(x(site_i)-gps_upd(prn_i));
            if(abs(tmpupd(k,1)-round(tmpupd(k,1)))<0.25)
                Y_TMPamb(site_i).amb(ep,prn_i)=round(Y_TMPamb(site_i).amb(ep,prn_i)-(x(site_i)-gps_upd(prn_i)));
                k1=k1+1;
            end            
        end
    end
%     fprintf('%5.3f ',k1/ka);
end
QSupd1=Y_TMPamb;
%  fprintf('\n');
end