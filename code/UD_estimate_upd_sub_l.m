function [gps_upd,QSupd]=UD_estimate_upd_sub_l(xi,FTamb23,FTamb12,FTamb43,rm_site_n,rm_site,v,B,P,v_f,prn_flag,rec_flag,site_num_fcb,cfg,ep,fid)
%% iterative calculation of UPD/FCB
% args:
%     FTamb23 : UWL combination ambiguity
%     FTamb12 :  WL combination ambiguity
%     FTamb43 :  NL combination ambiguity
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
xs=str2num(cfg.com3);
x1=xs(1);
x2=xs(2);
x3=xs(3);

lan=length(B(1,:));
lnc=length(prn_flag(:,1));
fix_yz=str2num(cfg.fix_amb);

tmp_fix_n1=0; 
tmp_fix_n2=0;
tmp_fix_n3=0; 

nuse_sat=0;
suse_sat=zeros(lnc,2);
QSupd=zeros(1,lnc);
%% update FTamb23,FTamb12,FTamb43,Aelev23,Aelev12,Aelev43 based on the excluded stations
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
    UDFTamb23(nsite)=FTamb23(site);
    UDFTamb12(nsite)=FTamb12(site);
    UDFTamb43(nsite)=FTamb43(site);
end
site_num_fcb=nsite;
%% iterative estimation
% store the satellite FCB in an array ---------------------------------(1)
x=xi;
gps_upd=prn_flag;
m=1;
for i=site_num_fcb*3+1:length(x)    
    for j=m:lnc*3+1
        ii=fix(j/3)+1; % row
        if(rem(j,3)==0); ii=fix(j/3); end
        ij=rem(j,3); % column
        if(ij==0); ij=3; end
        if(gps_upd(ii,ij)==999.99); continue; end
        gps_upd(ii,ij)=x(i);
        m=j+1;
        break;
    end
end
% initial update ---------------------------------(2)
pk=0;
for i=1:length(v_f)
    if(i>length(v_f)-3)
       pk=pk+1;
       v(pk,1)=0;
       continue;
    end
    site_i=rec_flag(i,1);
    prn_i=rec_flag(i,2);
    uwn_i=rec_flag(i,3);
    % UWL
    if(uwn_i==1)
        namb=v_f(i)-( x((site_i-1)*3+2)-gps_upd(prn_i,2)-(x((site_i-1)*3+3)-gps_upd(prn_i,3)) );
        pk=pk+1;
        v(pk,1)=namb-round(namb);
        if(namb>fix_yz(1)); P(pk,pk)=fix_yz(3)*P(pk,pk); end
    end
    % WL
    if(uwn_i==2)
        namb=v_f(i)-( x((site_i-1)*3+1)-gps_upd(prn_i,1)-(x((site_i-1)*3+2)-gps_upd(prn_i,2)) );
        pk=pk+1;
        v(pk,1)=namb-round(namb);
        if(namb>fix_yz(1)); P(pk,pk)=fix_yz(3)*P(pk,pk); end
    end
    % NL
    if(uwn_i==3)
        namb=v_f(i)-( x1*(x((site_i-1)*3+1)-gps_upd(prn_i,1))+x2*(x((site_i-1)*3+2)-gps_upd(prn_i,2))+x3*(x((site_i-1)*3+3)-gps_upd(prn_i,3)) );
        pk=pk+1;
        v(pk,1)=namb-round(namb);
        if(namb>fix_yz(1)); P(pk,pk)=fix_yz(3)*P(pk,pk); end
    end
end
% continuous iterative estimation ---------------------------------(3)
iteration_i=0; % total number of iterations
fprintf('    Numbers of iterations for %2dst epoch: ',ep);
outf=strcat(cfg.out_dir,'/ResFcb_ep',num2str(ep));
fout=fopen(outf,'w');
while(1)
    iteration_i=iteration_i+1;   
    fprintf(fout,'iteration %3d\n',iteration_i);
    % least squares solution
    dx=pinv((B'*P*B))*B'*P*v;
    x=x+dx;
    % interpret satellite UPD/FCB
    gps_upd=prn_flag;
    m=1;
    for i=site_num_fcb*3+1:length(x)    
        for j=m:lnc*3+1
            ii=fix(j/3)+1; % row
            if(rem(j,3)==0); ii=fix(j/3); end
            ij=rem(j,3); % column
            if(ij==0); ij=3; end
            if(gps_upd(ii,ij)==999.99); continue; end
            gps_upd(ii,ij)=x(i);
            m=j+1;
            break;
        end
    end
    % fixed ambiguity
    ka1=0;  % total number of station-satellite£¨0 1 -1£©
    k1=0;   % fixed number of station-satellite£¨0 1 -1£©->0.25c
    k11=0;  % fixed number of station-satellite£¨0 1 -1£©->0.15c
    ka2=0;  % total number of station-satellite£¨1 -1 0£©
    k2=0;   % fixed number of station-satellite£¨1 -1 0£©->0.25c
    k22=0;  % fixed number of station-satellite£¨0 1 -1£©->0.15c
    ka3=0;  % total number of station-satellite£¨4 -3 0£©
    k3=0;   % fixed number of station-satellite£¨4 -3 0£©->0.25c
    k33=0;  % fixed number of station-satellite£¨0 1 -1£©->0.15c
    pk=0;   
    for i=1:length(v_f)
       if(i>length(v_f)-3)
           pk=pk+1;
           v(pk,1)=0;
           continue;
       end
        site_i=rec_flag(i,1);
        prn_i=rec_flag(i,2);
        prn=prn_i;
        uwn_i=rec_flag(i,3);
        res_u=zeros(lnc,1);
        res_w=zeros(lnc,1);
        res_n=zeros(lnc,1);
        % UWL
        if(uwn_i==1)
            namb=v_f(i)-( x((site_i-1)*3+2)-gps_upd(prn_i,2)-(x((site_i-1)*3+3)-gps_upd(prn_i,3)) );
            ka1=ka1+1;
            pk=pk+1;
            v(pk,1)=namb-round(namb);
            res_u(prn,1)=namb-round(namb);
            namb=abs(round(namb)-namb);
            if(namb<fix_yz(2)) % fix_yz(2)
                k1=k1+1;
                P(pk,pk)=1*P(pk,pk);
                if(namb<0.15); k11=k11+1; end
            elseif(namb+1<fix_yz(2)) % fix_yz(2)
                k1=k1+1;
                P(pk,pk)=1*P(pk,pk);
                v_f(i)=v_f(i)+1;
                if(namb+1<0.15); k11=k11+1; end
            elseif(abs(namb-1)<fix_yz(2)) % fix_yz(2)
                k1=k1+1;
                P(pk,pk)=1*P(pk,pk);
                v_f(i)=v_f(i)-1;
                if(namb+1<0.15); k11=k11+1; end
            elseif(namb>fix_yz(1)) % fix_yz(1)
                 P(pk,pk)=fix_yz(3)*P(pk,pk);
            end 
        end
        % WL
        if(uwn_i==2)
            namb=v_f(i)-( x((site_i-1)*3+1)-gps_upd(prn_i,1)-(x((site_i-1)*3+2)-gps_upd(prn_i,2)) );
            ka2=ka2+1;
            pk=pk+1;
            v(pk,1)=namb-round(namb);
            res_w(prn,1)=namb-round(namb);
            namb=abs(round(namb)-namb);
            if(namb<fix_yz(2)) % fix_yz(2)
                k2=k2+1;
                P(pk,pk)=1*P(pk,pk);
                if(namb<0.15); k22=k22+1; end
            elseif(namb+1<fix_yz(2)) % fix_yz(2)
                k1=k1+1;
                P(pk,pk)=1*P(pk,pk);
                v_f(i)=v_f(i)+1;
                if(namb+1<0.15); k11=k11+1; end
            elseif(abs(namb-1)<fix_yz(2)) % fix_yz(2)
                k1=k1+1;
                P(pk,pk)=1*P(pk,pk);
                v_f(i)=v_f(i)-1;
                if(namb+1<0.15); k11=k11+1; end
            elseif(namb>fix_yz(1)) % fix_yz(1)
                P(pk,pk)=fix_yz(3)*P(pk,pk);
            end
        end
        % NL
        fix_yz_nl=0.35;
        if(uwn_i==3)
            namb=v_f(i)-( x1*(x((site_i-1)*3+1)-gps_upd(prn_i,1))+x2*(x((site_i-1)*3+2)-gps_upd(prn_i,2))+x3*(x((site_i-1)*3+3)-gps_upd(prn_i,3)) );
            ka3=ka3+1;
            pk=pk+1;
            v(pk,1)=namb-round(namb);
            res_n(prn,1)=namb-round(namb);
            namb=abs(round(namb)-namb);
            if(namb<fix_yz_nl)
                k3=k3+1;
                P(pk,pk)=1*P(pk,pk);
                if(namb<0.25); k33=k33+1; end
            elseif(namb+1<fix_yz_nl) % fix_yz(2)
                k1=k1+1;
                P(pk,pk)=1*P(pk,pk);
                UDFTamb43(nsite).amb(ep,prn)=UDFTamb43(nsite).amb(ep,prn)+1;
                if(namb+1<0.25); k11=k11+1; end
            elseif(abs(namb-1)<fix_yz_nl) % fix_yz(2)
                k1=k1+1;
                P(pk,pk)=1*P(pk,pk);
                UDFTamb43(nsite).amb(ep,prn)=UDFTamb43(nsite).amb(ep,prn)-1;
                if(namb+1<0.25); k11=k11+1; end
            elseif(namb>fix_yz_nl) % fix_yz(1)
                P(pk,pk)=fix_yz(3)*P(pk,pk);
            end
        end
        % residual output for each iteration
        fprintf(fout,'UW_%s ',UDFTamb23(nsite).name(1:4)); % UWL
        for prn=1:lnc
            fprintf(fout,'%6.3f ',res_u(prn,1));
        end
        fprintf(fout,'\n');
        %
        fprintf(fout,'WL_%s ',UDFTamb23(nsite).name(1:4)); % WL
        for prn=1:lnc
            fprintf(fout,'%6.3f ',res_w(prn,1));
        end
        fprintf(fout,'\n');
        %
        fprintf(fout,'NL_%s ',UDFTamb23(nsite).name(1:4)); % NL
        for prn=1:lnc
            fprintf(fout,'%6.3f ',res_n(prn,1));
        end
        fprintf(fout,'\n');
    end
    %% fixed ambiguity rate
    if((k1==tmp_fix_n1 && k2==tmp_fix_n2 && k3==tmp_fix_n3) || iteration_i>25)
        fprintf(fid,'Numbers of iterations and fix rate of three liner combination %3d %6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f\n',iteration_i,k1/ka1,k2/ka2,k3/ka3,k11/ka1,k22/ka2,k33/ka3);
        for s=1:lnc
            if(prn_flag(s,1)<100)
                nuse_sat=nuse_sat+1; 
                suse_sat(nuse_sat)=prn_flag(s,1);
            end
        end
        break;
    else
        tmp_fix_n1=k1;
        tmp_fix_n2=k2;
        tmp_fix_n3=k3;
        fprintf('%3d ',iteration_i);
    end  
end
fclose(fout);
%% assign initial value to the next epoch
x0=x;
fprintf('\n');