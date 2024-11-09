function [stat]=out_upd_file_nl(yyyy,doy,final_upd,cfg,QSupd,ref_prn)
%% output and plot FCB/UPD 
% args:
%     final_upd : FCB/UPD values
%     QSupd     :  interface to be assigned
%     cfg       : configuration information
% return:
%     
% made by Caiya Yue @ CUMTB and CASM
% ----
%% interpreting control information
int_e=cfg.span_t;
prn=length(final_upd(1).upd(:,1));
nep=length(final_upd);
xs=str2num(cfg.com3);
x1=xs(1);
x2=xs(2);
x3=xs(3);
if(cfg.NavSystem==1); sys='G'; end
if(cfg.NavSystem==2); sys='E'; end
if(cfg.NavSystem==3); sys='C'; end
if(cfg.NavSystem==4); sys='C'; end
%% Output combined UPD/FCB
outf=strcat(cfg.out_dir,'\cgs',yyyy,doy,'_LC_',sys,'.fcb');
fid=fopen(outf,'w');
[cmonth,cday]=doy2mon_day(yyyy,doy); % month,day
can_t=strcat('* ',yyyy,'_',cmonth,'_',cday);
[final_upd_a]=UD_process_upd_LC(final_upd,cfg);
% output FCB file header
st=cfg.span_t; 
if(cfg.NavSystem==1); fprintf(fid,'      1.00             FCB DATA                   G         VERSION / TYPE\n'); end
if(cfg.NavSystem==2); fprintf(fid,'      1.00             FCB DATA                   E         VERSION / TYPE\n'); end
if(cfg.NavSystem==3); fprintf(fid,'      1.00             FCB DATA                   C2        VERSION / TYPE\n'); end
if(cfg.NavSystem==3); fprintf(fid,'      1.00             FCB DATA                   C3        VERSION / TYPE\n'); end
fprintf(fid,'      CGS                                                   RUN BY / DATE\n');
fprintf(fid,'CASM,                                                       ANALYSIS CENTER\n');
fprintf(fid,'Email:                                                      COMMENT\n');
fprintf(fid,'GEC2C3    igs      igs14.atx                                SYS/EXT PROD APPLIED\n');
fprintf(fid,'    1                                                       # OF SOLN STA\n');
fprintf(fid,'****                                                        STA NAME LIST\n');
fprintf(fid,'****     UWL WL NL (cycle)                                 FCB type\n');
fprintf(fid,' %4.1f ',st);
fprintf(fid,'min                                                   SAMPLING INTERVAL\n');
% output file body
for i=1:nep
    time=int_e*(i-1);
    h=floor(time/60);
    m=time-h*60;
    s=0;
    ln=length(can_t);
    fprintf(fid,'%s ',can_t(1:ln));
    fprintf(fid,'%2d %2d %2d\n',h,m,s); 
    for j=1:prn
        % UWL WL NL UPD
        if(j<=9); fprintf(fid,'%s0%1d ',sys,j); end
        if(j>9); fprintf(fid,'%s%2d ',sys,j); end
        if(final_upd_a(i).upd(j,1)==0)
            final_upd_a(i).upd(j,1)=99999.9990; 
        else
        end
        if(final_upd_a(i).upd(j,2)==0)
            final_upd_a(i).upd(j,2)=99999.9990; 
        end
        if(final_upd_a(i).upd(j,3)==0)
            final_upd_a(i).upd(j,3)=99999.9990;
        end
        fprintf(fid,'%15.4f %15.4f %15.4f\n',final_upd_a(i).upd(j,1),final_upd_a(i).upd(j,2),final_upd_a(i).upd(j,3));
    end 
end
%% output raw frequency FCB
outf=strcat(cfg.out_dir,'\cgs',yyyy,doy,'_raw_',sys,'.fcb');
fid=fopen(outf,'w');
[cmonth,cday]=doy2mon_day(yyyy,doy); % month,day
can_t=strcat('* ',yyyy,'_',cmonth,'_',cday);
[final_upd_a]=UD_process_upd_raw(final_upd,ref_prn,cfg);
% output FCB file header
st=cfg.span_t; 
if(cfg.NavSystem==1); fprintf(fid,'      1.00             FCB DATA                   G         VERSION / TYPE\n'); end
if(cfg.NavSystem==2); fprintf(fid,'      1.00             FCB DATA                   E         VERSION / TYPE\n'); end
if(cfg.NavSystem==3); fprintf(fid,'      1.00             FCB DATA                   C2        VERSION / TYPE\n'); end
if(cfg.NavSystem==3); fprintf(fid,'      1.00             FCB DATA                   C3        VERSION / TYPE\n'); end
fprintf(fid,'      CGS                                                   RUN BY / DATE\n');
fprintf(fid,'CASM,                                                       ANALYSIS CENTER\n');
fprintf(fid,' Email:                                                     COMMENT\n');
fprintf(fid,'GEC2C3    igs      igs14.atx                                SYS/EXT PROD APPLIED\n');
fprintf(fid,'    1                                                       # OF SOLN STA\n');
fprintf(fid,'****                                                        STA NAME LIST\n');
fprintf(fid,'****     L1 L2 L3 (cycle)                                   FCB type\n');
fprintf(fid,' %4.1f ',st);
fprintf(fid,'min                                                   SAMPLING INTERVAL\n');
% output file body
for i=1:nep
    time=int_e*(i-1);
    h=floor(time/60);
    m=time-h*60;
    s=0;
    ln=length(can_t);
    fprintf(fid,'%s ',can_t(1:ln));
    fprintf(fid,'%2d %2d %2d\n',h,m,s); 
    for j=1:prn
        % raw frequency UPD
        if(j<=9); fprintf(fid,'%s0%1d ',sys,j); end
        if(j>9); fprintf(fid,'%s%2d ',sys,j); end
        if(final_upd_a(i).upd(j,1)==0)
            final_upd_a(i).upd(j,1)=99999.9990; 
        end
        if(final_upd_a(i).upd(j,2)==0)
            final_upd_a(i).upd(j,2)=99999.9990; 
        end
        if(final_upd_a(i).upd(j,3)==0)
            final_upd_a(i).upd(j,3)=99999.9990;
        end
        fprintf(fid,'%15.4f %15.4f %15.4f\n',final_upd_a(i).upd(j,1),final_upd_a(i).upd(j,2),final_upd_a(i).upd(j,3));
    end
end
stat=0;