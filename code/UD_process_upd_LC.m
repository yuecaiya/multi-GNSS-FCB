function [final_upd_a]=UD_process_upd_LC(final_upd,cfg)
%% plot FCB/UPD 
% args:
%     final_upd : FCB/UPD values
%     cfg       : configuration information
% return:
%     final_upd_a: standard output format
% made by Caiya Yue @ CUMTB and CASM
%% interpreting control information
prn=length(final_upd(1).upd(:,1));
nep=length(final_upd);
xs=str2num(cfg.com3);
x1=xs(1);
x2=xs(2);
x3=xs(3);
%% extract the original frequency UPD
for s=1:prn
   upd_t=zeros(nep,3);
   for e=1:nep
       upd_t(e,1)=final_upd(e).upd(s,1);
       upd_t(e,2)=final_upd(e).upd(s,2);
       upd_t(e,3)=final_upd(e).upd(s,3);
   end
   tmp_upd(s).upd=upd_t;
end
% form a linear combination
for s=1:prn
   for e=1:nep
       % UWL
       if(tmp_upd(s).upd(e,2)>800 || tmp_upd(s).upd(e,3)>800)
           tmp_upd_uwn(s).upd(e,1)=999.99;
       else
           tmp_upd_uwn(s).upd(e,1)=tmp_upd(s).upd(e,2)-tmp_upd(s).upd(e,3);
       end
       % WL
       if(tmp_upd(s).upd(e,1)>800 || tmp_upd(s).upd(e,2)>800)
           tmp_upd_uwn(s).upd(e,2)=999.99;
       else
           tmp_upd_uwn(s).upd(e,2)=tmp_upd(s).upd(e,1)-tmp_upd(s).upd(e,2);
       end
       % NL
       if(tmp_upd(s).upd(e,1)>800 || tmp_upd(s).upd(e,2)>800)
           tmp_upd_uwn(s).upd(e,3)=999.99;
       else
           tmp_upd_uwn(s).upd(e,3)=x1*tmp_upd(s).upd(e,1)+x2*tmp_upd(s).upd(e,2)+x3*tmp_upd(s).upd(e,3);
       end
   end
end
%% include data between [-1 1] by adding or subtracting one week
for s=1:prn
    %% UWL
    upd_u=tmp_upd_uwn(s).upd(:,1);
    for n=1:nep
       if(upd_u(n)>800); continue; end
       for i=1:100
          if(abs(upd_u(n))<1)
              break;
          else
              if(upd_u(n)<-1); upd_u(n)=upd_u(n)+1; end
              if(upd_u(n)> 1); upd_u(n)=upd_u(n)-1; end
          end
       end
    end
    % based on the median
    med_u=median(upd_u);
    for n=1:nep
        if(upd_u(n)>800); continue; end
        if((upd_u(n)-med_u)>0.5); upd_u(n)=upd_u(n)-1; end
        if((upd_u(n)-med_u)<-0.5); upd_u(n)=upd_u(n)+1; end
    end
    %% WL
    upd_w=tmp_upd_uwn(s).upd(:,2);
    for n=1:nep
       if(upd_w(n)>800); continue; end
       for i=1:100
          if(abs(upd_w(n))<1)
              break;
          else
              if(upd_w(n)<-1); upd_w(n)=upd_w(n)+1; end
              if(upd_w(n)> 1); upd_w(n)=upd_w(n)-1; end
          end
       end
    end
    % based on the median
    med_u=median(upd_w);
    for n=1:nep
        if(upd_w(n)>800); continue; end
        if((upd_w(n)-med_u)>0.5); upd_w(n)=upd_w(n)-1; end
        if((upd_w(n)-med_u)<-0.5); upd_w(n)=upd_w(n)+1; end
    end
    %% NL
    upd_n=tmp_upd_uwn(s).upd(:,3);
    for n=1:nep
       if(upd_n(n)>800); continue; end
       for i=1:100
          if(abs(upd_n(n))<1)
              break;
          else
              if(upd_n(n)<-1); upd_n(n)=upd_n(n)+1; end
              if(upd_n(n)> 1); upd_n(n)=upd_n(n)-1; end
          end
       end
    end
    % based on the median
    med_u=median(upd_n);
    for n=1:nep
        if(upd_u(n)>800); continue; end
        if((upd_n(n)-med_u)>0.5); upd_n(n)=upd_n(n)-1; end
        if((upd_n(n)-med_u)<-0.5); upd_n(n)=upd_n(n)+1; end
    end
    tmp_upd_uwn_0(s).upd=[upd_u upd_w upd_n];
end
%% remove and interpolate outliers
fprintf('UWL WL NL LC FCB STD\n');
for s=1:prn
    % UWL
    upd_uwn=tmp_upd_uwn_0(s).upd(:,1);
    med_uwn=median(upd_uwn);
    upd_uwn_std=[];
    kstd=0;
    std1=999.99;
    for ii=1:nep
       if(upd_uwn(ii)<800) 
           kstd=kstd+1;
           upd_uwn_std(kstd)=upd_uwn(ii);
       end
    end
    if(kstd<nep/3)
        for e=1:nep
            tmp_upd_uwn_0(s).upd(e,1)=999.99;
        end
    else    
        std_uwn=std(upd_uwn_std);    
        for e=1:nep
            if(abs(upd_uwn(e)-med_uwn)>2*std_uwn)
                tmp_upd_uwn_0(s).upd(e,1)=med_uwn;
            end
        end
        tmp_upd_uwn_0(s).upd(:,1) = smooth(tmp_upd_uwn_0(s).upd(:,1),8,'rlowess'); 
        figure(1)
        tmpxy=[];
        kxy=0;
        for e=1:nep
            if(tmp_upd_uwn_0(s).upd(e,1)~=0)
                kxy=kxy+1;
                tmpxy(kxy,1)=e;
                tmpxy(kxy,2)=tmp_upd_uwn_0(s).upd(e,1);
            end
        end
        if(kxy==0); continue; end
        gcf1=plot(tmpxy(:,1),tmpxy(:,2),'o');
        hold on
        ylim([-1.2 1.2]);
        title('UWL FCB');
        std1=std(tmpxy(:,2));
        sv_fig1=strcat(cfg.out_dir,'/figure1');  
        saveas(gcf,sv_fig1,'png');
    end
    % WL
    upd_uwn=tmp_upd_uwn_0(s).upd(:,2);
    med_uwn=median(upd_uwn);
    upd_uwn_std=[];
    kstd=0;
    std2=999.99;
    for ii=1:nep
       if(upd_uwn(ii)<800) 
           kstd=kstd+1;
           upd_uwn_std(kstd)=upd_uwn(ii);
       end
    end
    if(kstd<nep/3)
        for e=1:nep
            tmp_upd_uwn_0(s).upd(e,2)=999.99;
        end
    else    
        std_uwn=std(upd_uwn_std);    
        for e=1:nep
            if(abs(upd_uwn(e)-med_uwn)>2*std_uwn)
                tmp_upd_uwn_0(s).upd(e,2)=med_uwn;
            end
        end
        tmp_upd_uwn_0(s).upd(:,2) = smooth(tmp_upd_uwn_0(s).upd(:,2),8,'rlowess'); % Gaussian window method
        figure(2)
        tmpxy=[];
        kxy=0;
        for e=1:nep
            if(tmp_upd_uwn_0(s).upd(e,2)~=0)
                kxy=kxy+1;
                tmpxy(kxy,1)=e;
                tmpxy(kxy,2)=tmp_upd_uwn_0(s).upd(e,2);
            end
        end
        if(kxy==0); continue; end        
        std2=std(tmpxy(:,2));
        if(std2<0.11)
            plot(tmpxy(:,1),tmpxy(:,2),'o');
            hold on
            ylim([-1.2 1.2]);
            title('WL FCB');              
            sv_fig2=strcat(cfg.out_dir,'/figure2');
            saveas(gcf,sv_fig2,'png');
        end
    end
    % NL
    upd_uwn=tmp_upd_uwn_0(s).upd(:,3);
    med_uwn=median(upd_uwn);
    upd_uwn_std=[];
    kstd=0;
    std3=999.99;
    for ii=1:nep
       if(upd_uwn(ii)<800) 
           kstd=kstd+1;
           upd_uwn_std(kstd)=upd_uwn(ii);
       end
    end
    if(kstd<nep/3)
        for e=1:nep
            tmp_upd_uwn_0(s).upd(e,3)=999.99;
        end
    else    
        std_uwn=std(upd_uwn_std);    
        for e=1:nep
            if(abs(upd_uwn(e)-med_uwn)>3*std_uwn)
                tmp_upd_uwn_0(s).upd(e,3)=med_uwn;
            end
        end
        tmp_upd_uwn_0(s).upd(:,3) = smooth(tmp_upd_uwn_0(s).upd(:,3),8,'rlowess'); % 高斯窗方法\
        figure(3)
        tmpxy=[];
        kxy=0;
        for e=1:nep
            if(tmp_upd_uwn_0(s).upd(e,3)~=0)
                kxy=kxy+1;
                tmpxy(kxy,1)=e;
                tmpxy(kxy,2)=tmp_upd_uwn_0(s).upd(e,3);
            end
        end
        if(kxy==0); continue; end
        std3=std(tmpxy(:,2));
        if(std3<0.21)
            plot(tmpxy(:,1),tmpxy(:,2),'o');
            hold on
            ylim([-1.2 1.2]);
            title('NL FCB');
            sv_fig3=strcat(cfg.out_dir,'/figure3');
            saveas(gcf,sv_fig3,'png');  
        end

    end
    if(s<10); sprn=strcat('0',num2str(s)); end
    if(s>9); sprn=num2str(s); end
    fprintf('PRN%s %10.4f %10.4f %10.4f\n',sprn,std1,std2,std3);
end
%% convert to standard output format
for e=1:nep
    for s=1:prn
        final_upd_a(e).upd(s,1)=tmp_upd_uwn_0(s).upd(e,1);
        final_upd_a(e).upd(s,2)=tmp_upd_uwn_0(s).upd(e,2);
        final_upd_a(e).upd(s,3)=tmp_upd_uwn_0(s).upd(e,3);
    end
end
end