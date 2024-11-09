function [amb,P]=UD_combine_amb_sub_arc(site,prn,k,data,cfg,flag_uwn,DTamb,PAelev)
%% merge ambiguity for continuous data without cycle slip
% made by Caiya Yue @ CUMTB and CASM
% ----
% 
jg=cfg.comb_arc*(60/cfg.inter);
int_m=cfg.span_t*2;
%% data mergence and weight calculation
result_data=0;
com_p=0;
dat=0;
% equal weight
if(cfg.weight==0);
    dat=mean(data(k-jg:k,2)); % ambiguity
    com_p=1;                  % weight
end
% weight calculation based on STD
if(cfg.weight==1);
    dat=mean(data(k-jg:k,2)); % ambiguity
    com_p=std(data(k-jg:k,2));
    if(com_p==0); com_p=1e-05; end
end
% weighted mean 
if(cfg.weight==2);
    sum_p=sum(1./data(k-jg:k,3));
    for i=k-jg:k
        dat=dat+((1/data(i,3))/sum_p)*data(i,2);
    end
    for i=k-jg:k
        com_p=com_p+((1/data(i,3))/sum_p)^2*data(i,3);
    end
end
result_data=dat; 
%% remove arc segments with STD greater than 0.2 cycle, from 200 epochs to the end
if(strcmp(flag_uwn,'n'))
    if(std(data(k-fix(k/3):k,2))>0.2); result_data=0; end
end
%% assignment within the time period
if(k>int_m)
  for ep=1:k
    nep=fix((data(ep,1)-1)/(int_m))+1;
    DTamb(site).amb(nep,prn)=result_data;
    PAelev(site).elev(nep,prn)=com_p;
  end
end
amb=DTamb;
P=PAelev;
end