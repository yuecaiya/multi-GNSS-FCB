function tmp_upd_uwn_0= UD_estimate_upd_process_smooth(tmp_upd_uwn,stepi,pw)
ln=length(tmp_upd_uwn);
% std_uwn=std(tmp_upd_uwn);
% med_uwn=median(tmp_upd_uwn);
% for e=1:ln
%     if(abs(tmp_upd_uwn(e)-med_uwn)>1.5*std_uwn)
%         tmp_upd_uwn(e)=med_uwn;
%     end
% end
        
if(ln<=stepi)
    tmp_upd_uwn_0=tmp_upd_uwn;
else
    med_t=median(tmp_upd_uwn(1:stepi));
    tmp_upd_uwn_0=[];
    tmp_upd_uwn_0(1,1)=med_t;
    for i=2:ln-stepi;
        tmp_upd_uwn_0(i,1)=tmp_upd_uwn_0(i-1,1)*pw+median(tmp_upd_uwn(i:stepi+i))*(1-pw);
    end
    tstepi=stepi;
    k=0;
    for i=ln-tstepi+1:ln
        k=k+1;
        tmp_upd_uwn_0(i,1)=tmp_upd_uwn_0(i-1,1)*pw+median(tmp_upd_uwn(i:stepi+i-k))*(1-pw);
    end
end
end