function tmp_upd_uwn = UD_process_upd_smooth_median(tmp_upd,ln_step,pw)
%%
% time series, stride, weight
ln=length(tmp_upd);
if(ln<ln_step)
    tmp_upd_uwn=tmp_upd;
else
    tmp_upd_uwn=zeros(ln,1);
    for i=1:ln

    end
end
end