function [upd_w_m]=UD_estimate_upd_tq_fix_amb(FTamb43,upd_wm,ep,prn,site,cfg)
%% wide lane ambiguity extraction
[a,nsite_tc]=size(upd_wm(ep).amb);
upd_w_m=9999.999;
for s1=1:nsite_tc
    if(strcmp(FTamb43(site).name,upd_wm(ep).amb(s1).name))
        upd_w_m=upd_wm(ep).amb(s1).amb(ep,prn);
        break;
    end
end
%% determine whether the current width ambiguity is fixed
% if(rem(upd_w_m,1)~=0)
%     upd_w_m=9999.999; % un-fixed
% end
end