function [num_sats]=UD_total_sat_num(gps_fcb1,site_num_fcb)
%% Reference satellite adaptive judgment
% args:
%     gps_fcb1    : ambiguity
%     site_num_fcb: site number
% return:
%     num_sats    : reference satellite
% made by Caiya Yue @ CUMTB and CASM
% ----
ln=length(gps_fcb1(1).amb(:,1));
lc=length(gps_fcb1(1).amb(1,:));
num_sats=zeros(lc-1,1);
for i=2:lc 
    for j=1:site_num_fcb   
        for k=2:ln
           if(gps_fcb1(j).amb(k,i)~=0); num_sats(i-1)=num_sats(i-1)+1; end
        end
    end
end