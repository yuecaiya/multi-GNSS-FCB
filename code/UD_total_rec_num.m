function [num_recs]=UD_total_rec_num(gps_fcb1,site_num_fcb)
%% Reference stations and satellite adaptive judgment
% args:
%     gps_fcb1    : ambiguity
%     site_num_fcb: site number
% return:
%     num_sats    : reference satellite
% made by Caiya Yue @ CUMTB and CASM
% ----
ln=length(gps_fcb1(1).amb(:,1));
lc=length(gps_fcb1(1).amb(1,:));
num_recs=zeros(site_num_fcb,1);
for i=1:site_num_fcb    
    for j=2:lc
        for k=2:ln
           if(gps_fcb1(i).amb(k,j)~=0); num_recs(i)=num_recs(i)+1; end
        end
    end
end