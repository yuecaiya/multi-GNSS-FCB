function [FDTamb23,FDTamb12,FDTamb43]=UD_combine_Fra_amb(DTamb23,DTamb12,DTamb43,site_num_fcb,cfg)
%% extract the merged floating ambiguity fractional part
% args:
%     DTamb23 : UWL combination ambiguity
%     DTamb12 :  WL combination ambiguity
%     DTamb43 :  NL combination ambiguity
%     site_num_fcb: site number
%     cfg         : configuration information
% return:
%     FDTamb23 : UWL ambiguity fractional part
%     FDTamb12 :  WL ambiguity fractional part
%     FTamb43 :   NL ambiguity fractional part
% made by Caiya Yue @ CUMTB and CASM
% ----
FDTamb23=DTamb23;
FDTamb12=DTamb12;
FDTamb43=DTamb43;
for i=1:site_num_fcb
    ln=length(DTamb23(i).amb(:,1));
    lc=length(DTamb23(i).amb(1,:));
    for ep=1:ln
        for prn=1:lc

            fcb=DTamb23(i).amb(ep,prn)-fix(DTamb23(i).amb(ep,prn));% UWL£¨0 1 -1£©
            FDTamb23(i).amb(ep,prn)=fcb; %-------- UWL£¨0 1 -1£©
            fcb=DTamb12(i).amb(ep,prn)-fix(DTamb12(i).amb(ep,prn));%-------- WL£¨1 -1 0£©
            FDTamb12(i).amb(ep,prn)=fcb; %-------- WL£¨1 -1 0£©
            % Raw FCB
            if(cfg.FCB_MOD==1)
                fcb=DTamb43(i).amb(ep,prn)-fix(DTamb43(i).amb(ep,prn));%-------- NL£¨4 -3 0£©
                FDTamb43(i).amb(ep,prn)=fcb; %--------NL£¨4 -3 0£©
            end
            % UWL-WL-NL
            if(cfg.FCB_MOD==2)
                FDTamb43(i).amb(ep,prn)=DTamb43(i).amb(ep,prn); %--------IF combination
            end
        end        
    end
end