function [final_upd_w_day]=UD_compose_day_FCB(final_upd_w1,cfg)
%% Single day UPD/FCB synthesis
% made by Caiya Yue @ CUMTB and CASM
% ----
%% ---
[a,nep]=size(final_upd_w1);
[a,nsite]=size(final_upd_w1(1).amb);
lc=length(final_upd_w1(1).amb(1).amb(1,:));
for prn=1:lc
    for s=1:nsite
        % extract UPD/FCB
        upd_w=zeros(nep,1);
        k=0;
        for e=1:nep
            [a,nsite1]=size(final_upd_w1(e).amb);
            if(nsite1<s); break; end
            if(final_upd_w1(e).amb(s).amb(e,prn)~=0);
                k=k+1;
                upd_w(k)=final_upd_w1(e).amb(s).amb(e,prn);
            end
        end
        % best integer value
        if(k==0); continue; end
        upd_max=max(upd_w(1:k));
        upd_min=min(upd_w(1:k));
        upd_w1=[];
        kw1=0;
        upd_w2=[];
        kw2=0;
        for i=1:k
            if(abs(upd_w(i)-upd_max)<2)
                kw1=kw1+1;
                upd_w1(kw1,1)=upd_w(i);
            end
            if(abs(upd_w(i)-upd_min)<2)
                kw2=kw2+1;
                upd_w2(kw2,1)=upd_w(i);
            end
        end
        med_w1=median(upd_w1);
        med_w2=median(upd_w2);
        % reassigning WL integer ambiguity
        for e=1:nep
            [a,nsite1]=size(final_upd_w1(e).amb);
            if(nsite1<s); break; end
            if(final_upd_w1(e).amb(s).amb(e,prn)~=0);
                if(abs(final_upd_w1(e).amb(s).amb(e,prn)-med_w1)<2); final_upd_w1(e).amb(s).amb(e,prn)=med_w1; end 
                if(abs(final_upd_w1(e).amb(s).amb(e,prn)-med_w2)<2); final_upd_w1(e).amb(s).amb(e,prn)=med_w2; end 
            end
        end
    end  
end
final_upd_w_day=final_upd_w1;
%% include data between [-1 1] by adding or subtracting one week
% final_upd_w_day=zeros(lc,1);
% for s=1:lc
%     % [-1 1]
%     upd_w1=upd_w(:,s);
%     for n=1:ln
%        if(upd_w1(n)>800); continue; end
%        for i=1:100
%           if(abs(upd_w1(n))<1)
%               break;
%           else
%               if(upd_w1(n)<-1); upd_w1(n)=upd_w1(n)+1; end
%               if(upd_w1(n)> 1); upd_w1(n)=upd_w1(n)-1; end
%           end
%        end
%     end
%     % based on the median
%     med_u=median(upd_w1);
%     for n=1:ln
%         if(upd_w1(n)>800); continue; end
%         if((upd_w1(n)-med_u)>0.5); upd_w1(n)=upd_w1(n)-1; end
%         if((upd_w1(n)-med_u)<-0.5); upd_w1(n)=upd_w1(n)+1; end
%     end
%     final_upd_w_day(s)=median(upd_w1);
% end
end