function [x0]=UD_estimate_UPD_x0(v,B,P,prn_flag,rec_flag,ref_prn,nsite,FTamb23,cfg)
%% UPD/FCB initial value determination
% ---
%
%%
lnv=length(v);
lnc=length(prn_flag(:,1));
%% determine the receiver-end FCB from the reference satellite
rec_sat.rfcb=zeros(nsite,1);
for s=1:nsite
    rec_sat(s).rfcb=0;
    rec_sat(s).rflag=0;
    rec_sat(s).sfcb=zeros(lnc,2);
end
for i=1:lnv
   if(B(i,nsite+(ref_prn-1)+1)==-1)
      for j=1:nsite
         if(B(i,(j-1)+1)==1)
             rec_sat(j).rfcb=v(i); % FCB
             rec_sat(j).rflag=1; % update flag 
         end
      end
   end
end
%% Determine the satellite FCB according to the receiver end FCB
for s=1:nsite
   if(rec_sat(s).rflag==1) 
       for p=1:lnc
           if(B((s-1)*lnc+(p-1)+1,nsite+(p-1)+1)==-1)
               rec_sat(s).sfcb(p,1)=rec_sat(s).rfcb-v((s-1)*lnc+(p-1)+1);
               rec_sat(s).sfcb(p,2)=1;
           end
       end
   end
end
%% calculate mean value for the common satellite of multiple stations
ex_flag=0;
i_kt=0;
while (1)
    i_kt=i_kt+1;
    if(i_kt>nsite); break; end
    for p=1:lnc
        % extract of the common satellite FCB
        tmp_fcb=[];
        kt=0;
        for s=1:nsite
            if(rec_sat(s).sfcb(p,2)==1)
                kt=kt+1;
                tmp_fcb(kt,1)=rec_sat(s).sfcb(p,1);
            end
        end
        % calculate mean and reassignment 
        if(kt==0); continue; end
        ktt=0;
        tmp_fcb=r_tongyi_qj(tmp_fcb);
        med_fcb=median(tmp_fcb);
        for s=1:nsite
            if(B((s-1)*lnc+(p-1)+1,nsite+(p-1)+1)==-1)
                rec_sat(s).sfcb(p,1)=med_fcb;
                rec_sat(s).sfcb(p,2)=1;
                ktt=ktt+1;
            end
        end
    end
    %% calculate other stations' FCB based on the new satellite FCB
    for s=1:nsite
        if(rec_sat(s).rflag==1); continue; end
        tmp_fcb=[];
        kt=0;
        for p=1:lnc
            if(rec_sat(s).sfcb(p,2)==1)
                kt=kt+1;
                tmp_fcb(kt,1)=v((s-1)*lnc+(p-1)+1)+rec_sat(s).sfcb(p,1);
            end
        end
        if(kt==0); continue; end
        tmp_fcb=r_tongyi_qj(tmp_fcb);
        rec_sat(s).rfcb=median(tmp_fcb);
        rec_sat(s).rflag=1;
    end
    %% determine the satellite end FCB according to the receiver end FCB
    for s=1:nsite
       if(rec_sat(s).rflag==1) 
           for p=1:lnc
               if(B((s-1)*lnc+(p-1)+1,nsite+(p-1)+1)==-1)
                   rec_sat(s).sfcb(p,1)=rec_sat(s).rfcb-v((s-1)*lnc+(p-1)+1);
                   rec_sat(s).sfcb(p,2)=1;
               end
           end
       end
    end
    %% exit condition
    kt=0;
    for s=1:nsite
        if(rec_sat(s).rflag==1); kt=kt+1; end
    end
    if(kt==nsite && ex_flag==1); break; end
    if(kt==nsite && ex_flag==0); ex_flag=1; end
end
%% update the receiver and satellite FCB
xr=[];
xs=[];
kr=0;
ks=0;
% satellite FCB
for p=1:lnc
    tmp_fcb=[];
    kt=0;
    for s=1:nsite
        if(rec_sat(s).sfcb(p,2)==1)
            kt=kt+1;
            tmp_fcb(kt,1)=rec_sat(s).sfcb(p,1);
        end
    end
    if(kt==0); continue; end
    tmp_fcb=r_tongyi_qj(tmp_fcb);
    ks=ks+1;
    xs(ks)=median(tmp_fcb);
end
% receiver FCB
for s=1:nsite
    tmp_fcb=[];
    kt=0;
    for p=1:lnc
        if(B((s-1)*lnc+(p-1)+1,nsite+(p-1)+1)==-1)
            kt=kt+1;
            tmp_fcb(kt,1)=rec_sat(s).sfcb(p,1)+v((s-1)*lnc+(p-1)+1);
        end
    end
    if(kt==0); continue; end
    tmp_fcb=r_tongyi_qj(tmp_fcb);
    kr=kr+1;
    xr(kr)=median(tmp_fcb);
end
x0=[xr';xs'];
end
%% data processing was implemented by adding or subtracting for 1 cycle
function tmp_fcb=r_tongyi_qj(tmp_fcb)
    ln=length(tmp_fcb);
    for i=1:ln
        if(tmp_fcb(i)>0.5)
            tmp_fcb(i)=tmp_fcb(i)-1;
        end
        if(tmp_fcb(i)<-0.5)
            tmp_fcb(i)=tmp_fcb(i)+1;
        end
    end
end