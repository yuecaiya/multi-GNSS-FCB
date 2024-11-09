function [DTamb23,DTamb12,DTamb43,PAelev23,PAelev12,PAelev43]=UD_combine_amb_arc(DNamb23,DNamb12,DNamb43,PNelev23,PNelev12,PNelev43,gps_slip,site_num_fcb,cfg)
%% merge within arc segments based on variance and cycle slip information
% args:
%     DNamb23 : UWL combination ambiguity
%     DNamb12 :  WL combination ambiguity
%     DNamb43 :  NL combination ambiguity
%     PNelev23 : weight for UWL
%     PNelev12 : weight for  WL
%     PNelev43 : weight for  NL
%     gps_slip    : cycle slip\
%     site_num_fcb: site number
%     cfg         : configuration information
% return:
%     DTamb23 : UWL ambiguity after merging arc segments
%     DTamb12 :  WL ambiguity after merging arc segments
%     DTamb43 :  NL ambiguity after merging arc segments
%     PAelev23: weight for UWL after merging arc segments
%     PAelev12: weight for  WL after merging arc segments
%     PAelev43: weight for  NLafter merging arc segments
% made by Caiya Yue @ CUMTB and CASM
% ----
%% merge within arc segments, and weighting based on STD
time=cfg.span_t*2;
elev=cfg.elev;
yz=cfg.comb_time*(60/cfg.inter);
lc=length(DNamb23(1).amb(1,:));
stda=zeros(3000*40,lc);
kk=zeros(lc,1);
for i=1:site_num_fcb
    ln=length(DNamb23(i).amb(:,1));
    lc=length(DNamb23(i).amb(1,:));
    DTamb23(i).name=DNamb23(i).name;
    tmp_name=DNamb23(i).name(1:4);
    fprintf('     %s ---%3dst\n',DTamb23(i).name(1:4),i);
    DTamb12(i).name=DNamb12(i).name;
    DTamb43(i).name=DNamb43(i).name;
    lnn=ln/time;
    DTamb23(i).amb=zeros(lnn,lc);
    DTamb12(i).amb=zeros(lnn,lc);
    DTamb43(i).amb=zeros(lnn,lc);
    PAelev23(i).elev=zeros(lnn,lc);
    PAelev12(i).elev=zeros(lnn,lc);
    PAelev43(i).elev=zeros(lnn,lc);    
    % merge based on variance, cycle jumps, and quality control strategies
    for ii=1:lc
      %% extract data for ultra wide lane
       data=zeros(ln,3); % time N23 Pelev23
       k=0;
       for j=1:ln
           % extract continuous data without cycle slip
           if(PNelev23(i).elev(j,ii)~=0 && gps_slip(i).amb(j+1,ii+1)~=1 && j~=2880)
               k=k+1;
               data(k,1)=j;
               data(k,2)=DNamb23(i).amb(j,ii);
               data(k,3)=PNelev23(i).elev(j,ii);               
           else
              % merge data
               if(k<=yz) 
                   k=0;
                   data=zeros(ln,3); % time N23 Pelev23
                   continue;
               else
                   [amb23,P23]=UD_combine_amb_sub_arc(i,ii,k,data,cfg,'u',DTamb23,PAelev23); 
                   DTamb23=amb23;
                   PAelev23=P23;
                   k=0;
                   data=zeros(ln,3); % time N23 Pelev23
               end
           end
       end
      %% extract data for ultra wide lane
       data=zeros(ln,3); % time N23 Pelev23
       k=0;
       for j=1:ln
           % extract continuous data without cycle slip
           if(PNelev12(i).elev(j,ii)~=0 && gps_slip(i).amb(j+1,ii+1)~=1 && j~=2880)
               k=k+1;
               data(k,1)=j;
               data(k,2)=DNamb12(i).amb(j,ii);
               data(k,3)=PNelev12(i).elev(j,ii);               
           else
              % merge data
               if(k<=yz)
                   k=0;
                   data=zeros(ln,3); % time N23 Pelev23
                   continue;
               else
                   [amb12,P12]=UD_combine_amb_sub_arc(i,ii,k,data,cfg,'w',DTamb12,PAelev12);
                   DTamb12=amb12;
                   PAelev12=P12;
                   k=0;
                   data=zeros(ln,3); % time N23 Pelev23
               end
           end
       end
       %% extract data for ultra wide lane
       data=zeros(ln,3); % time N23 Pelev23
       k=0;
       for j=1:ln
           % extract continuous data without cycle slip
           if(PNelev43(i).elev(j,ii)~=0 && gps_slip(i).amb(j+1,ii+1)~=1 && j~=2880)
               k=k+1;
               data(k,1)=j;
               data(k,2)=DNamb43(i).amb(j,ii);
               data(k,3)=PNelev43(i).elev(j,ii); 
           else
              % merge data
               if(k<yz)
                   k=0;
                   data=zeros(ln,3); % time N23 Pelev23
                   continue;
               else
                   [amb43,P43]=UD_combine_amb_sub_arc(i,ii,k,data,cfg,'n',DTamb43,PAelev43); 
                   DTamb43=amb43;
                   PAelev43=P43;
                   k=0;
                   data=zeros(ln,3); % time N23 Pelev23
                   kk(ii)= kk(ii)+1;
%                    stda(kk(ii),ii)=std(data(k-100:k,2));
%                    k=0;
%                    data=zeros(ln,3); % time N23 Pelev23
               end
           end
       end
    end   
end
%% output the total number of arc segments observed by each satellite
for s=1:lc
   ss=num2str(s);
   if(s<10); ss=strcat('0',num2str(s)); end
   fprintf('PRN%s %3d\n',ss,kk(s));
end
end