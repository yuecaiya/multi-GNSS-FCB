function [month,day]=doy2mon_day(yyyy,doy)
%% Date conversion
% made by Caiya Yue @ CUMTB and CASM
% ----
yy=str2num(yyyy);
d=str2num(doy);
data1=[1 32 60 91 121 152 182 213 244 274 305  335 366]; % February has 28 days 
data2=[1 32 61 92 122 153 183 214 245 275 306  336 367]; % February has 29 days
if(rem(yy,4)==0); data=data2; end
if(rem(yy,4)~=0); data=data1; end
for i=1:12
   if(d>=data(i) && d<data(i+1))
       m=i; 
       day1=d-data(i)+1;
   end 
end
if(m<9)
    month=strcat('0',num2str(m));
else
    month=num2str(m);
end
if(day1<9)
    day=strcat('00',num2str(day1));
elseif(day1<100 && day1>9)
    day=strcat('0',num2str(day1));
else
    day=num2str(day1);
end