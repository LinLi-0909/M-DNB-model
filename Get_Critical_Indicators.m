function [topmCI,QI]=Get_Critical_Indicators(timeIdx,CI,m) 
% get top m Critical Indicators in each timepoint,e.g. m=50
topmCI=cell(size(timeIdx,1),1);
for t=1:size(timeIdx,1)
    [temp,idd]=sort(CI{t}(:,1),'descend');
    topmCI{t}(:,1)=temp(1:m,1);
    topmCI{t}(:,2)=idd(1:m);
    avetopCI(t,1)=mean(topmCI{t}(:,1));
end
QI= avetopCI;
end