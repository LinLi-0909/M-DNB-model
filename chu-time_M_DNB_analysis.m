clc
close all
ex=importdata('9606.protein.links.v10.txt');
links=ex.textdata(2:end,1:2);
k=1;
whole{1,1}=links{1,1};
whole{1,2}{1,1}=links{1,2};
j=1;
for i=1:size(links,1)-1
    if strcmpi(links{i,1},links{i+1,1})
        whole{k,2}{j+1,1}=links{i+1,2};
        j=j+1;
    else
        j=1;
        k=k+1;
        whole{k,1}=links{i+1,1};
        whole{k,2}{1,1}=links{i+1,2};
    end
end
whole(:,3)=importdata('towhole.csv');
clear ex;clear links;
alldata=importdata('GSE75748_sc_time_course_ec.csv');
data=alldata.data;
feature=alldata.textdata(2:end,1);
timeIdx=[{'Time1'},{[1:92]};
         {'Time2'},{[93:194]};
         {'Time3'},{[195:260]};
         {'Time4'},{[261:432]};
         {'Time5'},{[433:570]};
         {'Time6'},{[571:758]}];
CI = Get_CI(data,time_Idx,feature,whole);
allCI=[CI{1}(:,1),CI{2}(:,1),CI{3}(:,1),CI{4}(:,1),CI{5}(:,1),CI{6}(:,1)];
[b,c]=sort(allCI(:),'descend');
[m,n]=find(allCI>=b(100),100);
m=unique(m);
topCI=allCI(m,:);
usefulgene=modules(m,1);

%apply to compare with different gene
[data2,txt]=xlsread('table s3.xlsx',1);
difffeature=txt(2:end,1);
clear txt;clear data2
TF=ismember(usefulgene,difffeature);
x=[1:size(m)];y=[1:6];
mesh(y,x,allCI(m,:))
[top50CI,QI]=Get_Critical_Indicators(timeIdx,CI,50);


%plot figure
    lw = 1.2;
    mk = '.';
    ms = 16;
    fs = 14;
    % DNB behavior
figure('Name','selections behavior');
plot(avetopCI,'-k','linewidth',lw,'marker',mk,'MarkerSize',ms)
hold on
plot([2,4],[avetopCI(2,1),avetopCI(4,1)],'*r','MarkerSize',16)
set(gca,'xlim',[1,size(timeIdx,1)],'xtick',1:size(timeIdx,1))
ylabel('Mean of CI in DNB','FontSize',16)
box on;
xlabel('Progression \rightarrow','FontSize',16);
title('Quantitative Indicator','FontSize',18)
