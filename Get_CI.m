function [CI] = Get_CI(data,time_Idx,feature,whole)
%% deal with the data
for t=1:size(timeIdx,1)
    data(:,timeIdx{t,2})=zscore(data(:,timeIdx{t,2}),0,2); 
end
ssd=zeros(size(data,1),size(timeIdx,1));
for i=1:size(timeIdx,1)
    ssd(:,i)=std(data(:,timeIdx{i,2}),0,2);
end

%%calculate all pcc
i=0;
for j=1:size(feature,1)
    k=0;
    k=find(strcmp(feature{j,1},whole(:,3)));
    if k~=0
        i=i+1;
        idx(i,1)=j;
    end
end
allpcc=cell(size(timeIdx,1),1);
for t=1:size(timeIdx,1)
    expr_m=data(idx,timeIdx{t,2});
    allpcc{t,1}=abs(corrcoef(expr_m'));
    allpcc{t,1}(isnan(allpcc{t,1}))=0;
end
for i=size(allpcc{1,1},1):-1:1
    if allpcc{1,1}(i,:)<=0.02
        idx(i)=[];
        for t=1:size(timeIdx,1)
            allpcc{t,1}(i,:)=[];
            allpcc{t,1}(:,i)=[];
        end
    end
end
%% find modules
for i=1:size(idx,1)
    k=0;
    k=find(strcmp(feature{idx(i),1},whole(:,3)));
    modules{i,1}=feature{idx(i),1}; 
    modules{i,2}{1,1}=feature{idx(i),1};
    modidx{i,1}(1,1)=idx(i);
    modidx{i,1}(1,2)=find(strcmp(modules{i,2}{1,1},feature(idx)));
    n=1;
    p=0;
    for m=1:size(whole{k,2},1)
        p=find(strcmp(whole{k,2}{m,1},whole(:,1)));
        q=find(strcmp(whole{p,3},feature(idx)));
        if ~isempty(q)
            if allpcc{1,1}(q,i)<=0.2&&allpcc{2,1}(q,i)<=0.2&&allpcc{3,1}(q,i)<=0.2...
                    &&allpcc{4,1}(q,i)<=0.2&&allpcc{5,1}(q,i)<=0.2&&...
                    allpcc{6,1}(q,i)<=0.2
            else
                n=n+1;
                modules{i,2}{n,1}=whole{p,3};
                modidx{i,1}(n,1)=find(strcmp(modules{i,2}{n,1},feature));
                modidx{i,1}(n,2)=find(strcmp(modules{i,2}{n,1},feature(idx)));
            end
        end
    end
end

%{

for i=1:size(modules,1)
     modules{i,2}(cellfun('length',modules{i,2})==0)=[];
end

%}
%% get CI
CI=cell(size(timeIdx,1),1);
for t=1:size(timeIdx,1)
    for i=1:size(modules,1)
        %sd
        expr_m=modidx{i,1}(:,1);
        CI{t}(i,2)=sum(ssd(expr_m,t))/size(expr_m,1);
        
        % PCC_in
        expr_m=data(expr_m,timeIdx{t,2});
        temp=abs(corrcoef(expr_m'));
        temp(isnan(temp))=0;
        pccall=sum(temp(1,2:end));
        m=size(temp,1)-1;
        for j=2:size(modules{i,2},1)-1
            q=find(strcmp(modules{i,2}{j,1},whole(:,3)));
            TF=ismember(modules{i,2}(j+1:size(modules{i,2},1),1),whole{q,2});
            p=find(TF==1);
            if ~isempty(p)
               pccall=pccall+sum(temp(j,p+j));
               m=m+size(p,1);
            end
        end
        CI{t}(i,3)=pccall/m;
        
       %pcc_out
        pccall=0;
        m=0;
        for j=2:size(modules{i,2})
            p=find(strcmp(modules{i,2}{j,1},modules(:,1)));
            expr_m=setdiff(modidx{p,1}(:,2),modidx{i,1}(:,2));
            expr_m=[modidx{i,1}(j,2);expr_m];
            expr_m=data(expr_m,timeIdx{t,2});
            temp=abs(corrcoef(expr_m'));
            temp(isnan(temp))=0;
            pccall=pccall+sum(temp(1,2:end));
            m=m+size(temp,1)-1;
        end
        CI{t}(i,4)=pccall/m;
        CI{t}(i,1)=sqrt(size(modules{i,2},1))*CI{t}(i,2)*CI{t}(i,3)/CI{t}(i,4);  
    end
end
for t=1:size(timeIdx,1)
    CI{t}(isnan(CI{t}))=0;
    CI{t}(isinf(CI{t}))=0;
end
end	
