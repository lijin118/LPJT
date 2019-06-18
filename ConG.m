function [Ww,Wb] = ConG(gnd,options,fea,method)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Label_b = unique(gnd(:,1));
nLabel_b = length(Label_b);
Label_s = unique(gnd(:,2));
nLabel_s = length(Label_s);

nSmp=length(fea(:,1));

Ww=zeros(nSmp);
Wb=zeros(nSmp);

if(options.k>nLabel_b-1)
    options.k=nLabel_b-1;
end

options.gnd=gnd(:,1);
%%%%% construct Wb
switch lower(method)
    case {lower('CGE')}
        %within sub-cluster
        for idx=1:nSmp
            fea_Indx=intersect(find(gnd(:,1)==gnd(idx,1)), find(gnd(:,2)==gnd(idx,2)));
            if(idx==fea_Indx)
                continue;
            else
                D=EuDist2(fea(idx,:),fea(fea_Indx,:));
                [dump,idx_s]=sort(D,2);
                %disp '111';
                Ww(idx,fea_Indx(idx_s(1,2)))=1;%exp(-dump(1,2)/(2*options.t^2));
            end
            clear D dump fea_Indx
            
            %between sub-cluster
            for i=1:nLabel_s
                if(Label_s(i)==gnd(idx,2))
                    continue;
                else
                    %[idx i]
                    fea_Indx=intersect(find(gnd(:,1)==gnd(idx,1)), find(gnd(:,2)==Label_s(i)));
                    if(isempty(fea_Indx))
                        continue;
                    else
                        D=EuDist2(fea(idx,:),fea(fea_Indx,:));
                        %[dump,idx_s]=sort(D,2);
                        [dump,idx_s]=min(D);
                        Ww(idx,fea_Indx(idx_s,1))=1;%exp(-dump/(2*options.t^2));
                        clear D dump fea_Indx idx_s
                    end
                end
            end
        end
        
        %between class
        for s=1:nLabel_s
            for b=1:nLabel_b
                curr_Indx=intersect(find(gnd(:,1)==Label_b(b)), find(gnd(:,2)==Label_s(s)));
                if(isempty(curr_Indx))
                    continue;
                else
                    for j=1:nLabel_b
                        if(j==b)
                            continue;
                        else
                            com_Indx=intersect(find(gnd(:,1)==Label_b(j)), find(gnd(:,2)==Label_s(s)));
                            D=EuDist2(fea(curr_Indx,:),fea(com_Indx,:));
                            dump=min(min(D));
                            [row,col]=find(D==dump);
                            %com_Indx
                            if(isempty(com_Indx))
                                continue;
                            else
                                Wb(curr_Indx(row,1),com_Indx(col,1))=exp(-dump/(2*options.t^2));
                                clear D dump com_Indx row col
                            end
                        end
                    end
                end
            end
        end
        
    case {lower('LFDA')}
        Ww = constructW(fea,options);
        for idx=1:nSmp
            fea_Indx=find(gnd(:,1)~=gnd(idx,1));
            D=EuDist2(fea(idx,:),fea(fea_Indx,:));
            Wb(idx,fea_Indx)=exp(-D/(2*options.t^2));
            clear D fea_Indx;
        end
        clear idx;
    case {lower('SRGE')}
        Ww = constructW(fea,options);
        tmp=zeros(nLabel_b,3);
        for b=1:nLabel_b
            f1=find(gnd(:,1)==Label_b(b));
            for j=1:nLabel_b
                f2=find(gnd(:,1)==Label_b(j));
                if(j==b)
                    tmp(j,:)=[0 0 0];
                    continue;
                else
                    D=EuDist2(fea(f1,:),fea(f2,:));
                    %[dump,idx_s] = sort(D,2); % sort each row
                    dump=min(min(D));
                    [row,col]=find(D==dump);
                    tmp(j,:)=[dump f1(row,1) f2(col,1)];
                    clear dump f2 D
                end
            end
            tmp=sortrows(tmp,1);
            for k=2:options.k+1
                Wb(tmp(k,2),tmp(k,3))=exp(-tmp(k,1)/(2*options.t^2));
                clear dump row col
            end
            clear f1 tmp k
        end
        
    case {lower('LDE')}
        Ww = constructW(fea,options);
        for idx=1:nSmp
            fea_Indx=find(gnd(:,1)~=gnd(idx,1));
            D=EuDist2(fea(idx,:),fea(fea_Indx,:));
            [dump,idx_s] = sort(D,2); % sort each row
            clear D;
            idx_s = idx_s(1,2:options.k+1);
            dump = dump(1,2:options.k+1);
            dump = exp(-dump/(2*options.t^2));
            % for i=1:size(idx_s,2)
            %     Wb(idx,fea_Indx(idx_s(1,i)))=dump;
            %  end
            Wb(idx,fea_Indx(idx_s))=dump;
            clear  idx_s dump fea_Indx;
        end
        clear idx;
    otherwise
        error('not supported method!')
end


Ww = sparse(max(Ww,Ww'));
Wb = sparse(max(Wb,Wb'));

Db = full(sum(Wb,2));
Wb = -Wb;
for i=1:size(Wb,1)
    Wb(i,i) = Wb(i,i) + Db(i);
end

Dw = full(sum(Ww,2));
Ww = -Ww;
for i=1:size(Ww,1)
    Ww(i,i) = Ww(i,i) + Dw(i);
end

end

