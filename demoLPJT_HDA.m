clear;
close all;
clc;
warning off;
%%%%% addpath %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Dataset Setting!!!\n');
srcStr = {'Caltech10','Caltech10','Caltech10','amazon',   'amazon',  'amazon',  'webcam',  'webcam', 'webcam'};
tgtStr = {'Caltech10', 'amazon',   'webcam',   'Caltech10', 'amazon', 'webcam', 'Caltech10','amazon','webcam'};
surfPostFix = '_SURF_L10';
decafPostFix = '_decaf';
vggPostFix = '_VGG-FC6';

srcPostFix = vggPostFix;
tgtPostFix = decafPostFix;

fprintf('\n%s -------> %s\n',srcPostFix,tgtPostFix);
DataLen = 9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Reading Dataset!!!\n');
srcDataPath = '../source/';
tgtDataPath = '../target/';
TransferStr = repmat(struct('str',''),1,DataLen);
DataArr = repmat(struct('Xs',[],'Xu',[],'Ys',[],'Yu',[],'Xl',[],'Yl',[]),1,DataLen);

number = zeros(10,1);
cumuNumber = zeros(10,1);
for iData = 1:DataLen
    src = char(srcStr{iData});
    tgt = char(tgtStr{iData});
	TransferStr(iData).str = strcat(src,'  --->  ',tgt);   
    
    load([srcDataPath,srcPostFix,'/',src,'/1.mat']);
    Xs = sdata ./ repmat(sqrt(sum(sdata.^2,2)),1,size(sdata,2));
    Ys = slbl;           
    
    load([tgtDataPath,tgtPostFix,'/',tgt,'/1.mat']);
    Xl = ldata ./ repmat(sqrt(sum(ldata.^2,2)),1,size(ldata,2));
    Yl = llbl;
    Xu = udata ./ repmat(sqrt(sum(udata.^2,2)),1,size(udata,2));
    Yu = ulbl;
    
    Xs = Xs';
    Xl = Xl';
    Xu = Xu';
  
    DataArr(1,iData).Xs = Xs;
    DataArr(1,iData).Xu = Xu;
    DataArr(1,iData).Ys = Ys;
    DataArr(1,iData).Yu = Yu;
    DataArr(1,iData).Xl = Xl;
    DataArr(1,iData).Yl = Yl;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Options Setting!!!\n');
options = [];
options.T =5; 
options.nu = 0.5; 
options.ReducedDim = 120;
options.interK = 5; 
options.intraK = 5;
options.delta = 0.5; 
options.gamma = 0.00001;  
options.mu = 0.00001;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Executing!!!\n');
iterDataset = 1:DataLen;
iterLen = length(iterDataset);

for iData = iterDataset
    warning off;
    fprintf('\n------->Iteration Count: %d/%d\n',iData,iterLen);
    disp(TransferStr(iData).str);                                
    [accVec]=LPJT_HDA(DataArr(1,iData),options, 0);
end                         
fprintf('\n');