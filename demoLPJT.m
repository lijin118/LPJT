clear;
close all;
clc;
warning off;
%%%%% addpath %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%dataset settings%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Dataset Setting!!!\n');
srcStrSURF12 = {'Caltech10','Caltech10','Caltech10','amazon',   'amazon','amazon','webcam',  'webcam', 'webcam','dslr',    'dslr',   'dslr'};
tgtStrSURF12 = {'amazon',   'webcam',   'dslr',     'Caltech10','webcam','dslr',  'Caltech10','amazon','dslr',  'Caltech10','amazon','webcam'};

srcStrDecaf12 = {'Caltech10','Caltech10','Caltech10','amazon','amazon','amazon','webcam','webcam','webcam','dslr','dslr','dslr'};
tgtStrDecaf12 = {'amazon','webcam','dslr','Caltech10','webcam','dslr','Caltech10','amazon','dslr','Caltech10','amazon','webcam'};

srcStrUSPS = {'MNIST','USPS'};
tgtStrUSPS = {'USPS','MNIST'};

srcStrPie = {'PIE05','PIE05','PIE05','PIE05','PIE07','PIE07','PIE07','PIE07','PIE09','PIE09','PIE09','PIE09','PIE27','PIE27','PIE27','PIE27','PIE29','PIE29','PIE29','PIE29'};
tgtStrPie = {'PIE07','PIE09','PIE27','PIE29','PIE05','PIE09','PIE27','PIE29','PIE05','PIE07','PIE27','PIE29','PIE05','PIE07','PIE09','PIE29','PIE05','PIE07','PIE09','PIE27'};

dataset = 'SURF12';
% dataset = 'DECAF12';
% dataset = 'USPS';
% dataset = 'PIE';

disp(dataset);
if(strcmp(dataset,'SURF12'))
	srcStr = srcStrSURF12;
	tgtStr = tgtStrSURF12;
	DataLen = 12;
elseif(strcmp(dataset,'DECAF12'))
	srcStr = srcStrDecaf12;
	tgtStr = tgtStrDecaf12;
	DataLen = 12;
elseif(strcmp(dataset,'USPS'))
	srcStr = srcStrUSPS;
	tgtStr = tgtStrUSPS;
	DataLen = 2;
elseif(strcmp(dataset,'PIE'))
    srcStr = srcStrPie;
    tgtStr = tgtStrPie;
    DataLen = 20;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Reading Dataset!!!\n');
datapath = '../data/';
TransferStr = repmat(struct('str',''),1,DataLen);
DataArr = repmat(struct('Xs',[],'Xt',[],'Ys',[],'Yt',[]),1,DataLen);
for iData = 1:DataLen
    src = char(srcStr{iData});
    tgt = char(tgtStr{iData});
	TransferStr(iData).str = strcat(src,'  --->  ',tgt);
	if(strcmp(dataset,'SURF12'))       
        load([datapath src '_SURF_L10.mat']);
        Xs = fts ./ repmat(sum(fts,2),1,size(fts,2)); 
        Ys = labels;
        Xs = zscore(Xs);
        Xs = normr(Xs)';

		load([datapath tgt '_SURF_L10.mat']);
        Xt = fts ./ repmat(sum(fts,2),1,size(fts,2)); 
        Yt = labels;
        Xt = zscore(Xt);
        Xt = normr(Xt)';
             
	elseif(strcmp(dataset,'DECAF12'))
        load([datapath src '_decaf.mat']);
        Xs = fea ./ repmat(sum(fea,2),1,size(fea,2)); 
        Ys = gnd;
        Xs = zscore(Xs);
        Xs = normr(Xs)';

        load([datapath tgt '_decaf.mat']);
        Xt = fea ./ repmat(sum(fea,2),1,size(fea,2)); 
        Yt = gnd;
        Xt = zscore(Xt);
        Xt = normr(Xt)';
                
	elseif(strcmp(dataset,'USPS'))
		load([datapath src '_vs_' tgt '.mat']);
        Xs = X_src;
        Xt = X_tar;
        Xs = normc(Xs);                       
		Ys = Y_src;
        Xt = normc(Xt);
		Yt = Y_tar;
    elseif(strcmp(dataset,'PIE'))
        load([datapath src '.mat']);
        Xs = fea ./ repmat(sum(fea,2),1,size(fea,2)); 
        Ys = gnd;
        Xs = zscore(Xs);
        Xs = normr(Xs)';
        
        load([datapath tgt '.mat']);
        Xt = fea ./ repmat(sum(fea,2),1,size(fea,2)); 
        Yt = gnd;
        Xt = zscore(Xt);
        Xt = normr(Xt)';
    end
    DataArr(1,iData).Xs = Xs;
    DataArr(1,iData).Xt = Xt;
    DataArr(1,iData).Ys = Ys;
    DataArr(1,iData).Yt = Yt;
end

fprintf('Options Setting!!!\n');
options =[];
options.T = 5;
options.interK = 10; 
options.intraK = 5; 
options.ReducedDim = 25; 
options.delta = 0.5;  
options.mu = 1;  
options.lambda = 1; 
options.gamma = 0.05;  
options.nu = 0.5; 
options.dim = 40;

fprintf('Executing!!!\n');
iterDataset = 1:DataLen;
iterLen = length(iterDataset);
for iData = iterDataset 
    warning off;
    fprintf('\nIteration Count: %d/%d\n',iData,iterLen);
    disp(TransferStr(iData).str);                                
    if (strcmp(dataset,'SURF12') || strcmp(dataset,'DECAF12'))
        if iData > 9
           options.nu = 0.3; 
        elseif iData < 10
            options.nu = 0.5;                                               
        end                                            
    end
    [accVec]=LPJT(DataArr(1,iData),options);
end                               
fprintf('\n');
clear DataArr options Results srcStr tgtStr