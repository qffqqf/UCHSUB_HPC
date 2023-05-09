clear all
clc
addpath([pwd,'/Functions'])
addpath([pwd,'/plate_models'])
addpath([pwd,'/input_model_normal'])
warning('off','all')

tStart = tic;
%% Load
model_file = 'input_model_10x10.mat';
load(model_file)
savename = 'vbplate_10x10_hbms';

%% Analysis setup
options.nModeI = 20;
options.nModeA = 30;
options.nEig = 1000;
options.exppt = 2000;

%% Hierarchical reduction
uz_hbms = hss(uc_model, param, tree_model, options);

figure
freq = 0:2:1000;
semilogy(freq, abs(uz_hbms))

timing = toc(tStart);
save([savename,'.mat'],'freq','uz_hbms','timing')



