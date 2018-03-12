%This script creates the variables through which the required parameters 
%and files are inputted to the metagenomic pipeline (MgPipe). Input 
%variables should be changed by the user according to what specified in the 
%documentation. Running this script will automatically launch the pipeline. 

% Federico Baldini, 2017-2018


%REQUIRED INPUT VARIABLES
modPath='Y:\Federico\Eldermet\191017\models\'; %path to microbiota models
infoPath='P:\Documenti\GitLab\Microbiome_Toolbox\MgPipe\examples\'; %path to the information files
resPath='Y:\Federico\testingMgPipe\fede\fede2\'
abunFileName='normCoverage.csv'
resPath='Y:\Federico\testingMgPipe\fede\'; %path to results directory 
objre={'EX_biomass(e)'}; %name of objective function of organisms
sDiet='EUAverageDiet' %standard diet type
figForm = '-depsc' %the output is vectorized picture, change to '-dpng' for .png
nWok = 3; %number of cores dedicated for parallelization 
autoFix = 1 %autofix for names mismatch
compMod = 0; % if outputs in open formats should be produced for each section (1=T)
patStat = 0; %if documentations on patient health status is provided (0 not 1 yes)
rDiet = 0 %to enable also rich diet simulations 
extSolve = 0 
fvaType = 0;
%END OF REQUIRED INPUT VARIABLES

%%
%PIPELINE LAUNCHER 
[init,modPath,toolboxPath,resPath,dietFilePath,abunFilePath,objre,figForm,solver,numWorkers,autoFix,compMod,patStat,rDiet,extSolve,fvaType,autorun]= initMgPipe(modPath, toolboxPath, resPath, dietFilePath, abunFilePath, objre, figForm, solver, numWorkers, autoFix, compMod, patStat, rDiet,extSolve,fvaType,autorun);
