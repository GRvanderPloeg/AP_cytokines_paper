% PARAFAC functionality
addpath("../N-way shell/Scripts/"); % own scripts
addpath("../N-way shell/N-way toolbox/"); % from Rasmus Bro

%%
% Load raw cytokines data
cytokines_data = readmatrix(".\input_deduplicated_RvdP.csv", Filetype="delimitedtext", Delimiter=" ");
cytokines_meta = readmatrix(".\input_deduplicated_metadata_RvdP.csv", Filetype="delimitedtext", Delimiter=" ", OutputType="string");

onlyCases = cytokines_meta(:,6) == "control";
cytokines_data = cytokines_data(onlyCases,:);
cytokines_meta = cytokines_meta(onlyCases,:);

% NOTE: LOG TRANSFORM WAS REMOVED BECAUSE OF DERIVATIVE COMPUTATION.
% Yes, the data is already normally distributed because of the derivative.

%%
% Reshape to 3-way matrix
num_timepoints = 6;
keepIndividuals = true;
cytokines_cube = rawDataToCube(cytokines_data, cytokines_meta(:,1), cytokines_meta(:,2), keepIndividuals);

%%
% Convert to derivative data
cytokines_derivative = cytokines_cube;

for i=1:(num_timepoints-1)
    cytokines_derivative(:,:,i) = cytokines_cube(:,:,i+1) - cytokines_cube(:,:,i);
end

cytokines_derivative(:,:,6) = [];

%%
% Center and scale
[cytokines_cnt, cytokines_means] = centerData(cytokines_derivative, 1);
[cytokines_cnt_scl, cytokines_stds] = scaleData(cytokines_cnt, 2);

%%
% Might want to add editing the metadata here
[I,J,K] = size(cytokines_cnt_scl);

subjectMeta = cytokines_meta(:,[1 6]);
subjectMeta = sortrows(subjectMeta, 1);
subjectMeta = unique(subjectMeta, "rows", "stable"); % stable stops resorting rows

%subjectMeta(subjectMeta(:,2) == "A", 2) = "Asymptomatic";
%subjectMeta(subjectMeta(:,2) == "S", 2) = "Symptomatic";
%subjectMeta(subjectMeta(:,2) == "C", 2) = "Control";

%featureMeta = [1:J; ones(1,J,1)]';
featureMeta = ["VEGF" "CRP" "GM-CSF" "IL1alpha" "IL1beta" "IL4" "IL6" "IL8" "IL10" "IL12p70" "IL17A" "IFNgamma" "MIP1alpha" "OPG" "TNFalpha" "RANKL"]';

%%
% Initialize PARAFAC options
Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)
const = [0 0 0];

%%
% Filter out outliers
df_cnt_scl_filtered = cytokines_cnt_scl;
subjectMeta_filtered = subjectMeta;
featureMeta_filtered = featureMeta;

% Remove outliers
df_cnt_scl_filtered(:,[4],:) = [];
featureMeta_filtered([4],:) = [];

df_cnt_scl_filtered([13 14],:,:) = [];
subjectMeta_filtered([13 14],:) = [];

%df_cnt_scl_filtered(:,:,5) = [];


%%
% Run PARAFACs
path_start = "./test_run/Figures/";
maxComponents=5;
numReps=25;
maxIterations=20;
timeMeta = (1:5)';
metaData = {subjectMeta_filtered, featureMeta_filtered, timeMeta};
legendIndex = [0 1 0];
numColsPerLegend = [0 3 0];
xlabels = ["Subject index", "Feature index", "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
resort = [true false false];

[Models, Cons, VarExps, Boots, BootVarExps, Tuckers] = quickReport(df_cnt_scl_filtered, Options, const, maxComponents, numReps, maxIterations, metaData, resort, legendIndex, numColsPerLegend, xlabels, titles, "Athina study bootstrapped", path_start+"Athina_cytokines");

%%
% Pick the number of components that seems appropriate based on the plots
% and pick the model that you like the best.
% Best practices: pick between best Corcondia or best variance explained.
% Alternative: choose a good tucker congruence model.

numFactors = 1;
%choice = find(VarExps{numFactors}==max(VarExps{numFactors}));
choice = find(Cons{numFactors}==max(Cons{numFactors}));
Model = pickModel(Models, numFactors, choice);

%%
% Save the models
model_path = "./20230905_run_firstDerivative_onlyControl/PARAFAC models/";

annotatedModel = annotateModel(df_cnt_scl_filtered, Model, metaData);
savePARAFAC(df_cnt_scl_filtered, Model, annotatedModel, model_path + "Tongue");

%%
% Plot PARAFAC model
legendIndex = [0 0 0];
numColsPerLegend = [0 0 0];
xlabels = ["Subject index", "Feature index", "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
resort = [true true false];
path_start = "./20230905_run_firstDerivative_onlyControl/Figures/";

plotPARAFAC4(annotatedModel, numFactors, VarExps, choice, resort, legendIndex, numColsPerLegend, xlabels, titles, "PARAFAC Athina", path_start + "PARAFAC_Athina.jpg");
