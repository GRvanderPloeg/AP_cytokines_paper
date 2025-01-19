% PARAFAC functionality
addpath("../N-way-shell/Scripts/"); % own scripts
addpath("../N-way-shell/N-way toolbox/"); % from Rasmus Bro

%%
% Load raw cytokines data
cytokines_data = readmatrix(".\input_deduplicated_RvdP.csv", Filetype="delimitedtext", Delimiter=" ");
cytokines_meta = readmatrix(".\input_deduplicated_metadata_RvdP.csv", Filetype="delimitedtext", Delimiter=" ", OutputType="string");

onlyCases = cytokines_meta(:,6) == "case";
cytokines_data = cytokines_data(onlyCases,:);
cytokines_meta = cytokines_meta(onlyCases,:);

%%
% Attach pain/no-pain metadata to these case subjects
meta2 = readmatrix("./Data/Root_meta_data_parafac.txt", Filetype="delimitedtext", Delimiter="\t", OutputType="string");
meta2 = meta2(:,[2,8]);

allIndividuals = unique(cytokines_meta(:,1));
controls = unique(cytokines_meta(cytokines_meta(:,6) == "control", 1));

for i=1:length(allIndividuals)
    subjectID = allIndividuals(i);
    if ~ismember(subjectID, controls)
        status = meta2(meta2(:,1) == subjectID, 2);
        if length(status) > 1 % somehow there are duplicate entries
            status = status(1);
        end
        cytokines_meta(cytokines_meta(:,1) == subjectID, 6) = status;
    end
end

%% Check the distributions because I'm scared
subplot(4,4,1); histogram(cytokines_data(:,1));
subplot(4,4,2); histogram(cytokines_data(:,2));
subplot(4,4,3); histogram(cytokines_data(:,3));
subplot(4,4,4); histogram(cytokines_data(:,4));
subplot(4,4,5); histogram(cytokines_data(:,5));
subplot(4,4,6); histogram(cytokines_data(:,6));
subplot(4,4,7); histogram(cytokines_data(:,7));
subplot(4,4,8); histogram(cytokines_data(:,8));
subplot(4,4,9); histogram(cytokines_data(:,9));
subplot(4,4,10); histogram(cytokines_data(:,10));
subplot(4,4,11); histogram(cytokines_data(:,11));
subplot(4,4,12); histogram(cytokines_data(:,12));
subplot(4,4,13); histogram(cytokines_data(:,13));
subplot(4,4,14); histogram(cytokines_data(:,14));
subplot(4,4,15); histogram(cytokines_data(:,15));
subplot(4,4,16); histogram(cytokines_data(:,16));

%%
% Log transform
[I,J] = size(cytokines_data);
vectorizedX = reshape(cytokines_data, I*J, 1);
pseudocount = min(vectorizedX(vectorizedX>0));

cytokines_data_log = log(cytokines_data + pseudocount);

%%
% Reshape to 3-way matrix
% Normal data or log transformed data?

keepIndividuals = true;
cytokines_cube = rawDataToCube(cytokines_data_log, cytokines_meta(:,1), str2double(cytokines_meta(:,2)), keepIndividuals);
%cytokines_cube = rawDataToCube(cytokines_data, cytokines_meta(:,1), str2double(cytokines_meta(:,2)), keepIndividuals);

%% Find the mean of t1 - 3 and remove that baseline
baseline = cytokines_cube(:,:,1) + cytokines_cube(:,:,2) + cytokines_cube(:,:,3);
baseline = baseline ./ 3;

cytokines_cube = cytokines_cube - baseline;
cytokines_cube = cytokines_cube(:,:,3:6);

%% Check distributions per feature
subplot(4,4,1); histogram(cytokines_cube(:,1,:));
subplot(4,4,2); histogram(cytokines_cube(:,2,:));
subplot(4,4,3); histogram(cytokines_cube(:,3,:));
subplot(4,4,4); histogram(cytokines_cube(:,4,:));
subplot(4,4,5); histogram(cytokines_cube(:,5,:));
subplot(4,4,6); histogram(cytokines_cube(:,6,:));
subplot(4,4,7); histogram(cytokines_cube(:,7,:));
subplot(4,4,8); histogram(cytokines_cube(:,8,:));
subplot(4,4,9); histogram(cytokines_cube(:,9,:));
subplot(4,4,10); histogram(cytokines_cube(:,10,:));
subplot(4,4,11); histogram(cytokines_cube(:,11,:));
subplot(4,4,12); histogram(cytokines_cube(:,12,:));
subplot(4,4,13); histogram(cytokines_cube(:,13,:));
subplot(4,4,14); histogram(cytokines_cube(:,14,:));
subplot(4,4,15); histogram(cytokines_cube(:,15,:));
subplot(4,4,16); histogram(cytokines_cube(:,16,:));

%%
% Center and scale
% Centering can't hurt in this case, scaling I'm not sure about
%[cytokines_cnt, cytokines_means] = centerData(cytokines_cube, 1);
%[cytokines_cnt_scl, cytokines_stds] = scaleData(cytokines_cnt, 2);
cytokines_cnt_scl = cytokines_cube;

%% Check distributions per feature
subplot(4,4,1); histogram(cytokines_cnt_scl(:,1,:));
subplot(4,4,2); histogram(cytokines_cnt_scl(:,2,:));
subplot(4,4,3); histogram(cytokines_cnt_scl(:,3,:));
subplot(4,4,4); histogram(cytokines_cnt_scl(:,4,:));
subplot(4,4,5); histogram(cytokines_cnt_scl(:,5,:));
subplot(4,4,6); histogram(cytokines_cnt_scl(:,6,:));
subplot(4,4,7); histogram(cytokines_cnt_scl(:,7,:));
subplot(4,4,8); histogram(cytokines_cnt_scl(:,8,:));
subplot(4,4,9); histogram(cytokines_cnt_scl(:,9,:));
subplot(4,4,10); histogram(cytokines_cnt_scl(:,10,:));
subplot(4,4,11); histogram(cytokines_cnt_scl(:,11,:));
subplot(4,4,12); histogram(cytokines_cnt_scl(:,12,:));
subplot(4,4,13); histogram(cytokines_cnt_scl(:,13,:));
subplot(4,4,14); histogram(cytokines_cnt_scl(:,14,:));
subplot(4,4,15); histogram(cytokines_cnt_scl(:,15,:));
subplot(4,4,16); histogram(cytokines_cnt_scl(:,16,:));

%%
% Might want to add editing the metadata here
[I,J,K] = size(cytokines_cnt_scl);

subjectMeta = cytokines_meta(:,[1 6]);
subjectMeta = sortrows(subjectMeta, 1);
subjectMeta = unique(subjectMeta, "rows", "stable"); % stable stops resorting rows

subjectMeta(subjectMeta(:,2) == "A", 2) = "Asymptomatic";
subjectMeta(subjectMeta(:,2) == "S", 2) = "Symptomatic";

%featureMeta = [1:J; ones(1,J,1)]';
featureMeta = ["VEGF" "CRP" "GM-CSF" "IL1alpha" "IL1beta" "IL4" "IL6" "IL8" "IL10" "IL12p70" "IL17A" "IFNgamma" "MIP1alpha" "OPG" "TNFalpha" "RANKL"]';

%%
% Initialize options
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

timeMeta =  ["Before extraction" "After extraction" "After extraction", "After extraction"]';
timeMeta_filtered = timeMeta;

% Remove outlier subjects
%df_cnt_scl_filtered(:,4,:) = [];
%featureMeta_filtered(4,:) = [];

% This subject really is borked
%df_cnt_scl_filtered(25,:,:) = [];
%subjectMeta_filtered(25,:) = [];

% Remove timepoint 6
%df_cnt_scl_filtered(:,:,6) = [];

%%
% Run PARAFACs
path_start = "./test_run/Figures/";
maxComponents=4;
numReps=25;
maxIterations=20;
metaData = {subjectMeta_filtered, featureMeta_filtered, timeMeta_filtered};
legendIndex = [2 1 1];
numColsPerLegend = [2 3 2];
resort = [true false false];
xlabels = ["Subject index", "Feature index", "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
balancedJackKnifing = false;
subjectGroupCol = 2;

[Models, Cons, VarExps, Boots, BootVarExps, Tuckers] = quickReport(df_cnt_scl_filtered, Options, const, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, metaData, resort, legendIndex, numColsPerLegend, xlabels, titles, "PARAFAC Athina", path_start+"Athina_");

%%
% Pick the number of components that seems appropriate based on the plots
% and pick the model that you like the best.
% Best practices: pick between best Corcondia or best variance explained.
% Alternative: choose a good tucker congruence model.

numFactors = 2;
choice = find(VarExps{numFactors}==max(VarExps{numFactors}));
Model = pickModel(Models, numFactors, choice);

%%
% Save the models
model_path = "./test_run/PARAFAC models/";

annotatedModel = annotateModel(df_cnt_scl_filtered, Model, metaData);
savePARAFAC(df_cnt_scl_filtered, Model, annotatedModel, model_path + "Athina");

%%
% Plot PARAFAC model
legendIndex = [4 3 3];
numColsPerLegend = [2 3 2];
xlabels = ["Subject index", "Feature index", "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
resort = [true false false];
path_start = "./test_run/Figures/";

plotPARAFAC4(annotatedModel, numFactors, VarExps, choice, resort, legendIndex, numColsPerLegend, xlabels, titles, "PARAFAC Athina", path_start + "PARAFAC_Athina.jpg");
