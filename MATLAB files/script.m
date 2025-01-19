% PARAFAC functionality
addpath("../N-way-shell/Scripts/"); % own scripts
addpath("../N-way-shell/N-way toolbox/"); % from Rasmus Bro

%%
% Load raw cytokines data
cytokines_data = readmatrix(".\input_deduplicated_RvdP.csv", Filetype="delimitedtext", Delimiter=" ");
cytokines_meta = readmatrix(".\input_deduplicated_metadata_RvdP.csv", Filetype="delimitedtext", Delimiter=" ", OutputType="string");

%%
% Attach pain/no-pain metadata to these case subjects
meta2 = readmatrix("./Data/Root_meta_data_parafac.txt", Filetype="delimitedtext", Delimiter="\t", OutputType="string");
meta2 = meta2(:,[2,8]);

cytokines_meta(:,7) = "C";
allIndividuals = unique(cytokines_meta(:,1));

for i=1:length(allIndividuals)
    subjectID = allIndividuals(i);
    if ismember(subjectID, meta2(:,1))
        status = meta2(meta2(:,1) == subjectID, 2);
        if length(status) > 1 % somehow there are duplicate entries
            status = status(1);
        end
        cytokines_meta(cytokines_meta(:,1) == subjectID, 7) = status;
    end
end

%%
% Log transform
[I,J] = size(cytokines_data);
vectorizedX = reshape(cytokines_data, I*J, 1);
pseudocount = min(vectorizedX(vectorizedX>0));

cytokines_data_log = log(cytokines_data + pseudocount);

%%
% Reshape to 3-way matrix
keepIndividuals = true;
cytokines_cube = rawDataToCube(cytokines_data_log, cytokines_meta(:,1), str2double(cytokines_meta(:,2)), keepIndividuals);

%%
% Center and scale
[cytokines_cnt, cytokines_means] = centerData(cytokines_cube, 1);
[cytokines_cnt_scl, cytokines_stds] = scaleData(cytokines_cnt, 2);

%%
% Might want to add editing the metadata here
[I,J,K] = size(cytokines_cnt_scl);

subjectMeta = cytokines_meta(:,[1 7]);
subjectMeta = sortrows(subjectMeta, 1);
subjectMeta = unique(subjectMeta, "rows", "stable"); % stable stops resorting rows

subjectMeta(subjectMeta(:,2) == "A", 2) = "Asymptomatic";
subjectMeta(subjectMeta(:,2) == "S", 2) = "Symptomatic";
subjectMeta(subjectMeta(:,2) == "C", 2) = "Control";

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

timeMeta = ["Before extraction" "Before extraction" "Before extraction" "After extraction" "After extraction" "After extraction"]';
timeMeta_filtered = timeMeta;

% Remove outlier subjects
df_cnt_scl_filtered([25 51 52],:,:) = [];
subjectMeta_filtered([25 51 52],:) = [];

% Remove timepoint 6
%df_cnt_scl_filtered(:,:,6) = [];
%timeMeta_filtered(6) = [];

%%
% Run PARAFACs
path_start = "./test_run/Figures/";
maxComponents=5;
numReps=25;
maxIterations=20;

metaData = {subjectMeta_filtered, featureMeta_filtered, timeMeta_filtered};
legendIndex = [2 1 1];
numColsPerLegend = [3 3 2];
xlabels = ["Subject index", "Feature index", "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
resort = [true false false];

[Models, Cons, VarExps, Boots, BootVarExps, Tuckers] = quickReport(df_cnt_scl_filtered, Options, const, maxComponents, numReps, maxIterations, metaData, resort, legendIndex, numColsPerLegend, xlabels, titles, "Athina study bootstrapped", path_start+"Athina_cytokines");

%%
% Pick the number of components that seems appropriate based on the plots
% and pick the model that you like the best.
% Best practices: pick between best Corcondia or best variance explained.
% Alternative: choose a good tucker congruence model.

numFactors = 3;
choice = find(VarExps{numFactors}==max(VarExps{numFactors}));
Model = pickModel(Models, numFactors, choice);

%%
% Save the models
model_path = "./test_run/PARAFAC models/";

annotatedModel = annotateModel(df_cnt_scl_filtered, Model, metaData);
savePARAFAC(df_cnt_scl_filtered, Model, annotatedModel, model_path + "Athina");

%%
% Plot PARAFAC model
legendIndex = [5 4 4];
numColsPerLegend = [3 3 2];
xlabels = ["Subject index", "Feature index", "Time index"];
titles = ["Subject mode", "Feature mode", "Time mode"];
resort = [true false false];
path_start = "./test_run/Figures/";

plotPARAFAC4(annotatedModel, numFactors, VarExps, choice, resort, legendIndex, numColsPerLegend, xlabels, titles, "PARAFAC Athina", path_start + "PARAFAC_Athina.jpg");
