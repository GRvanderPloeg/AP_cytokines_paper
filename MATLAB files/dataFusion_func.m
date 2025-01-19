% ACMTF
addpath("../CMTF_shell/tensor_toolbox-v3.1/");
addpath("../CMTF_shell/poblano_toolbox/")
addpath("../CMTF_shell/CMTF_Toolbox/");
addpath("../CMTF_shell/L-BFGS-B-C/");
addpath("../CMTF_shell/Scripts/");
addpath("../N-way-shell/Scripts/"); % own scripts
addpath("../N-way-shell/N-way toolbox/"); % by Rasmus Bro

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

cytokines_meta(:,7) = "s";
allIndividuals = unique(cytokines_meta(:,1));

for i=1:length(allIndividuals)
    subjectID = allIndividuals(i);
    status = meta2(meta2(:,1) == subjectID, 2);
    if length(status) > 1 % somehow there are duplicate entries
        status = status(1);
    end
    cytokines_meta(cytokines_meta(:,1) == subjectID, 7) = status;
end

%% Impute missing data?
%cytokines_data(isnan(cytokines_data)) = 0; % set to 0 for now

%% Remove low library size subjects here too
remove = ["A11-18" "A11-3"];
mask = ~ismember(cytokines_meta(:,1), remove);

cytokines_data = cytokines_data(mask,:);
cytokines_meta = cytokines_meta(mask,:);

%%
% Log transform
[I,J] = size(cytokines_data);
vectorizedX = reshape(cytokines_data, I*J, 1);
pseudocount = min(vectorizedX(vectorizedX>0));

cytokines_data_log = log(cytokines_data + pseudocount);

%%
% Reshape to 3-way matrix
keepIndividuals = true;
[cytokines_cube, cytokine_subjectMeta, cytokine_conditionMeta] = rawDataToCube(cytokines_data_log, cytokines_meta(:,1), cytokines_meta(:,2), keepIndividuals);

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

%featureMeta = [1:J; ones(1,J,1)]';
featureMeta = ["VEGF" "CRP" "GM-CSF" "IL1alpha" "IL1beta" "IL4" "IL6" "IL8" "IL10" "IL12p70" "IL17A" "IFNgamma" "MIP1alpha" "OPG" "TNFalpha" "RANKL"]';

%%
% Filter out outliers
cytokine_cnt_scl_filtered = cytokines_cnt_scl;
subjectMeta_filtered = subjectMeta;
featureMeta_filtered = featureMeta;

timeMeta =  ["Before extraction" "Before extraction" "Before extraction" "After extraction" "After extraction", "After extraction"]';
timeMeta_filtered = timeMeta;

%% Functionally annotated microbiome

% Load preprocessed data from R pipeline
path = "./20240724_Tax4fun2_output/";
suffix = "_functional_analysis_raw";

% Tongue
microb = readmatrix(path+"Athina" + suffix + ".csv");
microb_meta = readmatrix(path+"Athina" + suffix + "_subjects.csv", OutputType="string");
microb_feature_meta = readmatrix(path+"Athina" + suffix + "_features.csv", OutputType="string");

% Remove problematic duplicate samples
remove = ["A11-18" "A11-3" "A11-8 36" "A11-10 17" "A11-15 17"];
mask = ~ismember(microb_meta(:,3), remove);

microb = microb(mask,:);
microb_meta = microb_meta(mask,:);

% % Make selection masks based on sparsity
% sparsityThreshold = 0.5;
% 
% sparsity = sum(microb==0) / size(microb,1);
% %histogram(sparsity)
% 
% sparsity_selection = sparsity <= sparsityThreshold;
% 
% % Make selection masks based on variation
% variationThreshold = 5e-6 %1e-5; %5e-6;
% 
% ssq = sum(microb.^2);
% %histogram(ssq);
% 
% variation_selection = ssq >= variationThreshold;

% CLR
microb_clr = transformCLR(microb);

% Filter the data based on the selections
microb_clr_filtered = microb_clr(:,:); %sparsity_selection & variation_selection);

% Center and scale
microbiome_cnt = microb_clr_filtered - mean(microb_clr_filtered);
microbiome_cnt_scl = microbiome_cnt ./ std(microbiome_cnt);

% Feature metadata
microb_feature_meta_filtered = microb_feature_meta(:,:);

%% Homogenize subjects
subjectMask = ismember(subjectMeta(:,1), microb_meta(:,2));
homogenized_subjectMeta = subjectMeta(subjectMask,:);

cytokine_cnt_scl_filtered = cytokine_cnt_scl_filtered(subjectMask,:,:);
microbiome_cnt_scl_filtered = microbiome_cnt_scl;

%% Remove outliers
%cytokine_cnt_scl_filtered(25,:,:) = [];
%microbiome_cnt_scl_filtered(25,:,:) = [];
%homogenized_subjectMeta(25,:) = [];

%% Create X
sizes = [24 16 6 8121];
modes = {[1 2 3], [1 4]};
lambdas = {[1 1 1], [1 1 1]};

datasets = {cytokine_cnt_scl_filtered microbiome_cnt_scl_filtered};
Z = setupCMTFdata(datasets, sizes, modes, true);

%% Set options
options = ncg('defaults');
options.Display ='off';
options.MaxFuncEvals = 100000;
options.MaxIters     = 10000;
options.StopTol      = 1e-8;
options.RelFuncTol   = 1e-8;

maxComponents = 3;
numReps = 24;
balancedJackKnifing = true;
subjectGroupCol = 2;
metaData = {homogenized_subjectMeta, featureMeta_filtered, timeMeta_filtered, microb_feature_meta_filtered};
resort = [true false false false];
legendIndex = [2 1 1 0];
numColsPerLegend = [2 2 2 2];
xlabels = ["Subject index", "Feature index", "Time index", "Feature index"];
titles = ["A", "B", "C", "D"];
overallTitle = "ACMTF";
path_start="./test_run_ACMTF_func/";
maxIterations = 24;
betas = [1e-3 1e-3];
cld_threshold = 0.5;

%% Quick Report
[allModels, allErrors, allVarExps, allLambdas, allOuts, allBootstrappedModels, allBootErrors, allBootVarExps, allBootLambdas, allBootOuts, FMS_result, allHomogenizedModels, combinations] = ACMTF_quickReport(datasets, betas, sizes, modes, options, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, cld_threshold, metaData, resort, legendIndex, numColsPerLegend, xlabels, titles, overallTitle, path_start+"Athina");

%% Pick your model
numComponents = 1;
choice = find(mean(allVarExps{1}) == max(mean(allVarExps{1})));

%% Fit CMTF_OPT using one of the first-order optimization algorithms 
Fac = {allModels{numComponents}{choice}{1} allModels{numComponents}{choice}{2} allModels{numComponents}{choice}{3} allModels{numComponents}{choice}{4}};
subplot(1,4,1); bar(Fac{1});
subplot(1,4,2); bar(Fac{2});
subplot(1,4,3); bar(Fac{3});
subplot(1,4,4); bar(Fac{4});

%% Save the model
path = path_start;

writematrix(Fac{1}, path+"Athina_ACMTF_A.csv");
writematrix(Fac{2}, path+"Athina_ACMTF_B.csv");
writematrix(Fac{3}, path+"Athina_ACMTF_C.csv");
writematrix(Fac{4}, path+"Athina_ACMTF_D.csv");
writematrix(homogenized_subjectMeta, path+"Athina_homogenized_subjectMeta.csv");
writematrix(featureMeta_filtered, path+"Athina_cytokines.csv");
writematrix(timeMeta_filtered, path+"Athina_time.csv");
writematrix(microb_feature_meta_filtered, path+"Athina_taxonomy.csv");
writematrix(allLambdas{numComponents}{choice}, path+"Athina_lambdas.csv");