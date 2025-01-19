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

%% Microbiome
% Load raw microbiome data
microbiome_raw = readmatrix("./Data/20240429_microbiome_counts.csv", Filetype="delimitedtext", Delimiter=" ");
taxonomy = readmatrix("./Data/20240429_taxonomy.csv", Filetype="delimitedtext", Delimiter=" ", OutputType="string");
subjectMeta2 = readmatrix("./Data/20240429_microbiome_sampleMeta.csv", Filetype="delimitedtext", Delimiter=" ", OutputType="string");

remove = ["A11-10 17" "A11-15 25" "A11-8 36"];
mask = ~ismember(subjectMeta2(:,3), remove);

microbiome_raw = microbiome_raw(mask,:);
subjectMeta2 = subjectMeta2(mask,:);

% Select ASVs based on sparsity
sparsityThreshold = 0.5;

sparsity = sum(microbiome_raw==0) / size(microbiome_raw,1);
featureMask = sparsity < sparsityThreshold;

% CLR transformation
microbiome_clr = transformCLR(microbiome_raw);

% Reduce to previously selected ASVs
microbiome_selected = microbiome_clr(:,featureMask);
taxonomy_selected = taxonomy(featureMask,:);

% Center and scale
microbiome_cnt = microbiome_selected - mean(microbiome_selected);
microbiome_cnt_scl = microbiome_cnt ./ std(microbiome_cnt);

%% Homogenize subjects
homogenized_subjectMeta = subjectMeta;

cytokine_cnt_scl_filtered = cytokine_cnt_scl_filtered;
microbiome_cnt_scl_filtered = microbiome_cnt_scl;

%% Remove outliers
%cytokine_cnt_scl_filtered(25,:,:) = [];
%microbiome_cnt_scl_filtered(25,:,:) = [];
%homogenized_subjectMeta(25,:) = [];

%% Create X
sizes = [26 16 6 76];
modes = {[1 2 3], [1 4]};
lambdas = {[1 1 1], [1 1 1]};

X = cytokine_cnt_scl_filtered;
Y = microbiome_cnt_scl_filtered;
Z = setupCMTFdata(X, Y, sizes, modes, true);

%% Set options
options = ncg('defaults');
options.Display ='final';
options.MaxFuncEvals = 100000;
options.MaxIters     = 10000;
options.StopTol      = 1e-8;
options.RelFuncTol   = 1e-8;
%beta = 0.00001; for ACMTF

maxComponents = 3;
numReps = 25;
balancedJackKnifing = true;
subjectGroupCol = 2;
metaData = {homogenized_subjectMeta, featureMeta_filtered, timeMeta_filtered, taxonomy_selected};
resort = [true false false true];
legendIndex = [2 1 1 2];
numColsPerLegend = [2 2 2 2];
xlabels = ["Subject index", "Feature index", "Time index", "Feature index"];
titles = ["A", "B", "C", "D"];
overallTitle = "CMTF";
path_start="./test_run_ACMTF/";
maxIterations = 26;

%% ACMTF
% fit ACMTF-OPT
R = 3;
P = 2;
beta = [1e-3 1e-3];
% need norms

[Zhat,G,out]    = acmtf_opt(Z,R,'alg_options',options,'beta',beta);        
%data.Zhat  = Zhat;
%data.W     = W;
%data.Xorig = X;
%data.Init  = Init;
%data.out   = out;
%data.Atrue = Atrue;
l_rec = zeros(P, R);
for p = 1:P
    temp        = normalize(Zhat{p});
    l_rec(p,:)  = temp.lambda;    
end
tt = [];
for i  = 1:length(lambdas)
    tt = [tt; lambdas{i}];
end
data.lambdas       = tt;
data.lambda_rec    = l_rec;
data.adjlambda_rec =  khatrirao(norms,l_rec');
data.norms         = norms;

%% Quick Report
[allModels,allVarExpsX,allVarExpsY, bootstrappedModels, bootVarExpsX, bootVarExpsY, FMS_tensor_result, FMS_matrix_result, goodness_X] = ACMTF_quickReport(X, Y, sizes, options, maxComponents, numReps, maxIterations, balancedJackKnifing, subjectGroupCol, metaData, resort, legendIndex, numColsPerLegend, xlabels, titles, overallTitle, path_start+"Athina");

%% Pick your model
numComponents = 1;
choice = 1;

%% Fit CMTF_OPT using one of the first-order optimization algorithms 
Fac = {allModels{1,numComponents}(:,:,choice) allModels{2,numComponents}(:,:,choice) allModels{3,numComponents}(:,:,choice) allModels{4,numComponents}(:,:,choice)};
subplot(1,4,1); bar(Fac{1});
subplot(1,4,2); bar(Fac{2});
subplot(1,4,3); bar(Fac{3});
subplot(1,4,4); bar(Fac{4});

%% Save the model
path = "./test_run_ACMTF/";

writematrix(Fac{1}, path+"Athina_ACMTF_A.csv");
writematrix(Fac{2}, path+"Athina_ACMTF_B.csv");
writematrix(Fac{3}, path+"Athina_ACMTF_C.csv");
writematrix(Fac{4}, path+"Athina_ACMTF_D.csv");
writematrix(homogenized_subjectMeta, path+"Athina_homogenized_subjectMeta.csv");
writematrix(featureMeta_filtered, path+"Athina_cytokines.csv");
writematrix(timeMeta_filtered, path+"Athina_time.csv");
writematrix(taxonomy_selected, path+"Athina_taxonomy.csv");