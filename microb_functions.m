% PARAFAC functionality
home = ".";
cd(home)
addpath("..\N-way-shell\Scripts\"); % own scripts
addpath("..\N-way-shell\N-way toolbox/"); % from Rasmus Bro

%% Load preprocessed data from R pipeline
path = "./Tax4fun2_output/";
suffix = "_functional_analysis_raw";

% Tongue
microb = readmatrix(path+"Athina" + suffix + ".csv");
microb_meta = readmatrix(path+"Athina" + suffix + "_subjects.csv", OutputType="string");
microb_feature_meta = readmatrix(path+"Athina" + suffix + "_features.csv", OutputType="string");

%% Remove problematic duplicate samples
remove = ["A11-10 17" "A11-15 25" "A11-8 36"];
mask = ~ismember(microb_meta(:,3), remove);

microb = microb(mask,:);
microb_meta = microb_meta(mask,:);

%% Make selection masks based on sparsity
sparsityThreshold = 0.5;

sparsity = sum(microb==0) / size(microb,1);
%histogram(sparsity)

sparsity_selection = sparsity <= sparsityThreshold;

%% Make selection masks based on variation
variationThreshold = 1e-5;

ssq = sum(microb.^2);
%histogram(ssq);

variation_selection = ssq >= variationThreshold;

%% CLR
microb_clr = transformCLR(microb);

%% Filter the data based on the selections
microb_clr_filtered = microb_clr(:, sparsity_selection & variation_selection);

%% Center and scale
microbiome_cnt = microb_clr_filtered - mean(microb_clr_filtered);
microbiome_cnt_scl = microbiome_cnt ./ std(microbiome_cnt);

%% Feature metadata
microb_feature_meta_filtered = microb_feature_meta(sparsity_selection & variation_selection,:);