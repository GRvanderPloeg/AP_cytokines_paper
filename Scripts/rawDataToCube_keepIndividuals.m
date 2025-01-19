function Xcube=rawDataToCube_keepIndividuals(Xlong, Xmeta_subjects, Xmeta_timepoints, numTimepoints)

% Create metadata overview
Xmeta = [Xmeta_subjects Xmeta_timepoints];

% Identify how many samples are missing
subjects = unique(Xmeta_subjects);
expectedNumSamples = length(subjects) * numTimepoints;
realNumSamples = size(Xlong,1);

numMissingSamples = expectedNumSamples - realNumSamples;

% Identify which samples are missing
missingSamples = zeros(numMissingSamples, 2);
missingSamples = string(missingSamples);
numOccurences = groupcounts(Xmeta(:,1));
missingSampleIterator = 1;

% Identification of individuals missing one or more samples
for i=1:length(subjects)
    subject = subjects(i);
    discrepancy = numTimepoints - numOccurences(i);
    if (discrepancy > 0)
        for j=1:discrepancy
            missingSamples(missingSampleIterator, 1) = subject;
            missingSampleIterator = missingSampleIterator + 1;
        end
    end 
end

% Identification of the timepoints that are missing
missingSubjects = unique(missingSamples(:,1));
missingSampleIterator = 1;

for i=1:length(missingSubjects)
    subject = missingSubjects(i,1);
    timepoints = str2double(Xmeta(Xmeta(:,1) == subject,2));
    expectedTimepoints = 1:6;

    ii = ~ismember(expectedTimepoints,timepoints);
    values = expectedTimepoints(ii);

    if (length(values) > 1)
        for j=1:length(values)
            missingSamples(missingSampleIterator,2) = string(values(j));
            missingSampleIterator = missingSampleIterator + 1;
        end
    else
        missingSamples(missingSampleIterator,2) = string(values);
        missingSampleIterator = missingSampleIterator + 1;
    end
end

% Append them to the dataset
missingData = nan(numMissingSamples, size(Xlong,2));
Xlong = [Xlong; missingData];

% Sort the data
Xmeta = [Xmeta; missingSamples];
Xmeta(:,end+1) = 1:size(Xmeta,1);
Xmeta = sortrows(Xmeta, [2 1]); % sort on visit number (ascending), then subject name (alphabetical)
keepRowIndices = str2double(Xmeta(:,end)); % new sorting of data
Xlong = Xlong(keepRowIndices, :);

% Reshape into cube
I = size(Xlong, 1) / numTimepoints;
J = size(Xlong, 2);
K = numTimepoints;

Xcube = reshape(Xlong, I, K, J);
Xcube = permute(Xcube, [1 3 2]);

