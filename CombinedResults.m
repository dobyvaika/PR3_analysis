% Norādiet mapi, kur atrodas jūsu .csv faili
folderPath = 'C:\Users\...\XLSX formāts\results';
files = dir(fullfile(folderPath, '*.xlsx'));

dataTables = struct('postTal', [], 'postTuv', [], 'preTal', [], 'preTuv', []);

for k = 1:length(files)
    filename = fullfile(folderPath, files(k).name);
    
    if strcmp(files(k).name, 'CombinedResults.xlsx')
        continue;
    end

    fileInfo = regexp(files(k).name, '(\d+Q)(post|pre)(Tal|Tuv)Dalibnieks[-_]?(\d+)?-\d{8}-\d{6}(_\d+)?_res\.xlsx', 'tokens');
    if isempty(fileInfo)
        warning('No matches found in the file name: %s', files(k).name);
        continue;
    end
    
    fileInfo = fileInfo{1};
    
    try
        numChars = readtable(filename, 'Sheet', 'Numerical characteristics');
        fourierData = readtable(filename, 'Sheet', 'Fourier');
    catch me
        warning('Problem with the file: %s. Error: %s', files(k).name, me.message);
        continue;
    end

    numChars.Properties.VariableNames = strcat('NC_', numChars.Properties.VariableNames);
    fourierData.Properties.VariableNames = strcat('F_', fourierData.Properties.VariableNames);
    
    sheetKey = strcat(fileInfo{2}, fileInfo{3}); % Piemēram, postTal
    tempTable = table(string(fileInfo{1}), 'VariableNames', {'Participant'});
    tempTable = [tempTable, numChars, fourierData];
    
    if isfield(dataTables, sheetKey)
        dataTables.(sheetKey) = [dataTables.(sheetKey); tempTable];
    else
        dataTables.(sheetKey) = tempTable;
    end
end

outputFilename = fullfile(folderPath, 'CombinedResults.xlsx');
for fn = fieldnames(dataTables)'
    if ~isempty(dataTables.(fn{1}))
        writetable(dataTables.(fn{1}), outputFilename, 'Sheet', fn{1});
    end
end

disp('Data processing and writing to sheets completed.');
