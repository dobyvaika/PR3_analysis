% Norādiet mapi, kur atrodas jūsu .csv faili
sourceFolder = 'C:\Users\...\Original'; 
destinationFolder = 'C:\Users\...\XLSX formāts'; 

if ~exist(destinationFolder, 'dir')
    mkdir(destinationFolder);
end
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

csvFiles = dir(fullfile(sourceFolder, '*.csv'));


for i = 1:length(csvFiles)
    csvFileName = fullfile(sourceFolder, csvFiles(i).name);
       dataTable = readtable(csvFileName, 'PreserveVariableNames', true);
    
     [~, name, ~] = fileparts(csvFiles(i).name);
    xlsxFileName = fullfile(destinationFolder, [name, '.xlsx']);
    
   writetable(dataTable, xlsxFileName, 'FileType', 'spreadsheet');
end

disp('Konvertācija pabeigta.');