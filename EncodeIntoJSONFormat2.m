function EncodeIntoJSONFormat2(SourceFile)
    
    % Load matlab variables into current workspace
    load(SourceFile);

    % Extract path and file parts from SourceFile argument
    [Filepath, Name, Ext] = fileparts(SourceFile);

    % Construct name of json calcium file to be saved.
    SaveName = fullfile(Filepath, [Name, '_C.json']);

    %%%% Begin Calcium Data processing %%%%
    % Parse data into json format
    C_data = jsonencode(C);

    % Save parsed json data to a file in same directory as mat file.
    fid = fopen(SaveName, 'w');
    fwrite(fid, C_data, 'char');
    fclose(fid);
    %%%% End Calcium Data processing %%%%
    
    %%%% Begin Behavioral Data processing %%%%
    % Write table contents into struct as required by the json encoder
    Behavior_struct = struct;
    column_names = B.Properties.VariableNames;
    for i = 1:length(column_names)
        column = column_names{i};
        Behavior_struct.(column) = B.(column);
    end

    % Parse behavioral data into json format
    B_data = jsonencode(Behavior_struct);
    
    % Construct name of json behavioral file to be saved.
    SaveName = fullfile(Filepath, [Name, '_B.json']);
    
    % Save parsed json data to a file in same directory as mat file.
    fid = fopen(SaveName, 'w');
    fwrite(fid, B_data, 'char');
    fclose(fid);
    %%%% End Calcium Data processing %%%%
    
    % Clear variables from workspace
    clear('B','C')

end