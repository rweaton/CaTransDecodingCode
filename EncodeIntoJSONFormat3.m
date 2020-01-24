function EncodeIntoJSONFormat3(varargin)

    if (nargin == 0)
        
        % Obtain file to analyze using the gui file browser
        [Name, Filepath] = uigetfile;
        
        % Assemble full path to source file
        SourceFile = fullfile(Filepath, Name);
        
        % Unpack file parts to construct save name.
        [junk, Name, Ext] = fileparts(Name);

        % Load matlab variables into current workspace        
        load(SourceFile);
        
    end
        
    if (nargin == 1)
        
        SourceFile = varargin(1);
        
        % Extract path and file parts from SourceFile argument
        [Filepath, Name, Ext] = fileparts(SourceFile);
        
        % Load matlab variables into current workspace
        load(SourceFile);
        
    end    

    % Construct name of json calcium file to be saved.
    SaveName = fullfile(Filepath, [Name, '_new_unique_C.json']);

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
    SaveName = fullfile(Filepath, [Name, '_new_unique_B.json']);
    
    % Save parsed json data to a file in same directory as mat file.
    fid = fopen(SaveName, 'w');
    fwrite(fid, B_data, 'char');
    fclose(fid);
    %%%% End Calcium Data processing %%%%
    
    % Clear variables from workspace
    clear('B','C')

end