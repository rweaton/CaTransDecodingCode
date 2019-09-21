C_data = jsonencode(C);
dir_path = uigetdir();
[filepath, name, ext] = fileparts(dir_path);
filename = fullfile(dir_path, 'test_C.json');
fid = fopen(filename, 'w');
fwrite(fid, C_data, 'char');
fclose(fid)

test_struct = struct;
column_names = B.Properties.VariableNames;
for i = 1:length(column_names)
    column = column_names{i};
    test_struct.(column) = B.(column);
end

B_data = jsonencode(test_struct);
dir_path = uigetdir();
filename = fullfile(dir_path, 'test_B.json');
fid = fopen(filename, 'w');
fwrite(fid, B_data, 'char');
fclose(fid)