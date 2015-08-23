function name = findLastInDir(filedir,isDir)
files = dir(filedir);
if length(files)<1
    name = '';
    return;
end

if isDir
    valid_names= find([files.isdir]);
else
    valid_names= find(~[files.isdir]);
end

file_date = [files.datenum];
[~, idx]=sort(file_date);

%keep index of file names only
idx = idx(ismember(idx,valid_names)); 

name = files(idx(end)).name;