function saveStringArray(filepath,strArray)
% Write string array strArray to a text file specified by filepath

fid = fopen(filepath,'w');
fmtString = [repmat('%s\t',1,size(strArray,2)-1),'%s\n'];
fprintf(fid,fmtString,strArray{:});
fclose(fid);
