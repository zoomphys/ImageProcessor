function replace_in_file(filepath_in, filepath_out, repl_array)
% Replace strings specified by repl_array in source file and copy to
% destination folder
fid  = fopen(filepath_in,'r');
f=fread(fid,'*char')';
fclose(fid);

for iPair = 1:size(repl_array,1)
    f = strrep(f,repl_array{iPair,1},repl_array{iPair,2});
end

fid  = fopen(filepath_out,'w');
fprintf(fid,'%s',f);
fclose(fid);
