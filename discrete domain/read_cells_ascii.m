
function [n, connected_nodes] = read_cells_ascii(filepath)

hdr = {};
out_cellarray = {}; % store the information in a cell array first
[fid, msg] = fopen(filepath, 'rt'); % open the file for reading in text mode
assert(fid >= 3, msg) % ensure the file opened correctly.
fgetl(fid); % read and ignore the very first line.
while ~feof(fid)
    hdr{end+1} = fgetl(fid);
    temp = textscan(hdr{end}, '%f %*s %f %f %f %f', 'CollectOutput', true); % see textscan() documentation
    out_cellarray{end+1} = transpose(temp{:,:});
end
fclose_status = fclose(fid); % fclose_status = 0 if operation is successful, otherwise = -1
out_mat = transpose(cell2mat(out_cellarray)); % convert the content in the cell array to matrix

n = out_mat(:,1);
connected_nodes = out_mat(:,2:end);
end