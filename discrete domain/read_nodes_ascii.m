
function [n, xyz_coord] = read_nodes_ascii(filepath, delimiter)

hdr = {};
out_cellarray = {}; % store the information in a cell array first
[fid, msg] = fopen(filepath, 'rt'); % open the file for reading in text mode
assert(fid >= 3, msg) % ensure the file opened correctly.
fgetl(fid); % read and ignore the very first line.
while ~feof(fid)
    hdr{end+1} = fgetl(fid);
    out_cellarray(end+1) = textscan(hdr{end}, '%f', 'CollectOutput', true, 'Delimiter', delimiter); % see textscan() documentation
end
fclose_status = fclose(fid); % fclose_status = 0 if operation is successful, otherwise = -1
out_mat = transpose(cell2mat(out_cellarray)); % convert the content in the cell array to matrix

n = out_mat(:,1);
xyz_coord = out_mat(:,2:end);
end