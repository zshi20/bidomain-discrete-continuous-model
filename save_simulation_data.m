function save_simulation_data(folder_name, varargin)
%SAVE_SIMULATION_DATA saves the variables into a new .mat file in the specified folder. The filename will 
%follow the pattern 'simulationX.mat' where X is the next index in the series.
%The first argument is the directory where the folder is being saved. 
%The second argument is the model parameter
%The third argument is the output of simulation (cell location, etc)
%The fourth argument, which is optional, is the setting of the parameter sweep.

    p = inputParser;
    addRequired(p, 'folder_name', @(x) isstring(x)||ischar(x));
    addRequired(p, 'modelParam', @isstruct);
    addRequired(p, 'outputMatrix', @isstruct);
    addOptional(p, 'swpParam', struct, @isstruct);
    parse(p, folder_name, varargin{:});

    modelParam = p.Results.modelParam;
    swpParam   = p.Results.swpParam;
    outputMatrix = p.Results.outputMatrix;

    file_list = dir(fullfile(folder_name, 'simulation*.mat'));

    max_index = 0;

    for i = 1:length(file_list)
        
        index_string = regexp(file_list(i).name, 'simulation(\d+)\.mat', 'tokens');

        if ~isempty(index_string)
            current_index = str2double(index_string{1}{1});

            if current_index > max_index
                max_index = current_index;
            end
        end
    end

    new_file_name = fullfile(folder_name, sprintf('simulation%d.mat', max_index + 1));

    save(new_file_name, 'modelParam', 'outputMatrix', 'swpParam', '-v7.3');

    fprintf('Variables saved to %s\n', new_file_name);
end
