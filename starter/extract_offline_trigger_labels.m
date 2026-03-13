%% Extract trigger labels (codes) from all offline GDF files
% Outputs:
%   - fig_data/offline_trigger_inventory.csv
%   - fig_data/offline_unique_trigger_labels.txt

clear; clc;

repo_root = resolve_repo_root_from_pwd(pwd);
addpath(genpath(fullfile(repo_root, 'starter', 'functions')));

dataset_root = resolve_dataset_root(repo_root);
offline_root = fullfile(dataset_root, 'offline_data');
files = dir(fullfile(offline_root, 'Sub_*', '*offline*', '*.gdf'));
assert(~isempty(files), 'No offline GDF files found under %s', offline_root);

[~, ord] = sort(cellfun(@(a,b) fullfile(a,b), {files.folder}, {files.name}, 'UniformOutput', false));
files = files(ord);

subject_id = zeros(numel(files),1);
run_id = zeros(numel(files),1);
file_relpath = strings(numel(files),1);
unique_trigger_codes = strings(numel(files),1);
n_unique = zeros(numel(files),1);
n_events = zeros(numel(files),1);

all_codes = [];

for i = 1:numel(files)
    fp = fullfile(files(i).folder, files(i).name);
    [~, HDR] = sload(fp);

    if isfield(HDR, 'EVENT') && isfield(HDR.EVENT, 'TYP') && ~isempty(HDR.EVENT.TYP)
        codes = unique(double(HDR.EVENT.TYP(:)));
        n_events(i) = numel(HDR.EVENT.TYP);
    else
        codes = [];
        n_events(i) = 0;
    end

    all_codes = [all_codes; codes(:)]; %#ok<AGROW>
    n_unique(i) = numel(codes);
    unique_trigger_codes(i) = strjoin(string(codes.'), ';');

    m_sub = regexp(files(i).name, 'Sub_(\d+)_run_(\d+)_', 'tokens', 'once');
    if isempty(m_sub)
        m_sub = regexp(fullfile(files(i).folder, files(i).name), 'Sub_(\d+).*_run_(\d+)_', 'tokens', 'once');
    end
    if ~isempty(m_sub)
        subject_id(i) = str2double(m_sub{1});
        run_id(i) = str2double(m_sub{2});
    else
        subject_id(i) = NaN;
        run_id(i) = NaN;
    end

    file_relpath(i) = string(strrep(fp, [repo_root filesep], ''));
end

T = table(subject_id, run_id, file_relpath, n_events, n_unique, unique_trigger_codes);
T = sortrows(T, {'subject_id', 'run_id'});

out_dir = fullfile(repo_root, 'fig_data');
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end

csv_out = fullfile(out_dir, 'offline_trigger_inventory.csv');
writetable(T, csv_out);

global_unique = unique(all_codes);
txt_out = fullfile(out_dir, 'offline_unique_trigger_labels.txt');
fid = fopen(txt_out, 'w');
fprintf(fid, 'Unique trigger codes across all offline files (decimal):\n');
fprintf(fid, '%s\n', strjoin(string(global_unique.'), ', '));
fprintf(fid, '\nCount: %d\n', numel(global_unique));
fclose(fid);

fprintf('Wrote:\n  %s\n  %s\n', csv_out, txt_out);
fprintf('Global unique trigger codes:\n  %s\n', strjoin(string(global_unique.'), ', '));

%% Local helpers
function repo_root = resolve_repo_root_from_pwd(start_dir)
repo_root = start_dir;
for k = 1:10
    has_starter = exist(fullfile(repo_root, 'starter', 'functions'), 'dir') == 7;
    has_data = exist(fullfile(repo_root, 'BCI_Harmony_ExperimentalData'), 'dir') == 7 || ...
               exist(fullfile(repo_root, 'BCI_course_EXP'), 'dir') == 7;
    if has_starter && has_data
        return;
    end
    parent = fileparts(repo_root);
    if strcmp(parent, repo_root)
        break;
    end
    repo_root = parent;
end
error('Could not locate repository root from current directory.');
end

function dataset_root = resolve_dataset_root(repo_root)
new_name = fullfile(repo_root, 'BCI_Harmony_ExperimentalData');
old_name = fullfile(repo_root, 'BCI_course_EXP');
if exist(new_name, 'dir') == 7
    dataset_root = new_name;
elseif exist(old_name, 'dir') == 7
    dataset_root = old_name;
else
    error('Could not find dataset folder under %s', repo_root);
end
end
