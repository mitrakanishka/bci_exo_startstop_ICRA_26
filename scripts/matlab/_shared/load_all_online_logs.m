function T = load_all_online_logs(log_root)
% Parse all online Python logs into a tidy trial-level table.

files = dir(fullfile(log_root, 'Sub_*', '*.txt'));
assert(~isempty(files), 'No online python log files found under %s', log_root);

parts = cell(numel(files), 1);
for i = 1:numel(files)
    parts{i} = parse_online_log_file(fullfile(files(i).folder, files(i).name));
end
parts = parts(~cellfun(@isempty, parts));
assert(~isempty(parts), 'No trials parsed from online logs.');

T = vertcat(parts{:});
T = sortrows(T, {'subject', 'session', 'run', 'trial'});
end

function T = parse_online_log_file(path)
lines = readlines(path);
name_pat = 'sub_(\d+)_session_(\d+)_run_(\d+)_online_python_log\.txt$';
head_pat = 'Subject\s+(\d+)\s*,\s*Session\s+(\d+)\s*,\s*Run\s+(\d+)';
trial_pat = ['Trial:\s*(-?\d+).*?' ...
    'Onset\s+Success:\s*(-?\d+).*?' ...
    'Onset\s+Decode\s*time:\s*(-?\d*\.?\d+)s.*?' ...
    'Offset\s+Success:\s*(-?\d+).*?' ...
    'Offset\s+Decode\s*time:\s*(-?\d*\.?\d+)s'];

sub_file = []; ses_file = []; run_file = [];
m_name = regexp(path, name_pat, 'tokens', 'once', 'ignorecase');
if ~isempty(m_name)
    sub_file = str2double(m_name{1});
    ses_file = str2double(m_name{2});
    run_file = str2double(m_name{3});
end

sub_head = []; ses_head = []; run_head = [];
for i = 1:min(3, numel(lines))
    m_head = regexp(lines(i), head_pat, 'tokens', 'once', 'ignorecase');
    if ~isempty(m_head)
        sub_head = str2double(m_head{1});
        ses_head = str2double(m_head{2});
        run_head = str2double(m_head{3});
        break;
    end
end

subject = choose_first(sub_file, sub_head);
session = choose_first(ses_file, ses_head);
run_id = choose_first(run_file, run_head);
assert(~isempty(subject) && ~isempty(session) && ~isempty(run_id), ...
    'Could not parse subject/session/run from %s', path);

rows = [];
for i = 1:numel(lines)
    m_trial = regexp(lines(i), trial_pat, 'tokens', 'once', 'ignorecase');
    if isempty(m_trial)
        continue;
    end
    onset_code = str2double(m_trial{2});
    offset_code = str2double(m_trial{4});
    row = [subject, session, run_id, str2double(m_trial{1}), ...
        onset_code, clean_rt(m_trial{3}), offset_code, clean_rt(m_trial{5}), double(offset_code ~= -1)];
    rows = [rows; row]; %#ok<AGROW>
end

if isempty(rows)
    T = table();
    return;
end

n = size(rows, 1);
file_col = repmat(string(path), n, 1);
onset_label = arrayfun(@(x) decode_label(x, 'onset'), rows(:, 5));
offset_label = arrayfun(@(x) decode_label(x, 'offset'), rows(:, 7));
T = table( ...
    rows(:, 1), rows(:, 2), rows(:, 3), rows(:, 4), ...
    rows(:, 5), onset_label, rows(:, 6), ...
    rows(:, 7), offset_label, rows(:, 8), logical(rows(:, 9)), file_col, ...
    'VariableNames', {'subject', 'session', 'run', 'trial', ...
    'onset_code', 'onset_label', 'onset_rt_s', ...
    'offset_code', 'offset_label', 'offset_rt_s', 'offset_attempted', 'file'});
end

function value = choose_first(a, b)
if ~isempty(a)
    value = a;
else
    value = b;
end
end

function label = decode_label(code, phase)
if code == 1
    label = "hit";
elseif code == 0
    label = "miss";
elseif code == -2
    label = "timeout";
elseif strcmpi(phase, 'offset') && code == -1
    label = "not_attempted";
else
    label = "unknown";
end
end

function value = clean_rt(token)
value = str2double(token);
if value < 0
    value = NaN;
end
end
