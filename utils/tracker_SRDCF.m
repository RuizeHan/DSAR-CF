% Set the path to your trackers
path_to_your_trackers = '';

if isempty(path_to_your_trackers)
    error('Set the path to your trackers!');
end

tracker_label = 'SRDCF';
tracker_command = generate_matlab_command('vot_wrapper(''SRDCF'', ''SRDCF_VOT_settings'')', {[path_to_your_trackers '/SRDCF/']});
tracker_interpreter = 'matlab';
tracker_trax = false;