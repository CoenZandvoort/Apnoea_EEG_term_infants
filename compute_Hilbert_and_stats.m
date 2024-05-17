% script computes EEG amplitude during breathing pauses and statistically
% compare them.
%
% Before running this script, it is assumed that (unique) apnoeas have been
% identified, and EEG data are time-locked to these apnoeas.
%
% see scripts GET_IBI, GET_APNOEAS, TIMELOCK_TO_APNOEAS,
% FIND_UNIQUE_APNOEAS, and CHECK_EEG_DATA.
%
% CZ, Feb-2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(fullfile(pwd, 'functions'))
addpath('~/Documents/MATLAB/toolboxes/fieldtrip-20210325/')
addpath('~/Documents/MATLAB/toolboxes/cmocean/')
addpath('~/Documents/MATLAB/toolboxes/brewermap/')


% EEG parameters
eeg_channels = {'Fp1', 'Fp2', 'C3', 'C4', 'T3', 'T4', 'O1', 'O2', 'Avg'};
freq_range = 1.5 : 29.5; % Hz

for time_window = 1

    switch time_window

        case 1
            ibi_type_label = 'ibi_apnoeas';
            t_window = 30; % sec
            t_last_ref = -15; % sec

            timelocking_type = 'start';
            post_st = -5; % sec
            post_en = 15; % sec

        case 2
            ibi_type_label = 'ibi_apnoeas';
            t_window = 30; % sec
            t_last_ref = -15; % sec

            timelocking_type = 'end';
            post_st = -15; % sec
            post_en = 10; % sec

        case 3
            ibi_type_label = 'ibi_apnoeas';
            t_window = 30; % sec
            t_last_ref = -15; % sec

            timelocking_type = 'end';
            post_st = -5; % sec
            post_en = 5; % sec

        case 4
            ibi_type_label = 'ibi_5_15_sec';
            t_window = 10; % sec
            t_last_ref = -5; % sec

            timelocking_type = 'start';
            post_st = -5; % sec
            post_en = 5; % sec

        case 5
            ibi_type_label = 'ibi_5_15_sec';
            t_window = 10; % sec
            t_last_ref = -5; % sec

            timelocking_type = 'end';
            post_st = -5; % sec
            post_en = 5; % sec

    end


    files = dir(fullfile(pwd, 'data/*.EDF'));
    for p = 1 : numel(files)

        % load time-locked EEG and all ibi data
        sub_id = files(p).name(1 : strfind(files(p).name, '.EDF') - 1);

        try
            load(fullfile('~/Documents/MATLAB/EEG_apnoea_Leuven/output/timelock', ibi_type_label, [sub_id, '.mat'])); % loads "data"
            ibi_data = load(fullfile(pwd, 'output/ibi', [sub_id, '_ibi.mat'])); % loads "ibi"
        catch
            continue
        end
        fprintf('process\t"%s"\n', sub_id)


        if isempty(data); continue; end
        for pause = 1 : numel(data)

            % check if apnoea should be considered (unique_apnoea is set to
            % 0 in case two apnoeas were detected in numerous breathing
            % modalities). 
            if strcmp(ibi_type_label, 'ibi_apnoeas')
                if ~data{pause}.info.unique_apnoea
                    continue
                end
            end


            output_file = fullfile('~/Documents/MATLAB/EEG_apnoea_Leuven/output/hilbert_amp', ibi_type_label, sprintf('hilbert_%s_%d.mat', sub_id, pause));
            if ~exist(output_file, 'file')

                % find the optimal control breathing period (i.e., the
                % period of X sec in which the maximal ibi is lowest - we
                % scan up to four minutes prior to the apnoea).
                % first, find timing of apnoea
                ibi = ibi_data.ibi.(strcat('ibi_', data{pause}.info.channel{1}));
                ibi_time = ibi_data.ibi.(strcat('ibi_time_', data{pause}.info.channel{1}));
                idx_t_apnoea = find(ibi_time == data{pause}.info.timing);

                t_starts = ceil(data{pause}.time{1}(1)) + 1 : t_last_ref - t_window;

                [t_start, t_stop, max_ibi] = find_optimal_ibi_window(ibi, ibi_time, idx_t_apnoea, t_starts, t_window);


                % check if the pre-window before the pause also contains an
                % inter-breath interval of 5 sec or longer
                if time_window == 4 || time_window == 5
                    idx_ibi = find(ibi_time >= data{pause}.info.timing - 60 & ibi_time < data{pause}.info.timing + 90);
                    if sum(ibi(idx_ibi) > 5) > 2; continue; end
                end


                % find EEG channels
                idx_eeg_chan = NaN(numel(eeg_channels) - 1, 1);
                for e = 1 : numel(idx_eeg_chan)

                    if any(contains(data{pause}.labels, eeg_channels{e}) & ~contains(data{pause}.labels, 'Sa') & ~contains(data{pause}.labels, 'Sp'))
                        idx_eeg_chan(e, 1) = find(contains(data{pause}.labels, eeg_channels{e}) & ~contains(data{pause}.labels, 'Sa') & ~contains(data{pause}.labels, 'Sp'));
                    end

                end
                if sum(isnan(idx_eeg_chan)) == numel(eeg_channels); continue; end


                % compute Hilbert amplitude of time-locked EEG data
                amp = NaN(numel(eeg_channels) - 1, numel(freq_range), numel(data{pause}.data{1}));
                for e = 1 : numel(idx_eeg_chan)

                    % channel may not be available or should be discarded
                    % based on its quality
                    if isnan(idx_eeg_chan(e)) || ~data{pause}.sig_quality(idx_eeg_chan(e)); continue; end


                    % time series
                    data{pause}.data{idx_eeg_chan(e)} = detrend(data{pause}.data{idx_eeg_chan(e)}, 'constant');


                    % Hilbert amplitude
                    for f = 1 : numel(freq_range)

                        % additional pre-processing...
                        [B, A] = butter(2, (freq_range(f) + [-1, 1]) ./ (data{pause}.fsample(idx_eeg_chan(e)) / 2));
                        d_ff = filtfilt(B, A, data{pause}.data{idx_eeg_chan(e)}');


                        % get the instantaneous amplitude
                        d_ff = abs(hilbert(d_ff));
                        amp(e, f, :) = log10(d_ff);
                    end

                end


                % convert to FieldTrip-struct
                freq.hilspctrm = [amp; mean(amp, 'omitnan')];
                freq.time = data{pause}.time{1};
                freq.label = eeg_channels(:);
                freq.dimord = 'chan_freq_time';
                freq.freq = freq_range;
                freq.elec = [];
                freq.cfg = [];
                freq.t_control = [t_start, t_stop];
                t_length = data{pause}.info.length; % sec


                % save output
                save(output_file, 'freq', 't_length')
                clear freq t_length amp

            end

        end

    end


    % get Hilbert spectra ready for cluster-based permutation testing
    hil_files = dir(fullfile('~/Documents/MATLAB/EEG_apnoea_Leuven/output/hilbert_amp', ibi_type_label, '*.mat')); % hil_files = dir(fullfile(pwd, 'output/hilbert_amp', ibi_type_label, '*.mat'));

    fs_eeg = 250; % Hz
    stat = cell(numel(eeg_channels), 1);
    for e = 1 : numel(eeg_channels)

        freq = cell(1, 2);

        f_count = 1;
        for h = 1 : numel(hil_files)

            % load data
            X = load(fullfile(hil_files(h).folder, hil_files(h).name));


            if ~any(any(isnan(X.freq.hilspctrm(e, :, :))))
                freq(f_count, :) = {X.freq, X.freq};


                % store amplitude of pre-window
                diff_time_window = t_window - sum(abs([post_st, post_en]));
                idx_pre = (freq{f_count, 1}.time >= X.freq.t_control(1) + diff_time_window & freq{f_count, 1}.time <= X.freq.t_control(2));
                freq{f_count, 1}.hilspctrm = freq{f_count, 1}.hilspctrm(e, :, idx_pre);
                if isempty(freq{f_count, 1}.hilspctrm); continue; end
                freq{f_count, 1}.time = linspace(post_st, post_en, size(freq{f_count, 1}.hilspctrm, 3));
                freq{f_count, 1}.label = freq{f_count, 1}.label(e);


                % store amplitude of post-window
                if strcmp('start', timelocking_type)
                    idx_start_apnoea = find(X.freq.time == 0); % idx_start_apnoea = find(X.freq.time > -0.5 / 60 & X.freq.time < 0.5 / 60);
                    idx_post = idx_start_apnoea + post_st * fs_eeg : idx_start_apnoea + post_en * fs_eeg;
                elseif strcmp('end', timelocking_type)
                    [~, idx_end_apnoea] = min(abs(X.freq.time - X.t_length));
                    idx_post = (idx_end_apnoea + post_st * fs_eeg : idx_end_apnoea + post_en * fs_eeg) - 1;
                end

                if any(idx_post <= 0); continue; end
                freq{f_count, 2}.hilspctrm = freq{f_count, 2}.hilspctrm(e, :, idx_post);
                freq{f_count, 2}.time = freq{f_count, 1}.time;
                freq{f_count, 2}.label = freq{f_count, 2}.label(e);

                f_count = f_count + 1;

            end

        end
        freq = freq(~cellfun(@isempty, freq(:, 1)), :);


        % run the statistics on the time-frequency-amplitude
        % representations
        cfg = [];
        cfg.method = 'analytic';
        cfg.correctm = 'fdr';
        cfg.alpha = 0.01;
        stat_type = 'fdr';
        cfg.parameter = 'hilspctrm';
        cfg.statistic = 'ft_statfun_depsamplesT';
        cfg.channel = eeg_channels(e);
        cfg.correcttail = 'alpha';
        cfg.tail = 0;
        cfg.neighbours = [];

        design = zeros(2, size(freq, 1) * 2);
        design(1, :) = [1 : size(freq, 1), 1 : size(freq, 1)];
        design(2, :) = [ones(1, size(freq, 1)), repmat(2, [1, size(freq, 1)])];

        cfg.design = design;
        cfg.ivar = 2;
        cfg.uvar = 1;

        stat{e, 1} = ft_freqstatistics(cfg, freq{:});
        stat{e, 1}.cfg = [];


        % plot probability of Hilbert
        figure(1 + e * 10 + time_window * 100);
        po = get(gcf, 'position');
        set(gcf, 'position', [po(1:2), 400, 100], 'name', sprintf('prob_hilbert_%s_%s_%s_%d_%d_sec_%s', ibi_type_label, timelocking_type, eeg_channels{e}, post_st, post_en, stat_type))

        cfg = [];
        cfg.parameter = 'prob';
        cfg.channel = stat{e, 1}.label;
        cfg.colormap = cmocean('matter');
        cfg.ylim = [freq_range(1), freq_range(end)];
        cfg.zlim = [1e-5, 0.01];

        ft_singleplotTFR(cfg, stat{e, 1})
        title('')
        set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [], 'colorscale', 'log')
        colorbar off


        % plot mean Hilbert differences
        figure(2 + e * 10 + time_window * 100);
        po = get(gcf, 'position');
        set(gcf, 'position', [po(1:2), 400, 100], 'name', sprintf('diff_hilbert_mean_%s_%s_%s_%d_%d_sec_%s', ibi_type_label, timelocking_type, eeg_channels{e}, post_st, post_en, stat_type))

        cfg = [];
        cfg.parameter = 'hilspctrm';
        hil_pre = ft_freqgrandaverage(cfg, freq{:, 1});
        hil_post = ft_freqgrandaverage(cfg, freq{:, 2});
        hil_pre.mask = stat{e, 1}.mask;
        hil_pre.hilspctrm = (hil_post.hilspctrm - hil_pre.hilspctrm);

        cfg = [];
        cfg.parameter = 'hilspctrm';
        cfg.maskparameter = 'mask';
        cfg.colormap = flipud(brewermap(64, 'RdBu'));
        cfg.ylim = [freq_range(1), freq_range(end)];
        cfg.zlim = [-0.1, 0.1];

        ft_singleplotTFR(cfg, hil_pre)
        title('')
        set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
        colorbar off
        drawnow

        datetime
    end

    % save stat output
    output_file = fullfile(pwd, 'output/hilbert_stats', ibi_type_label, sprintf('hilbert_stats_%s_%d_%d_sec', timelocking_type, post_st, post_en));
    save(output_file, 'stat', '-v7.3')

end