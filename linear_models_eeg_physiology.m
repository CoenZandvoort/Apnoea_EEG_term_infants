% Linear mixed effects models on the heart rate, oxygen aturation, pause 
% duration, PMA, and sleep states predicting EEG amplitude changes.
%
% CZ, Mar-2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(fullfile('~/Documents/MATLAB/toolboxes/Violinplot-Matlab-master/'))

% get sheet overview
sheet = importdata('data_overview.xlsx');
subj_labels = sheet.textdata.to_include(3 : end, 1);
pma = sheet.data.to_include(:, 6) + sheet.data.to_include(:, 7) / 7;

for eeg_window = 1 : 2

    switch eeg_window

        case 1
            ibi_type_label = 'ibi_apnoeas';
            timelocking_type = 'end';
            post_st_eeg = -5; % sec
            post_en_eeg = 5; % sec
            t_crop_window = 20; % sec

        case 2
            ibi_type_label = 'ibi_5_15_sec';
            timelocking_type = 'end';
            post_st_eeg = -5; % sec
            post_en_eeg = 5; % sec
            t_crop_window = 0; % sec

    end


    post_st_vs = -5; % sec
    post_en_vs = 60; % sec
    post_type_vs = 'min';


    output_file = fullfile(pwd, 'output/linear_models', ibi_type_label, sprintf('lme_amp_eeg_%d_%d_sec_vs_%d_%d_sec_%s.mat', post_st_eeg, post_en_eeg, post_st_vs, post_en_vs, post_type_vs));
    if ~exist(output_file, 'file')

        breath_counter = 1;
        eeg = NaN(1, 1); eeg_pre = NaN(1, 1); eeg_post = NaN(1, 1);
        sats_m = NaN(1, 1); sats_pre = NaN(1, 1); sats_post = NaN(1, 1);
        hr_m = NaN(1, 1); hr_pre = NaN(1, 1); hr_post = NaN(1, 1);
        sleep_state_m = cell(1, 1);
        pause_length_m = NaN(1, 1);
        pma_m = NaN(1, 1);
        subj = NaN(1, 1);


        fs_eeg = 250; % Hz
        for s = 1 : numel(subj_labels)

            % load physiological data
            try
                data = load(fullfile('/Volumes/Coen/MATLAB/EEG_apnoea_Leuven/output/timelock', ibi_type_label, [subj_labels{s}(1 : end - 4), '.mat'])); % data = load(fullfile(pwd, 'output/timelock', ibi_type_label, [subj_labels{s}(1 : end - 4), '.mat']));
                sleep_states = readtable(fullfile(pwd, 'output/sleep_states', [subj_labels{s}(1 : end - 4), '_sleep.xlsx']));
            catch err
                warning(err.message)
            end


            for pause = 1 : numel(data.data)

                % load Hilbert amplitude and physiological data
                amp_folder = fullfile('/Volumes/Coen/MATLAB/EEG_apnoea_Leuven/output/hilbert_amp', ibi_type_label); % amp_folder = fullfile(pwd, 'output/hilbert_amp', ibi_type_label);
                try
                    X = load(fullfile(amp_folder, sprintf('hilbert_%s_%d.mat', subj_labels{s}(1 : end - 4), pause)));
                catch err
                    warning(err.message)
                    continue
                end


                % get pre- and post-amplitudes
                idx_pre_eeg = (X.freq.time >= X.freq.t_control(1) & X.freq.time <= X.freq.t_control(2) - t_crop_window);
                pre_eeg = X.freq.hilspctrm(end, :, idx_pre_eeg);

                if strcmp(timelocking_type, 'start')
                    [~, idx_end_apnoea] = min(abs(X.freq.time));
                elseif strcmp(timelocking_type, 'end')
                    [~, idx_end_apnoea] = min(abs(X.freq.time - X.t_length));
                end
                idx_post_eeg = (idx_end_apnoea + post_st_eeg * fs_eeg : idx_end_apnoea + post_en_eeg * fs_eeg) - 1;
                post_eeg = X.freq.hilspctrm(end, :, idx_post_eeg);


                % get EEG amplitude difference 
                eeg(breath_counter, 1) = mean(post_eeg(post_eeg ~= 0)) - mean(pre_eeg(pre_eeg ~= 0));

                eeg_post(breath_counter, 1) = mean(post_eeg(post_eeg ~= 0));
                eeg_pre(breath_counter, 1) = mean(pre_eeg(pre_eeg ~= 0));


                % define offset
                if strcmp(timelocking_type, 'start')
                    t_start = 0; % sec
                elseif strcmp(timelocking_type, 'end')
                    t_start = X.t_length; % sec
                end


                % get saturation
                idx_sheet_sats = find(startsWith(sheet.textdata.to_include(2, :), 'SATS'));
                sats_label = sheet.textdata.to_include(s + 2, idx_sheet_sats);
                idx_sats_label = find(strcmp(data.data{pause}.labels, sats_label));

                if ~isempty(idx_sats_label)
                    idx_extremes = (data.data{pause}.data{idx_sats_label} < 60 | data.data{pause}.data{idx_sats_label} > 100);
                    data.data{pause}.data{idx_sats_label}(idx_extremes) = NaN;

                    sats_pre(breath_counter, 1) = mean(data.data{pause}.data{idx_sats_label}(data.data{pause}.time{idx_sats_label} >= X.freq.t_control(1) + t_crop_window & data.data{pause}.time{idx_sats_label} <= X.freq.t_control(2)), 'omitnan');
                    if strcmp(post_type_vs, 'mean')
                        sats_post(breath_counter, 1) = mean(data.data{pause}.data{idx_sats_label}(data.data{pause}.time{idx_sats_label} >= t_start + post_st_vs & data.data{pause}.time{idx_sats_label} <= t_start + post_en_vs), 'omitnan');
                    elseif strcmp(post_type_vs, 'min')
                        sats_post(breath_counter, 1) = min(data.data{pause}.data{idx_sats_label}(data.data{pause}.time{idx_sats_label} >= t_start + post_st_vs & data.data{pause}.time{idx_sats_label} <= t_start + post_en_vs));
                    end

                    sats_m(breath_counter, 1) = sats_post(breath_counter, 1) - sats_pre(breath_counter, 1);
                else
                    sats_m(breath_counter, 1) = NaN;
                end


                % get heartrate
                idx_sheet_hr = find(startsWith(sheet.textdata.to_include(2, :), 'HR'));
                hr_label = sheet.textdata.to_include(s + 2, idx_sheet_hr);
                idx_hr_label = find(strcmp(data.data{pause}.labels, hr_label));


                if ~isempty(idx_hr_label)
                    % in some cases, the heart rate has some very sudden
                    % changes (e.g., changes greater than 50 bpm/s). We will
                    % therefore check for that.
                    idx_spline = find(abs(diff(data.data{pause}.data{idx_hr_label})) < 50);
                    if ~isempty(idx_spline)
                        data.data{pause}.data{idx_hr_label} = spline(idx_spline, data.data{pause}.data{idx_hr_label}(idx_spline), 1 : numel(data.data{pause}.data{idx_hr_label}));
                    end


                    idx_extremes = (data.data{pause}.data{idx_hr_label} < 40 | data.data{pause}.data{idx_hr_label} > 230);
                    data.data{pause}.data{idx_hr_label}(idx_extremes) = NaN;

                    hr_pre(breath_counter, 1) = mean(data.data{pause}.data{idx_hr_label}(data.data{pause}.time{idx_hr_label} >= X.freq.t_control(1) + t_crop_window & data.data{pause}.time{idx_hr_label} <= X.freq.t_control(2)), 'omitnan');
                    if strcmp(post_type_vs, 'mean')
                        hr_post(breath_counter, 1) = mean(data.data{pause}.data{idx_hr_label}(data.data{pause}.time{idx_hr_label} >= t_start + post_st_vs & data.data{pause}.time{idx_hr_label} <= t_start + post_en_vs), 'omitnan');
                    elseif strcmp(post_type_vs, 'min')
                        hr_post(breath_counter, 1) = min(data.data{pause}.data{idx_hr_label}(data.data{pause}.time{idx_hr_label} >= t_start + post_st_vs & data.data{pause}.time{idx_hr_label} <= t_start + post_en_vs));
                    end

                    hr_m(breath_counter, 1) = hr_post(breath_counter, 1) - hr_pre(breath_counter, 1);
                else
                    hr_m(breath_counter, 1) = NaN;
                end


                % get sleep state
                if exist('sleep_states', 'var')
                    
                    idx_sleep_state = find(sleep_states.start_time - data.data{pause}.info.timing <= 0, 1, 'last');
                    sleep_state_m(breath_counter, 1) = sleep_states.sleep_label_cnn(idx_sleep_state);

                end
                

                % get remaining variables
                pause_length_m(breath_counter, 1) = data.data{pause}.info.length;
                pma_m(breath_counter, 1) = pma(s);
                subj(breath_counter, 1) = s;
                breath_counter = breath_counter + 1;

            end

            clear sleep_states

        end


        % run linear mixed-effects model (here, we test an lme model with sats,
        % hr, apnoea length, pma, or sleep states as main effect and subj 
        % as random effect of subject for the slope, which will be 
        % something like:
        % eeg = Beta0 + Beta1 * main effect + Gamma1 * subj + error).
        lme = cell(5, 1);

        for mdl = 1 : 5

            switch mdl
                case 1 % 'eeg ~ 1 + hr + (1 + hr | infant)';
                    response_var = eeg;
                    pred_var = hr_m;
                    subj_var = subj;

                case 2 % 'eeg ~ 1 + sats + (1 + sats | infant)';
                    response_var = eeg;
                    pred_var = sats_m;
                    subj_var = subj;

                case 3 % 'eeg ~ 1 + breathing_pause_length + (1 + breathing_pause_length | infant)';
                    response_var = eeg;
                    pred_var = pause_length_m; 
                    subj_var = subj;

                case 4 % 'eeg ~ 1 + pma + (1 + pma | infant)';
                    
                    unique_subj = unique(subj);
                    response_var = NaN(numel(unique_subj), 1);
                    pred_var = NaN(numel(unique_subj), 1);
                    subj_var = NaN(numel(unique_subj), 1);
                    
                    for s = 1 : numel(unique_subj)
                        response_var(s, 1) = mean(eeg(unique_subj(s) == subj), 'omitnan');
                        pred_var(s, 1) = mean(pma_m(unique_subj(s) == subj), 'omitnan');
                        subj_var(s, 1) = unique_subj(s);
                    end

                case 5 % 'eeg ~ 1 + sleep_state + (1 + sleep_state | infant)';

                    response_var = eeg(~cellfun(@isempty, sleep_state_m));
                    pred_var = categorical(sleep_state_m(~cellfun(@isempty, sleep_state_m))); 
                    subj_var = subj(~cellfun(@isempty, sleep_state_m));

            end

            stat_model_label = 'EEG ~ 1 + phys + (1 + phys | infant)';
            tbl = table(response_var, pred_var, subj_var, 'VariableNames', {'EEG', 'phys', 'infant'});

            lme{mdl} = fitlme(tbl, stat_model_label);

        end


        % save linear model output
        save(output_file, 'lme', 'eeg', 'eeg_pre', 'eeg_post', 'hr_m', 'hr_pre', 'hr_post', 'sats_m', 'sats_pre', 'sats_post', 'pause_length_m', 'subj', 'pma_m', 'sleep_state_m');
        clear lme eeg eeg_pre eeg_post hr_m hr_pre hr_post sats_m sats_pre sats_post pause_length_m subj pma_m sleep_state_m

    end


    % plot linear mixed effect models
    mdl_main_effect = {'hr_m', 'sats_m', 'pause_length_m', 'pma_m_avg', 'sleep_state_m'};
    if strcmp(ibi_type_label, 'ibi_apnoeas')
        edges = {-160 : 80 : 80; -40 : 10 : 25; 15 : 5 : 70; 30 : 5 : 50; 1 : 4};
    elseif strcmp(ibi_type_label, 'ibi_5_15_sec')
        edges = {-160 : 80 : 80; -40 : 10 : 25; 5 : 15; 30 : 5 : 50; 1 : 4};
    end

    for m = 1 : numel(mdl_main_effect)

        % load linear mixed effect models
        lme = load(output_file);


        % sats/hr/breathing_pause_length --> eeg
        figure(1 + m * 10 + eeg_window * 100);
        po = get(gcf, 'position');
        set(gcf, 'position', [po(1:2), 400, 400], 'name', sprintf('glme_EEG_%s_%d_%d_sec_vs_%d_%d_sec_%s_%s', ...
            mdl_main_effect{m}, post_st_eeg, post_en_eeg, post_st_vs, post_en_vs, post_type_vs, ibi_type_label));


        % get average ages
        if strcmp(mdl_main_effect{m}, 'pma_m_avg')

            unique_subj = unique(lme.subj);
            lme.pma_m_avg = NaN(numel(unique_subj), 1);
            lme.eeg_avg = NaN(numel(unique_subj), 1);

            for s = 1 : numel(unique_subj)
                lme.eeg_avg(s, 1) = mean(lme.eeg(unique_subj(s) == lme.subj), 'omitnan');
                lme.pma_m_avg(s, 1) = mean(lme.pma_m(unique_subj(s) == lme.subj), 'omitnan');
            end

            lme.pma_m = lme.pma_m_avg;
            lme.eeg = lme.eeg_avg;
            
        elseif strcmp(mdl_main_effect{m}, 'sleep_state_m')

            response_var = [];
            pred_var = [];
            subj_var = [];
            sleep_state_2_include = {'ASI', 'LVI', 'QS HVS', 'QS TA'};
            for s = 1 : numel(sleep_state_2_include)
                idx_sleep = find(strcmp(lme.sleep_state_m , sleep_state_2_include{s}));

                response_var = [response_var; lme.eeg(idx_sleep)];
                pred_var = [pred_var; lme.sleep_state_m(idx_sleep)];
                subj_var = [subj_var; lme.subj(idx_sleep)];
            end


            stat_model_label = 'EEG ~ 1 + phys + (1 + phys | infant)';
            tbl = table(response_var, pred_var, subj_var, 'VariableNames', {'EEG', 'phys', 'infant'});
            lme.lme{m} = fitlme(tbl, stat_model_label);

        end


        if any(strcmp(mdl_main_effect{m}, {'hr_m', 'sats_m', 'pause_length_m', 'pma_m_avg'}))
            % plot data
            scatter(lme.(mdl_main_effect{m}), mean(lme.eeg, 2, 'omitnan'), 'Filled', 'SizeData', 16, ...
                'MarkerFaceColor', [0.83, 0.14, 0.14], 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.2)
            hold on
            ci_mean = glmval(lme.lme{m}.Coefficients.Estimate, sort(lme.(mdl_main_effect{m})), 'identity', lme.lme{m});
            plot(sort(lme.(mdl_main_effect{m})), ci_mean, 'k', 'linewidth', 1.5) % mean regression estimate
            xlim([edges{m}(1), edges{m}(end)])
            ylim([-0.7, 0.7])
            grid


            % plot histograms
            figure;
            po = get(gcf, 'position');
            set(gcf, 'position', [po(1 : 2), 400, 100], 'name', sprintf('hist_%s_%s', mdl_main_effect{m}, ibi_type_label))
            h = histogram(lme.(mdl_main_effect{m}), 'BinWidth', diff([edges{m}(1), edges{m}(end)]) / 50, 'BinLimits', [edges{m}(1), edges{m}(end)], 'FaceColor', [0.5, 0.5, 0.5]);
            hold on
            plot([median(lme.(mdl_main_effect{m}), 'omitnan'), median(lme.(mdl_main_effect{m}), 'omitnan')], [0, max(h.Values)], 'Color', [0.0, 0.0, 0.7], 'LineWidth', 1.5)
            set(gca, 'XTick', [], 'XTickLabels', [], 'YTick', [], 'YTickLabels', [], 'Visible', 'Off')


            figure;
            po = get(gcf, 'position');
            set(gcf, 'position', [po(1 : 2), 400, 100], 'name', sprintf('hist_eeg_%s_%s', mdl_main_effect{m}, ibi_type_label))
            h = histogram(lme.eeg, 'BinWidth', diff([-1, 1]) / 50, 'BinLimits', [-1, 1], 'FaceColor', [0.5, 0.5, 0.5]);
            hold on
            plot([median(lme.eeg, 'omitnan'), median(lme.eeg, 'omitnan')], [0, max(h.Values)], 'Color', [0.0, 0.0, 0.7], 'LineWidth', 1.5)
            set(gca, 'XTick', [], 'XTickLabels', [], 'YTick', [], 'YTickLabels', [], 'Visible', 'Off')
            xlim([-0.7, 0.7])
        
        elseif strcmp(mdl_main_effect{m}, 'sleep_state_m')
            
            num_idx_sleep = NaN(numel(pred_var), 1);
            vi = violinplot(response_var(response_var > -0.7), pred_var(response_var > -0.7));
            for v = 1 : numel(vi)
                vi(v).ViolinPlot.FaceColor = [0.83, 0.14, 0.14];
                vi(v).ScatterPlot.MarkerFaceColor = [0.50, 0.50, 0.50];
            end
            xlim([0.5, 4.5])
            ylim([-0.7, 0.7])
            
            for s = 1 : numel(sleep_state_2_include)
                idx_sleep = strcmp(sleep_state_2_include{s}, pred_var);
                num_idx_sleep(idx_sleep, 1) = s;
            end
            

            % plot histograms
            figure;
            po = get(gcf, 'position');
            set(gcf, 'position', [po(1 : 2), 400, 100], 'name', sprintf('hist_%s_%s', mdl_main_effect{m}, ibi_type_label))
            histogram(num_idx_sleep, 'FaceColor', [0.5, 0.5, 0.5]);
            set(gca, 'XTick', [], 'XTickLabels', [], 'YTick', [], 'YTickLabels', [], 'Visible', 'Off')

            figure;
            po = get(gcf, 'position');
            set(gcf, 'position', [po(1 : 2), 400, 100], 'name', sprintf('hist_eeg_%s_%s', mdl_main_effect{m}, ibi_type_label))
            h = histogram(response_var, 'BinWidth', diff([-1, 1]) / 50, 'BinLimits', [-1, 1], 'FaceColor', [0.5, 0.5, 0.5]);
            hold on
            plot([median(response_var, 'omitnan'), median(response_var, 'omitnan')], [0, max(h.Values)], 'Color', [0.0, 0.0, 0.7], 'LineWidth', 1.5)
            set(gca, 'XTick', [], 'XTickLabels', [], 'YTick', [], 'YTickLabels', [], 'Visible', 'Off')
            xlim([-0.7, 0.7])

        end

        fprintf('%s:\tBeta: %.08f; p: %.08f\n', mdl_main_effect{m}, lme.lme{m}.Coefficients.Estimate(2), lme.lme{m}.Coefficients.pValue(2))
        
    end

    set(findall(0, 'type', 'axes'), 'FontName', 'Times', 'Fontsize', 16, 'TickDir', 'out', 'box', 'off', 'ticklength', [0.010, 0.010], 'Linewidth', 2)

end