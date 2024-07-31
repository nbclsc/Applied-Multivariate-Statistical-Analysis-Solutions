% Exercise 4.31
% Multiple-Sclerosis Data
data_4_31 = readmatrix(fullfile(data_folder, 'Table1.6.xlsx'));
no_ms = data_4_31(data_4_31(:,6) == 0, 1:5);
yes_ms = data_4_31(data_4_31(:,6) == 1,1:5);

% Compute simulation of critical values for Q-Q correlation test.
[~, crit_no_ms] = ppcc_simulation(height(no_ms), 1000000, [0.01 0.05 0.10]);
[~, crit_yes_ms] = ppcc_simulation(height(yes_ms), 1000000, [0.01 0.05 0.10]);

% Store the Q-Q correlation values.
r_Q = zeros(5,2);
% Store the data related to Q-Q computation.
qq_no_ms = zeros(height(no_ms), width(no_ms)*3);
qq_yes_ms = zeros(height(yes_ms), width(yes_ms)*3);

for i=1:width(no_ms)
    offset = 3*(i-1) + 1;
    [qq_no_ms(:,offset:offset+2), r_Q(i,1)] = corr_q_q(no_ms(:,i));
    [qq_yes_ms(:,offset:offset+2), r_Q(i,2)] = corr_q_q(yes_ms(:,i));
end

% Check which MS negative are variables NOT normally distributed.
r_Q(:,1) < crit_no_ms(1)
% Check which MS positve are variables NOT normally distributed.
r_Q(:,2) < crit_yes_ms(1)

% MS negative: 1, 2, 3, and 5 needs univariate transform. 4 is okay.
% MS positive: 3 and 5 need univariate transform. 1, 2, and 4 are okay.

cols = ["Age", "S1L + S1R", "|S1L - S1R|", "S2L + S2R", "|S2L - S2R|"];

% Q-Q plots for MS negative.
for i=1:length(cols)
    offset = 3*(i-1) + 1;
    if i == 1
        figure1 = figure();
        % Want a narrow and tall output image.
        screenSize = get(0, 'ScreenSize');
        desiredWidth = 600;
        figure1.Position = [100, 1, desiredWidth, screenSize(4)]; % [left, bottom, width, height]
    end
    subplot1 = subplot(5,1,i,'Parent',figure1);
    hold(subplot1,'on');
    scatter(qq_no_ms(:, offset+2), qq_no_ms(:, offset))
    title(cols(i))
    xlabel("q_{(j)}")
    ylabel("x_{(j)}")
    hold(subplot1,'off');
    if i == length(cols)
        sgtitle("Q-Q Plots MS Negative");
    end
end
exportgraphics(gcf, append('.\', 'sol4.31.qq.msneg', '.png'))

% Q-Q plots for MS positive.
for i=1:length(cols)
    offset = 3*(i-1) + 1;
    if i == 1
        figure1 = figure();
        % Want a narrow and tall output image.
        screenSize = get(0, 'ScreenSize');
        desiredWidth = 600;
        figure1.Position = [100, 1, desiredWidth, screenSize(4)]; % [left, bottom, width, height]
    end
    subplot1 = subplot(5,1,i,'Parent',figure1);
    hold(subplot1,'on');
    scatter(qq_yes_ms(:, offset+2), qq_yes_ms(:, offset))
    title(cols(i))
    xlabel("q_{(j)}")
    ylabel("x_{(j)}")
    hold(subplot1,'off');
    if i == length(cols)
        sgtitle("Q-Q Plots MS Positive");
    end
end
exportgraphics(gcf, append('.\', 'sol4.31.qq.mspos', '.png'))

% If the Q-Q correlation is less than the critical value, perform the
% univariate power transformation.
ms_power = zeros(5,2);
for i=1:length(cols)
    % Power transformation for MS Negative group.
    if r_Q(i,1) < crit_no_ms(1)
        ms_power(i,1) = var_power_transform(no_ms(:,i), ...
            sprintf('MS Negative %s', cols(i)), ...
            append(".\", "sol4.31.msneg.power.", string(i), ".png"));
    end
    % Power transformation for MS Positive group.
    if r_Q(i,2) < crit_yes_ms(1)
        ms_power(i,2) = var_power_transform(yes_ms(:,i), ...
            sprintf('MS Positive %s', cols(i)), ...
            append(".\", "sol4.31.mspos.power.", string(i), ".png"));
    end
end

% Perform the power transformations.
no_ms_trans = [no_ms(:,1).^(-0.4) no_ms(:,2).^(-3.5) no_ms(:,3).^(0.25) no_ms(:,4) no_ms(:,5).^(0.2)];
yes_ms_trans = [yes_ms(:,1) yes_ms(:,2) yes_ms(:,3).^(0.25) yes_ms(:,4) yes_ms(:,5).^(0.2)];

% Store the Q-Q correlation values for the power transformations.
r_Q_trans = zeros(5,2);
% Store the data related to Q-Q computation for the power transformations.
qq_no_ms_trans = zeros(height(no_ms), width(no_ms)*3);
qq_yes_ms_trans = zeros(height(yes_ms), width(yes_ms)*3);

for i=1:width(no_ms)
    offset = 3*(i-1) + 1;
    [qq_no_ms_trans(:,offset:offset+2), r_Q_trans(i,1)] = corr_q_q(no_ms_trans(:,i));
    [qq_yes_ms_trans(:,offset:offset+2), r_Q_trans(i,2)] = corr_q_q(yes_ms_trans(:,i));
end

% Check which MS negative are variables NOT normally distributed.
r_Q_trans(:,1) < crit_no_ms(1)
% Check which MS positve are variables NOT normally distributed.
r_Q_trans(:,2) < crit_yes_ms(1)

% Check proportion of zeros.
sum(no_ms(:,3)==0) / height(no_ms)
sum(no_ms(:,5)==0) / height(no_ms)

sum(yes_ms(:,3)==0) / height(yes_ms)
sum(yes_ms(:,5)==0) / height(yes_ms)

% Apply (4-29). These are all true, so they're still not normally dist.
abs(sum((no_ms_trans(:,1) >= (mean(no_ms_trans(:,1)) - 1*var(no_ms_trans(:,1)))) & (no_ms_trans(:,1) <= (mean(no_ms_trans(:,1)) + 1*var(no_ms_trans(:,1))))) / height(no_ms_trans) - 0.683) > 1.396/sqrt(height(no_ms_trans))
abs(sum((no_ms_trans(:,3) >= (mean(no_ms_trans(:,3)) - 1*var(no_ms_trans(:,3)))) & (no_ms_trans(:,3) <= (mean(no_ms_trans(:,3)) + 1*var(no_ms_trans(:,3))))) / height(no_ms_trans) - 0.683) > 1.396/sqrt(height(no_ms_trans))
abs(sum((no_ms_trans(:,5) >= (mean(no_ms_trans(:,5)) - 1*var(no_ms_trans(:,5)))) & (no_ms_trans(:,5) <= (mean(no_ms_trans(:,5)) + 1*var(no_ms_trans(:,5))))) / height(no_ms_trans) - 0.683) > 1.396/sqrt(height(no_ms_trans))