
load ('DATA/shocks_agg1.mat') 
MA = shocks2;
MA = MA(1:end,1:2); % Barsky and Sims news shock and CGV news shock
%% Compute filtered variables
MA_HF = bandpass(MA,2,6) ;
MA_BC = bandpass(MA,6,40) ;
MA_MC = bandpass(MA,32,80) ;
MA_HBC= bandpass(MA,16,40) ;
%% TABLE 3:
temp = corrcoef(MA);
TABLE(1,:) = temp(1,:);
temp = corrcoef(MA_HF);
TABLE(2,:) = temp(1,:); % High Frequencies correlation
temp = corrcoef(MA_BC);
TABLE(3,:) = temp(1,:);
TABLE=round(TABLE,2);


disp ('TABLE 3 correlations ================================================')
[TABLE(1,2);TABLE(2,2);TABLE(3,2)] 
%% PLOT without HF (that is fluctiation above 6 quarters)

corr_shocks=figure('Color',[0.9412 0.9412 0.9412],'Position',[1 1 800-100 600-100],'Name','FIG8_CGV_RESTAT');   

figure(corr_shocks)
tt = 1962:0.25:2010.75;
subplot(2,1,1)
plot(tt, MA(:,1), tt, MA(:,2), '--',  'Linewidth',1.5)
ylim([-3.5 3.5])
xlim([1962 2010.5])
title ('Unfiltered Shock Series')
subplot(2,1,2)
MA_noHF=MA-MA_HF;
plot(tt, MA_BC(:,1), tt,MA_BC(:,2), '--',  'Linewidth',1.5)
ylim([-3.5 3.5])
xlim([1962 2010.75])
title ('Filtered Shock Series: Extended Business-Cycle Frequency')
legend('Patent-Based News Shocks','Barsky and Sims News Shocks')


%%