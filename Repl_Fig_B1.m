choose_vars   = [12 7 2 3 4 5 9 10 8 11];
[data_adj,vmed_adj,vars_adj] = est_irf_fev(nirf,bound,ndraws, nconst, prior, begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);
     

%% Repeat with unadjusted tfp. 
choose_vars = [12 13 2 3 4 5 9 10 8 11]; % replace adjusted TFP with unadjusted TFP, ordered 37. 
[data_unadj,vmed_unadj,vars_unadj] = est_irf_fev(nirf,bound,ndraws, nconst, prior, begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);
%% ==============================================================================================================================================================================================
%----------------------------------------------------------------------
% Plot the responses of KPSS and TFP adjusted. 
%----------------------------------------------------------------------
fonttype          ='Arial';
ftsizeaxis        = 11;
titlefontsize     = 10;
[nvars,temp]=size(vars_adj);
[nirf,temp,temp]=size(data_adj);
%% Plot both together as Robustness for the Appendix:
figB1_appendix=figure('Color',[0.9412 0.9412 0.9412],'Position',[1 1 800-100 600-100],'Name','FIGB1_CGV_REStat');

figure(figB1_appendix)
n = 2; % TFP ordered second.
hold on 
h(1) = plot(1:nirf,data_adj(:,nvars+n,1), 'b', 'LineWidth',2.5);hold on                                        
h(2) = plot(1:nirf,data_adj(:,nvars+n,2),'r--','LineWidth',2.5);hold on                            
h(3) = plot(1:nirf,data_adj(:,nvars+n,3),'r--','LineWidth',2.5);hold on
grpyat = [(1:nirf)' data_unadj(1:nirf,nvars+n,2); (nirf:-1:1)' data_unadj(nirf:-1:1,nvars+n,3)];
h(4)= patch(grpyat(:,1),grpyat(:,2),[0.4 0.3 0.2],'edgecolor','none'); 
h(5) = plot(1:nirf,data_unadj(:,nvars+n,1), '*-', 'color', [0.4 0.3 0.2], 'LineWidth',1.5);
hold on                                        
plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);
ylabel('percent','FontSize',12)
xlabel('Horizon (quarters)','FontSize',12)
set(gca,'XTick',[0;4;8;12;16;20])
set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
set(gca, 'FontName', fonttype);
set(gca, 'FontSize', ftsizeaxis);
set(gca,'Layer','top');
alpha(.1)
axis tight
hold off;
set(legend, 'NumColumns',2)
legend('Orientation','horizontal');%,'Location','southoutside')
hlegend=legend(h([1 5 4 3]),'Utilization-adjusted TFP', 'Unadjusted TFP',  'Confidence bands', 'Confidence bands');

 
 %%                                  
clear data_adj data_unadj