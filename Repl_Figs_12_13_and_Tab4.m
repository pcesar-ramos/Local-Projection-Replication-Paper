choose_vars       = [7 2 3 4 5 9 10 8 11];
[SVAR] = est_proxy(begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars,instrument_var,ext_shocks_ind,rep_proxy);
%% ==============================================================================================================================================================================================
%----------------------------------------------------------------------
% Plot the responses of TFP adjusted and other variables with a
% Bayesian Proxy SVAR (Caldara and Herbst (2019))
%----------------------------------------------------------------------
fonttype          ='Arial';
ftsizeaxis        = 11;
titlefontsize     = 10;
graph_opt.font_num = 10;
Horizon = size(SVAR.LtildeFull,2);
H = Horizon -1;
varSelec = SVAR.varSelec;

%% FIGURE 12 OF THE PAPER.
a=figure('Color',[0.9412 0.9412 0.9412],'Position',[1 1 800-100 600-100],'Name','FIG12_CGV_RESTAT');

figure(a);
plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(1),1:Horizon,3,1))*100,'LineWidth',2)
hline(0,':k')
hold on
plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(1),1:Horizon,1,1))*100,'r--','LineWidth',2)
plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(1),1:Horizon,5,1))*100,'r--','LineWidth',2)
ylabel('percent','FontSize',12)
xlabel('Horizon (quarters)','FontSize',12)
title(num2str(SVAR.vars(1,:)),'Interpreter','tex','FontSize',titlefontsize)

set(gca,'XTick',[0;4;8;12;16;20])
set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
set(gca, 'FontName', fonttype);
set(gca, 'FontSize', ftsizeaxis);
set(gca,'Layer','top');

box off
axis tight
hold off;

%% FIGURE 13 OF THE PAPER.
b=figure('Color',[0.9412 0.9412 0.9412],'Position',[1 1 800-100 600-100],'Name','FIG13_CGV_RESTAT');

figure(b);
for ii = 2:length(varSelec) % Variable
    subplot(4,2,ii-1)
    plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(ii),1:Horizon,3,1))*100,'LineWidth',2)
    hline(0,':k')
    hold on
    plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(ii),1:Horizon,1,1))*100,'r--','LineWidth',2)
    plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(ii),1:Horizon,5,1))*100,'r--','LineWidth',2)
    ylabel('percent','FontSize',12)
    xlabel('Horizon (quarters)','FontSize',12)
    title(num2str(SVAR.vars(varSelec(ii),:)),'Interpreter','tex','FontSize',titlefontsize)
    
    set(gca,'XTick',[0;4;8;12;16;20])
    set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
    set(gca, 'FontName', fonttype);
    set(gca, 'FontSize', ftsizeaxis);
    set(gca,'Layer','top');
    
    box off
    axis tight
    hold off;
end

%% CORRELATIONS WITH EXTERNAL SHOCKS - TABLE 4.
ext_shocks = SVAR.ext_shocks;
shocks_var = SVAR.shocks_var;
iv_mm = SVAR.mm;

corr_var = NaN(size(ext_shocks,2),1); pv_var = NaN(size(ext_shocks,2),1);
for jj=1:size(corr_var,1)
        validIndices = ~isnan(ext_shocks(:,jj));
        [corr_temp,pv_temp,~,~] = corrcoef([iv_mm(validIndices) ext_shocks(validIndices,jj)]);
        corr_var(jj,1) = corr_temp(2,1);
        pv_var(jj,1) = pv_temp(2,1);
end
[~, ~, ~, pv_var_adj]=fdr_bh(pv_var,0.05,'pdep','no');



VarNames = {'Correlation', 'Pvalue'};
disp('TABLE 4: ===============================================================================================')
 for i=1:size(ext_shocks,2)
 disp(['External shock:', num2str(shocks_var(i,:))]) %Interpreter','tex','FontSize',titlefontsize))
 disp(table(round(corr_var(i),2),round(pv_var_adj(i),2),'VariableNames',VarNames))
 end
disp('========================================================================================================')
