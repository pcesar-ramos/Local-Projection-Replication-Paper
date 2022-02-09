function vm_plot_irf_bvar(SVAR)

graph_opt.font_num = 10;

Horizon = size(SVAR.LtildeFull,2);
H = Horizon -1;
nIV = size(SVAR.i_var_instr,2);
nshockplot = nIV;
varSelec = SVAR.varSelec;

for jj = 1:nshockplot % Shock
    % Plot IRFs TFP
    figure('Units','centimeters','Position', [0,0, 20, 10],'Name','News_IRF_TFP');

        plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(1),1:Horizon,3,jj))*100,'LineWidth',2)
        hline(0,':k')
        hold on
        plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(1),1:Horizon,1,jj))*100,'r--', ...
             'LineWidth',2)
        plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(1),1:Horizon,5,jj))*100,'r--', ...
             'LineWidth',2)
        title(SVAR.i_var_str_names(:,varSelec(1)),'FontSize',graph_opt.font_num,'FontWeight','bold')
        axis tight
        ylabel('percent','FontSize',12)
        xlabel('Horizon (quarters)','FontSize',12)
        set(gca,'XTick',0:4:H)
        set(gca,'LineWidth',2)
        hold off
        box off
    
    % Plot IRFs others
    figure('Units','centimeters','Position', [0,0, 30, 40],'Name','News_IRF_IV_other_vars');
    for ii = 2:length(varSelec) % Variable
        subplot(4,2,ii-1)
        plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(ii),1:Horizon,3,jj))*100,'LineWidth',2)
        hline(0,':k')
        hold on
        plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(ii),1:Horizon,1,jj))*100,'r--', ...
             'LineWidth',2)
        plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(ii),1:Horizon,5,jj))*100,'r--', ...
             'LineWidth',2)
        title(SVAR.i_var_str_names(:,varSelec(ii)),'FontSize',graph_opt.font_num,'FontWeight','bold')
        axis tight
        ylabel('percent','FontSize',12)
        xlabel('Horizon (quarters)','FontSize',12)
        set(gca,'XTick',0:4:H)
        set(gca,'LineWidth',2)
        hold off
        box off
    end    

end
