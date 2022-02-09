% esty1 = 1961;     
% estq1 = 1;        
% esty2 = 2010;       
% estq2 = 4;        

%%
choose_vars  = [15 7  2 3 4 5 9 10 8 11];
[data_manufacturing,vmed_manufacturing,vars] = est_irf_fev(nirf,bound,ndraws,nconst, prior, begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);
%% Run industry. 
choose_vars  = [23 7  2 3 4 5 9 10 8 11];
[data_industry17,vmed_industry17,vars] = est_irf_fev(nirf,bound,ndraws,nconst, prior, begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);
%% Run Services
choose_vars  = [17 7  2 3 4 5 9 10 8 11];
[data_services,vmed_services,vars] = est_irf_fev(nirf,bound,ndraws,nconst, prior, begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);
%% Run Business Services
choose_vars  = [24 7  2 3 4 5 9 10 8 11];
[data_bservices,vmed_bservices,vars] = est_irf_fev(nirf,bound,ndraws,nconst, prior, begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);
%% ==============================================================================================================================================================================================
%% Figures:
fonttype          ='Arial';
ftsizeaxis        = 11;
titlefontsize     = 10;
[nvars,temp]=size(vars);
N = nvars;
[nirf,temp,temp]=size(data_manufacturing);
R=round(N/2);
%%

%% FIGURE 7:
fig7=figure('Color',[0.9412 0.9412 0.9412],'Position',[1 1 800-100 600-100],'Name','FIG7_CGV_RESTAT');   


figure(fig7)

subplot(1,3,1);           
    plot(1:nirf,data_manufacturing(:,nvars+2,1), 'LineWidth',2); hold on                                        
    plot(1:nirf,data_manufacturing(:,nvars+2,2),'r--','LineWidth',2);hold on                            
    plot(1:nirf,data_manufacturing(:,nvars+2,3),'r--','LineWidth',2); hold on
    plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);
    ylabel('percent','FontSize',12)
    xlabel('Horizon (quarters)','FontSize',12)
    title('TFP to Manufacturing News', 'FontSize',titlefontsize) %'FontWeight','bold',
    set(gca,'XTick',[0;4;8;12;16;20])
    set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
    ylim([-0.1 0.6])
    set(gca, 'FontName', fonttype);
    set(gca, 'FontSize', ftsizeaxis);
    set(gca,'Layer','top');
    axis tight  
    hold on 
 
subplot(1,3,2);            
    plot(1:nirf,data_industry17(:,nvars+2,1), 'LineWidth',2); hold on                                        
    plot(1:nirf,data_industry17(:,nvars+2,2),'r--','LineWidth',2); hold on                            
    plot(1:nirf,data_industry17(:,nvars+2,3),'r--','LineWidth',2); hold on
    plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);
    ylabel('percent','FontSize',12)
    xlabel('Horizon (quarters)','FontSize',12)
    title('TFP to Electronic/Electrical News','FontSize',titlefontsize) %
    set(gca,'XTick',[0;4;8;12;16;20])
    set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
    ylim([-0.1 0.6])
    set(gca, 'FontName', fonttype);
    set(gca, 'FontSize', ftsizeaxis);
    set(gca,'Layer','top');
    axis tight
    hold on

subplot(1,3,3);            
    plot(1:nirf,data_bservices(:,nvars+2,1), 'LineWidth',2); hold on                                        
    plot(1:nirf,data_bservices(:,nvars+2,2),'r--','LineWidth',2);hold on                            
    plot(1:nirf,data_bservices(:,nvars+2,3),'r--','LineWidth',2); hold on
    plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);
    ylabel('percent','FontSize',12)
    xlabel('Horizon (quarters)','FontSize',12)
    title('TFP to Business Services News','FontSize',titlefontsize) 
    set(gca,'XTick',[0;4;8;12;16;20])
    set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
    ylim([-0.1 0.6])
    set(gca, 'FontName', fonttype);
    set(gca, 'FontSize', ftsizeaxis);
    set(gca,'Layer','top');
    axis tight
    hold off
 
 
%% Figure C1 Appendix:

figc1=figure('Color',[0.9412 0.9412 0.9412],'Position',[1 1 800-100 600-100],'Name','FIGC1_CGV_RESTAT');   

figure(figc1);
for n=3:N
    subplot((R-1)*2,2,(n-2)+(n-3));            
    plot(1:nirf,data_manufacturing(:,nvars+n,1), 'LineWidth',2); hold on                                        
    plot(1:nirf,data_manufacturing(:,nvars+n,2),'r--','LineWidth',2); hold on                            
    plot(1:nirf,data_manufacturing(:,nvars+n,3),'r--','LineWidth',2); hold on
    plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);
    ylabel('percent','FontSize',12)
    xlabel('Horizon (quarters)','FontSize',12)
    title(num2str(vars(n,:)),'Interpreter','tex','FontSize',titlefontsize)       
    set(gca,'XTick',[0;4;8;12;16;20])
    set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
    set(gca, 'FontName', fonttype);
    set(gca, 'FontSize', ftsizeaxis);
    set(gca,'Layer','top');
    axis tight  
end
hold on 

for n=3:N
    subplot((R-1)*2,2,(n-2)+(n-2));            

    plot(1:nirf,data_industry17(:,nvars+n,1), 'LineWidth',2); hold on                                         
    plot(1:nirf,data_industry17(:,nvars+n,2),'r--','LineWidth',2);hold on                            
    plot(1:nirf,data_industry17(:,nvars+n,3),'r--','LineWidth',2);hold on
    plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);
    ylabel('percent','FontSize',12)
    xlabel('Horizon (quarters)','FontSize',12)
    title(num2str(vars(n,:)),'Interpreter','tex','FontSize',titlefontsize) 
    set(gca,'XTick',[0;4;8;12;16;20])
    set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
    set(gca, 'FontName', fonttype);
    set(gca, 'FontSize', ftsizeaxis);
    set(gca,'Layer','top');
    axis tight
    hold off
    
end
%%


%% TABLE 2:
disp('TABLE 2: ==============================================================')
hor = [1 5 9 17 21]';
for j=1:5;
    a= hor(j);
    for i=1:nvars
        FEV_table_manufacturing(j,:,i) = [data_manufacturing(a,i,2), data_manufacturing(a,i,1),data_manufacturing(a,i,3)];                   
    end
end
hor1 = [0 4 8 16 20]';
hor1 = round(hor1,0);
VarNames = {'horizon', 'lb', 'median', 'ub'};

%for i=1:nvars;
 i=2;
 disp(['Manufacturing:', num2str(vars(i,:))]) %Interpreter','tex','FontSize',titlefontsize))
 aa = [hor1, FEV_table_manufacturing(:,:,i)];
 aa = round(aa,1);
 disp(table(aa(:,1),aa(:,2),aa(:,3),aa(:,4), 'VariableNames',VarNames))
 %end
% -------------- Electronic/Electrical:
for j=1:5;
    a= hor(j);
    for i=1:nvars
        FEV_table_industry17(j,:,i) = [data_industry17(a,i,2), data_industry17(a,i,1),data_industry17(a,i,3)];                   
    end
end
i=2;
 disp(['Electronic/Electrical:', num2str(vars(i,:))]) %Interpreter','tex','FontSize',titlefontsize))
 aa = [hor1, FEV_table_industry17(:,:,i)];
 aa = round(aa,1);
 disp(table(aa(:,1),aa(:,2),aa(:,3),aa(:,4), 'VariableNames',VarNames))
 %end
 

% ------------  Business Services 
for j=1:5;
    a= hor(j);
    for i=1:nvars
        FEV_table_bservices(j,:,i) = [data_bservices(a,i,2), data_bservices(a,i,1),data_bservices(a,i,3)];                   
    end
end

i=2;
 disp(['Business Services:', num2str(vars(i,:))]) %Interpreter','tex','FontSize',titlefontsize))
 aa = [hor1, FEV_table_bservices(:,:,i)];
 aa = round(aa,1);
 disp(table(aa(:,1),aa(:,2),aa(:,3),aa(:,4), 'VariableNames',VarNames))
 
 
%% TABLE C1:
disp('TABLE C1: ==============================================================')

 for i=3:6;
 disp(['Manufacturing:', num2str(vars(i,:))]) %Interpreter','tex','FontSize',titlefontsize))
 aa = [hor1, FEV_table_manufacturing(:,:,i)];
 aa = round(aa,1);
 disp(table(aa(:,1),aa(:,2),aa(:,3),aa(:,4), 'VariableNames',VarNames))
 end

 for i=3:6;
 disp(['Electronic/Electrical:', num2str(vars(i,:))]) %Interpreter','tex','FontSize',titlefontsize))
 aa = [hor1, FEV_table_industry17(:,:,i)];
 aa = round(aa,1);
 disp(table(aa(:,1),aa(:,2),aa(:,3),aa(:,4), 'VariableNames',VarNames))
 end
 
 for i=3:6;
 disp(['Business Services:', num2str(vars(i,:))]) %Interpreter','tex','FontSize',titlefontsize))
 aa = [hor1, FEV_table_bservices(:,:,i)];
 aa = round(aa,1);
 disp(table(aa(:,1),aa(:,2),aa(:,3),aa(:,4), 'VariableNames',VarNames))
 end


clear data_manufacturing data_bservices data_industry17
