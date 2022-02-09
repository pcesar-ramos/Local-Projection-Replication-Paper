% Plot SP shock vs KPSS shock: 
% Here generate three surprise shocks:
% -- 1. Shock to Stock Prices. 
% -- 2. Shock to Consumer Confidence. 
% -- 3. Shock to KPSS measure. 
% Restricted = 0 (Variable ordered first, TFP second)
% Restricted = 1 (TFP ordered first, Variable ordered second)
%%
                           
%% Unrestricted:
restricted=0;
% ---------- SP shock: 
choose_vars =  [11 7  2 3 4 5 9  10 8  12];
[ data_SP,vmed_SP,vars_SP,datadraws_SP] = est_irf_fev(nirf,bound,ndraws,nconst, prior, begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);       
% ---------- CC shock: 
choose_vars =  [8 7  2 3 4 5 9  10  11  12]; 
[data_CC,vmed_CC,vars_CC,datadraws_CC] = est_irf_fev(nirf,bound,ndraws,nconst, prior, begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);
% ---------- KPSS shock: 
choose_vars =  [12 7  2 3 4 5 9  10 8  11]; 
[data_KPSS,vmed_KPSS,vars_KPSS,datadraws_KPSS] = est_irf_fev(nirf,bound,ndraws,nconst, prior, begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);

%% Plot Figures: 
fonttype          ='LM Roman 12';
ftsizeaxis        = 11;
titlefontsize     = 10;
[nvars,temp]=size(vars_KPSS);
N = nvars;
[nirf,temp,temp]=size(data_KPSS);
R=round(N/2)-1;

%%
% data_SP_0(:,[1 2 nvars+1 nvars+2],:)    = data_SP_0  (:,[2 1 nvars+2 nvars+1],:);
% data_CC_0(:,[1 2 nvars+1 nvars+2],:)    = data_CC_0  (:,[2 1 nvars+2 nvars+1],:);
% data_KPSS_0(:,[1 2 nvars+1 nvars+2],:)  = data_KPSS_0(:,[2 1 nvars+2 nvars+1],:);


%% Figure 9 of the paper!  TFP (1st row) +FEV (2nd row) only unrestricted. 
      
h9=figure('Color',[0.9412 0.9412 0.9412],'Position',[1 1 800-100 600-100],'Name','FIG9_CGV_RESTAT');   
figure(h9)
n=2;
subplot(2,3,1);                   
 
plot(1:nirf,data_KPSS(:,nvars+n,1), 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2); hold on                                         
plot(1:nirf,data_KPSS(:,nvars+n,2), 'r--','LineWidth',2); hold on                            
plot(1:nirf,data_KPSS(:,nvars+n,3),'r--','LineWidth',2); hold on
plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);   
ylabel('percent','FontSize',12)
xlabel('Horizon (quarters)','FontSize',12)
title('TFP to Patent-Based News Shock')                                
set(gca,'XTick',[0;4;8;12;16;20])
set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
set(gca, 'FontName', fonttype);
set(gca, 'FontSize', ftsizeaxis);
set(gca,'Layer','top');
alpha(.1)                                   
axis tight
ylim([-0.25, 0.5])
hold off;

subplot(2,3,2);           
grpyat = [(1:nirf)' data_SP(1:nirf,nvars+n,2); (nirf:-1:1)' data_SP(nirf:-1:1,nvars+n,3)];
patch(grpyat(:,1),grpyat(:,2),[0.4 0.3 0.2],'edgecolor','none'); hold on 
plot(1:nirf,data_SP(:,nvars+n,1), '-*', 'color', [0.4 0.3 0.2],'LineWidth',1.5); hold on
plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);   
ylabel('percent','FontSize',12)
xlabel('Horizon (quarters)','FontSize',12)
title('TFP to SP Shock')
set(gca,'XTick',[0;4;8;12;16;20])
set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
set(gca, 'FontName', fonttype);
set(gca, 'FontSize', ftsizeaxis);
set(gca,'Layer','top');
alpha(.1)                                 
axis tight
ylim([-0.25, 0.5])
hold off;

subplot(2,3,3);           
grpyat = [(1:nirf)' data_CC(1:nirf,nvars+n,2); (nirf:-1:1)' data_CC(nirf:-1:1,nvars+n,3)];
patch(grpyat(:,1),grpyat(:,2),[0.8 0.7 0.1],'edgecolor','none');  hold on 
plot(1:nirf,data_CC(:,nvars+n,1), '-o', 'color', [0.8 0.7 0.1],'LineWidth',1.5); hold on
plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);   
ylabel('percent','FontSize',12)
xlabel('Horizon (quarters)','FontSize',12)
title('TFP to Consumer-Confidence Shock')
set(gca,'XTick',[0;4;8;12;16;20])
set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
set(gca, 'FontName', fonttype);
set(gca, 'FontSize', ftsizeaxis);
set(gca,'Layer','top');
alpha(.1)                                  
axis tight
ylim([-0.25, 0.5])
hold off;

subplot(2,3,4);           
plot(1:nirf,data_KPSS(:,n,1),'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2); hold on                                         
plot(1:nirf,data_KPSS(:,n,2),'r--','LineWidth',2); hold on                            
plot(1:nirf,data_KPSS(:,n,3),'r--','LineWidth',2); hold on                 
plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);   
ylabel('Share of FEV','FontSize',12)
xlabel('Horizon (quarters)','FontSize',12)
title('TFP FEV share due to Patent-Based News Shock')                                  
set(gca,'XTick',[0;4;8;12;16;20])
set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
set(gca, 'FontName', fonttype);
set(gca, 'FontSize', ftsizeaxis);
set(gca,'Layer','top');
alpha(.1)                                   
axis tight
ylim ([0 25])
hold off;

subplot(2,3,5);                                     
grpyat = [(1:nirf)' data_SP(1:nirf,n,2); (nirf:-1:1)' data_SP(nirf:-1:1,n,3)];
patch(grpyat(:,1),grpyat(:,2),[0.4 0.3 0.2],'edgecolor','none'); 
hold on 
plot(1:nirf,data_SP(:,n,1), '-*', 'color', [0.4 0.3 0.2],'LineWidth',1.5);%, 'HandleVisibility', 'off' );
hold on
plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);   
ylabel('Share of FEV','FontSize',12)
xlabel('Horizon (quarters)','FontSize',12)
title(' TFP FEV share due to SP Shock')
set(gca,'XTick',[0;4;8;12;16;20])
set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
set(gca, 'FontName', fonttype);
set(gca, 'FontSize', ftsizeaxis);
set(gca,'Layer','top');
alpha(.1)                                 
axis tight
ylim ([0 25])
hold off;

subplot(2,3,6);                   
grpyat = [(1:nirf)' data_CC(1:nirf,n,2); (nirf:-1:1)' data_CC(nirf:-1:1,n,3)];
patch(grpyat(:,1),grpyat(:,2),[0.8 0.7 0.1],'edgecolor','none'); 
hold on 
plot(1:nirf,data_CC(:,n,1), '-o', 'color', [0.8 0.7 0.1],'LineWidth',1.5);%, 'HandleVisibility', 'off' );
hold on
plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);   
ylabel('Share of FEV','FontSize',12)
xlabel('Horizon (quarters)','FontSize',12)
title('TFP FEV share due to CC shock')
set(gca,'XTick',[0;4;8;12;16;20])
set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
set(gca, 'FontName', fonttype);
set(gca, 'FontSize', ftsizeaxis);
set(gca,'Layer','top');
alpha(.1)
                                    
axis tight
ylim ([0 25])
hold off;


%% For the historgram 
datadraws_KPSS  = datadraws_KPSS*100;
datadraws_SP    = datadraws_SP*100;
datadraws_CC    = datadraws_CC*100;

%% Figure 10 in the paper
h_aux=162;
t1_aux = (1:6:h_aux)';
v1=figure('Color',[0.9412 0.9412 0.9412],'Position',[1 1 800-100 600-100],'Name','FIG10_CGV_RESTAT');   
figure(v1);

        temp1 = datadraws_KPSS(1,nvars+n,:);
        temp2=  datadraws_SP(1,nvars+n,:);
        temp3 = datadraws_CC(1,nvars+n,:);


        h(1) = histogram(temp1,'facecolor','none','EdgeColor', 'r', 'LineStyle', '--');
        hold on
        h(2)= histogram(temp2,'facecolor', [0.4 0.3 0.2],'EdgeColor', [0.4 0.3 0.2],'LineStyle', '-.' );
        hold on
        h(3)= histogram(temp3,'facecolor', [0.8 0.7 0.1],'EdgeColor', [0.8 0.7 0.1],'LineStyle', ':' );
        hold on 

        h1_aux = repmat(median(temp1), h_aux/6,1);
        h2_aux = repmat(median(temp2), h_aux/6,1);
        h3_aux = repmat(median(temp3), h_aux/6,1);


       
      h(4)= plot(h1_aux,t1_aux, 'b', 'LineWidth', 2);
      hold on 
      h(5)=  plot(h2_aux,t1_aux, '-*', 'LineWidth', 2 , 'Color', [0.4 0.3 0.2]);
      hold on 
      h(6)= plot(h3_aux,t1_aux, '-o', 'LineWidth', 2, 'Color', [0.8 0.7 0.1]);
      hold off
       
  
         hlegend=legend(h([4 5 6 1 2 3]),'Median Impact Response (Patent Shock)', 'Median Imact Response (SP Shock)', 'Median Impact Response (CC Shock)', ...
     'Impact Distribution (Patent Shock)', 'Impact Distribution (SP Shock)', 'Impact Distribution (CC Shock)');
        set(gca, 'FontSize', ftsizeaxis);
        set(gca,'Layer','top');
        alpha(.1)
        axis tight
        hold off;

        
        %% Figure11:
        
   h1=figure('Color',[0.9412 0.9412 0.9412],'Position',[1 1 800-100 600-100],'Name','FIG11_CGV_RESTAT');   
    figure(h1)
        for n=3:8
            subplot(R-1,2,n-2);                                           
                                      hold on 
                                      h(1) = plot(1:nirf,data_KPSS(:,nvars+n,1), 'b-', 'LineWidth', 2);
                                      hold on                                         
                                      h(2) = plot(1:nirf,data_KPSS(:,nvars+n,2),'r--','LineWidth',2);%, 'HandleVisibility','off');
                                      hold on                            
                                      h(3) = plot(1:nirf,data_KPSS(:,nvars+n,3),'r--','LineWidth',2);% 'HandleVisibility','off');
                                      hold on
                          
                                      grpyat = [(1:nirf)' data_SP(1:nirf,nvars+n,2); (nirf:-1:1)' data_SP(nirf:-1:1,nvars+n,3)];   
                                        h(4) = patch(grpyat(:,1),grpyat(:,2),[0.4 0.3 0.2],'edgecolor','none');%, 'HandleVisibility','off'); 
                                        h(5) =plot(1:nirf,data_SP(:,nvars+n,1), '-*', 'color', [0.4 0.3 0.2],'LineWidth',1.5);%, 'HandleVisibility', 'off' );
                                       hold on
                                      
                                            %'-.'             
                                       
                                       
                                      grpyat = [(1:nirf)' data_CC(1:nirf,nvars+n,2); (nirf:-1:1)' data_CC(nirf:-1:1,nvars+n,3)];
                                       h(6) = patch(grpyat(:,1),grpyat(:,2),[0.8, 0.7, 0.1],'edgecolor','none');%, 'HandleVisibility','off'); 
                                      h(7) =plot(1:nirf,data_CC(:,nvars+n,1), '-o', 'color', [0.8, 0.7, 0.1],'LineWidth',1.5);%, 'HandleVisibility', 'off' );

                                      hold on
                                      
                                     plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);                                      
                                     ylabel('percent','FontSize',12)
                                     xlabel('Horizon (quarters)','FontSize',12)
                                     title(num2str(vars_KPSS(n,:)),'Interpreter','tex','FontSize',titlefontsize) %'FontWeight','bold',
                                    set(gca,'XTick',[0;4;8;12;16;20])
                                    set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
                                    %set(gca,'YTick',[-1.2;-0.6;0;0.6;1.2])
                                    set(gca, 'FontName', fonttype);
                                    set(gca, 'FontSize', ftsizeaxis);
                                    set(gca,'Layer','top');
                                    alpha(.1)
                                
                                    axis tight
                                    hold off;
        end
       set(legend, 'NumColumns' ,3)
       legend('Orientation','horizontal');%,'Location','southoutside')
       
 %hlegend=legend(h([1 2 5 4 7 6]),'Patent-based News Shock', 'Confidence bands', 'SP shock', 'Confidence bands', 'CC shock', 'Confidence bands' );
 hlegend=legend(h([1 5 7 2 4 6]),'Patent-based News Shock', 'SP shock',  'CC shock', 'Confidence bands', 'Confidence bands', 'Confidence bands');

 %%
clear data_KPSS data_CC data_SP datadraws_KPSS datadraws_CC datadraws_SP
