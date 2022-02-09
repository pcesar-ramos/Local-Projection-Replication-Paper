%'% '    GZ Credit Spread          '; %20
  % ' Excess Bond Premium          '; %21
  %' BAA-AAA Spread               '; %22
  
% (Change accordingly if you want BAA-AA  (22) and for GZ Credit spread
% (20)
choose_vars       = [12 7 2 3 4 5 9 10 8 11 21];
[data_credit,vmed_credit,vars_credit,datadraws_credit,y_credit] = est_irf_fev(nirf,bound,ndraws,nconst, prior, begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);

choose_vars       = [12 7 2 3 4 5 9 10 8 11];
[data_nocredit,vmed_nocredit,vars_nocredit] = est_irf_fev(nirf,bound,ndraws,nconst, prior, begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,DATA_MASTER,labels_vars);

%% ==============================================================================================================================================================================================
%----------------------------------------------------------------------
% Plot the responses of KPSS and TFP adjusted. 
%----------------------------------------------------------------------
fonttype                ='Arial';
ftsizeaxis              = 11;
titlefontsize           = 10;
[nvars_nocredit,temp]   =size(vars_nocredit);
[nvars_credit,temp]     =size(vars_credit);
N = nvars_nocredit;
[nirf,temp,temp]=size(data_credit);
%R=round(N/2);
R=round(N/2);

%%      
data_aux_nocredit               = zeros(size(data_credit));   
data_aux_nocredit(:,1:10,:)     = data_nocredit(:,1:10,:);
data_aux_nocredit(:,12:21,:)    = data_nocredit(:,11:20,:);

%% Figure B3 - comparing the two:
figb3=figure('Color',[0.9412 0.9412 0.9412],'Position',[1 1 800-100 600-100],'Name','FIGB3_CGV_REStat');   

 figure(figb3);
       for n=2:N+1
        subplot(R,2,n-1);                                       
                                    
             if n~= N+1
                h(1)= plot(1:nirf,data_aux_nocredit(:,nvars_credit+n,1), 'LineWidth',2); hold on 
                h(2)= plot(1:nirf,data_aux_nocredit(:,nvars_credit+n,2),'r--','LineWidth',2);hold on   
                h(3)= plot(1:nirf,data_aux_nocredit(:,nvars_credit+n,3),'r--','LineWidth',2);hold on
             end 
                h(4)=  plot(1:nirf,data_credit(:,nvars_credit+n,1), '-*','color', [0.4, 0.4, 0.4],'LineWidth',1.5);hold on 
                grpyat = [(1:nirf)' data_credit(1:nirf,nvars_credit+n,2); (nirf:-1:1)' data_credit(nirf:-1:1,nvars_credit+n,3)];
                h(5)=   patch(grpyat(:,1),grpyat(:,2),[0.5, 0.5, 0.5],'edgecolor','none'); 
                plot(1:nirf,zeros(nirf),':k','LineWidth',0.3);
                ylabel('percent','FontSize',12)
                xlabel('Horizon (quarters)','FontSize',12)
                title(num2str(vars_credit(n,:)),'Interpreter','tex','FontSize',titlefontsize) %'FontWeight','bold',
        set(gca,'XTick',[0;4;8;12;16;20])
        set(gca,'XTickLabel',['0 ';'4 ';'8 ';'12';'16';'20'])
        set(gca, 'FontName', fonttype);
        set(gca, 'FontSize', ftsizeaxis);
        set(gca,'Layer','top');
        alpha(.1)
        box off 
        axis tight
        hold off;                                
       end
      
     set(legend, 'NumColumns',2)
     legend('Orientation','horizontal');%,'Location','southoutside')
     hlegend=legend(h([1 4 2 5]), 'Benchmark',  'Controling for Credit Conditions', 'Confidence bands', 'Confidence bands');

     
%% Figure B2 of the Appendix:
x1     = y_credit(:,1);
x2     = y_credit(:,11);

% get HP component:
x1_cycle        = bandpass(x1,6,32);                % extract the cycle using BP filter. 
x1_cycle1       = x1 - hpfilter(x1,1600);           % extract the cycle using HP filter with LAMBDA=1600
std_x1          = std(x1_cycle)*100 ;               % calculate std (volatility) of output in the whole sample. 
% get HP component:
x2_cycle        = bandpass(x2,6,32);                % extract the cycle using BP filter. 
x2_cycle1       = x2 - hpfilter(x2,1600);           % extract the cycle using HP filter with LAMBDA=1600
std_x2          = std(x2_cycle)*100 ;               % calculate std (volatility) of output in the whole sample. 
corr_x1_x2_cycle = correl(x2_cycle, x1_cycle);
corr_x1_x2_raw  = correl(x2, x1);
t=1974:0.25:2010.5;

%%
figb2=figure('Color',[0.9412 0.9412 0.9412],'Position',[1 1 800-100 600-100],'Name','FIGB2_CGV_REStat');   
figure(figb2);
hold on 
yyaxis left
x  = t; 
y  = x1; 
xlim([1973 2010])
xlabel('r') 
ylabel('Patent-Based Innovation Index','Interpreter','tex', 'FontSize',10)
NBERbc(x, y, '-', 2, 'blue');
hold on 
yyaxis right
ylabel('Excess Bond Premium','Interpreter','tex', 'FontSize',10);
hold on 
y2 =  x2;
plot(x,y2, 'd-', 'LineWidth', 1.5);
hold off

clear data_aux_nocredit datadraws_credit data_credit 