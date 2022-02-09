% Construct real quarterly patent-based innovation index:
% MV, November 2018
% ---------------------

clear all; close all; clc;
% Raw file used to produce these measures is: 
%/DATA/patents.csv, run STATA QuarterlyXI.do to obtain xiq.csv)

sumxiq       = csvread('DATA/xiq.csv', 1,1); % quarterly sjm of all xi's: deflated in 1982 (million) dollars using CPI (KPSS) 
GDPrq        = csvread('DATA/GDPC1-2.csv',1,1); %

sumxiq = sumxiq((1947-1926)*4+1:end,:); % cut to 1947q1-2010q4
GDPrq = GDPrq(1:256);
GDPrq = GDPrq*1000; % millions 

lxiq = log(sumxiq./GDPrq);

XX = xlsread('DATA/POPq');  %number in thousands
popq =  XX(:,3);       
popq = popq(1:252,:);  %from 1948Q1 to 2010Q4
popq = popq./1000000;

lxiq_pc = [NaN; NaN; NaN; NaN; log(sumxiq(5:end)./GDPrq(5:end)./popq)];


save patentindex.mat lxiq_pc 





% save xiq_share.mat 

% 
% 
% X0 = lxiq_pc;
% T0 = size(X0,1);
% 
% %Convert the dates to serial date numbers, as required by recessionplot. 
% dates = datenum([dates,ones(T0,2)]);
% 
% %Create a time series plot of the four credit default predictors. 
% figure;
% plot(dates,X0,'LineWidth',2);
% ax = gca;
% ax.XTick = dates(1:2:end);
% datetick('x','yyyy','keepticks')
% xlabel('Year');
% ylabel('Level');
% axis tight;

% recessionplot

