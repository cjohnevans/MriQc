function [RuleMatrix] = shewhart(qadata, qatitle, qaaxis, InvisibleFig, ...
    PlotDate, PlotLastN)
% function [RuleMatrix] = shewhart(qadata, qatitle, qaaxis, InvisibleFig, ...
%    PlotDate, PlotLastN)
%     
%   This version copied from ParkPlace archive on 8/6/20
% qadata = Vector data to be analysed
% qatitle = graph title
% qaaxis = plot axis
% InvisibleFig = 1; Don't display figure in matlab
% PlotDate = Plot against date (for EPIQA)
% PlotLastN = only plot the last N measurements (but analyse them all)
%

qadata_DM = qadata - mean(qadata);
qadata_SD = std(qadata);
qadata_CV=100*qadata_SD./mean(qadata);
qadata_normSD = qadata_DM ./ qadata_SD;

% Shewhart Warning Rules
% Rule 1 (Warning). Measurement > 2SD from mean
rule1 = (qadata_normSD>2) + (qadata_normSD<-2);
rule1plot=rule1;
rule1plot(rule1plot==0)=NaN; 
rule1_failed = find(rule1>0);

% Rule 2.  Measurement > 3SD from mean
rule2 = (qadata_normSD>3) + (qadata_normSD<-3);
rule2plot=rule2;
rule2plot(rule2plot==0)=NaN;
rule2_failed = find(rule2>0);

% Rule 3. Two consecutive measurements > 2SD. 
% copy and offset rule1, then sum.  Watch for 2's.
rule3 = [ rule1; 0] + [0; rule1];
rule3=rule3(1:(end-1));
rule3(rule3<2)=0;
rule3=rule3/2;
rule3plot=rule3;
rule3plot(rule3plot==0)=NaN;
rule3_failed = find(rule3>0);

% Rule4 Jump between two measurements is 4SD
rule4 = abs(diff(qadata_normSD));
rule4diff=rule4;
rule4(rule4>4)=-1; %set to -1 number if fails
rule4(rule4>-1)=0; % zero all passes.
rule4=abs(rule4); 
rule4 = [rule4; 0] + [0; rule4]; %catch both sides of the jump
rule4(rule4>0)=1;
rule4_failed = find(rule4>0);

% Rule 5: Four +/-1 SDs on the same side.
high1SD=(qadata_normSD>1);
low1SD=(qadata_normSD<-1);
for jj = 4:length(qadata_normSD)
    mvavg4hi(jj)=sum(high1SD(jj-3:jj));
    mvavg4lo(jj)=sum(low1SD(jj-3:jj));
end
mvavg4hi=mvavg4hi./4;
mvavg4hi(mvavg4hi<1)=0;
mvavg4lo=mvavg4lo./4;
mvavg4lo(mvavg4lo<1)=0;
mvavg4 = or(mvavg4hi, mvavg4lo);
% for jj = 4:length(qadata_normSD) % fail all four
%     if(mvavg4(jj)==1)
%         mvavg4((jj-3):jj)=1;
%     end
% end
rule5=mvavg4';
rule5_failed = find(rule5>0);

% Rule 6; 10 measures on same side of mean
mhigh=(qadata_normSD>0);
mlow=(qadata_normSD<0);
for jj = 10:length(qadata_normSD)
    mvavghi(jj)=sum(mhigh(jj-9:jj));
    mvavglo(jj)=sum(mlow(jj-9:jj));
end
mvavghi=mvavghi/10;
mvavghi(mvavghi<1)=0;
mvavglo=mvavglo/10;
mvavglo(mvavglo<1)=0;
rule6=or(mvavghi,mvavglo)';
rule6_failed = find(rule6>0);
% 
% figure(100);
% plot([1:236], qadata_normSD, 'k-o', ...
%     [1:236], mvavghi, '^r-', ...
%     [1:236], mvavglo, 'vb-', ...
%     [1:236], rule6, 'xk', ...
%     [1 236], [0 0], 'k-')
% 

warningplot = rule1plot;
failureplot = rule2 + rule3 + rule4 + rule5 + rule6;
failureplot(failureplot>1)=1;
failureplot(failureplot==0)=NaN;

%linxx = [1:length(qadata)];
linxx = qaaxis;
%figure
if(InvisibleFig) == 1
    set(gcf,'Visible','off');
end

plot(linxx, qadata_normSD,'b-', ...
    linxx, (qadata_normSD.*warningplot),'bo', ...
    linxx, (qadata_normSD.*failureplot),'rv', ...
    [linxx(1) linxx(end)], [1 1], 'k:', ...
    [linxx(1) linxx(end)], [-1 -1], 'k:', ...
    [linxx(1) linxx(end)], [2 2], 'k-', ...
    [linxx(1) linxx(end)], [-2 -2], 'k-', ...
    [linxx(1) linxx(end)], [3 3], 'k--', ...
    [linxx(1) linxx(end)], [-3 -3], 'k--', ...
    [linxx(1) linxx(end)], [0 0], 'k-')
title([ qatitle  ' Shewhart Chart (mean,CV%)']);
ylabel('Deviation from mean (SD)');
legend(['(' num2str(mean(qadata),'%.1f') ' ,'  num2str(qadata_CV,'%.1f') '%)'], ...
    'Warning', ...
    'Failure', ...
    'Location','SouthWest');

hold on 
plot(linxx, (qadata_normSD.*failureplot), 'v', 'Color','red', 'MarkerFaceColor', 'red', 'MarkerSize',7)
hold off

plotaxis=axis;

% if(PlotLastN)
%     plotaxis(1)= qaaxis(end-PlotLastN);
%     plotaxis(2)= qaaxis(end);
% else
%     plotaxis(1)= qaaxis(1);
%     plotaxis(2)= qaaxis(end);
% end

axis(plotaxis)


%NumTicks = 2;
%L = get(gca,'XLim');
%set(gca,'XTick',linspace(L(1),L(2),NumTicks))

if(PlotDate) 
    datetick('x','mm/yy');
end

%set(gca, 'XTickLabel',num2str(get(gca, 'XTick')))


RuleMatrix = [ rule1'; rule2'; rule3'; rule4'; rule5'; rule6'];
fail = sum(RuleMatrix,1);
idx = [1:length(qadata)];
RuleMatrix = [idx;  RuleMatrix]';

%return all, so I can do summary plot for EPIQA_shewhart
%for jj = length(qadata):-1:1
%    if fail(jj)==0
%        RuleMatrix(jj,:)=[];
%    end
%end

%figure(104);imagesc(RuleMatrix(:,2:end))

