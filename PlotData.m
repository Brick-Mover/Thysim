function PlotData( color )

global fitIndex;
global T3conv T4conv TSHconv
global dataTimeForPlots
global dataPointsT3Red errorT3Red
global dataPointsT3Blue errorT3Blue
global dataPointsT3Pooled errorT3Pooled
global dataPointsT4Red errorT4Red
global dataPointsT4Blue errorT4Blue
global dataPointsT4Pooled errorT4Pooled
global wenzelDataTime wenzelDataT3 wenzelDataT4

%% Plot T4 Plasma Error Bars
subplot(3,2,1,'align');
hold on
z=[0 0 0 0 0 0 0 0 0];
if fitIndex == 3
    if strcmp(color,'red') || strcmp(color,'')
        T3 = errorbar(dataTimeForPlots/24, dataPointsT4Pooled,errorT4Pooled,errorT4Pooled,'. red','markersize',30);
        set(T3,'linewidth',3)
    end
else
    if strcmp(color,'red') || strcmp(color,'')
        T4 = errorbar(dataTimeForPlots/24, dataPointsT4Red,errorT4Red,errorT4Red,'. red ','markersize',30);
        set(T4,'linewidth',3)
    end
    if strcmp(color,'blue') || strcmp(color,'')
        T1 = errorbar(dataTimeForPlots/24, dataPointsT4Blue,errorT4Blue,errorT4Blue,'. blue','markersize',30);
        set(T1,'linewidth',3)
    end
end
if strcmp(color,'wenzel') || strcmp(color,'')
    plot(wenzelDataTime/24, wenzelDataT4, '. k', 'markersize',15);
end
hold off

%% Plot T3 Plasma Error Bars
subplot(3,2,2,'align');
hold on
if fitIndex == 3
    if strcmp(color,'red') || strcmp(color,'')
        T5 = errorbar(dataTimeForPlots/24, dataPointsT3Pooled,errorT3Pooled,errorT3Pooled,'. red','markersize',30);
        set(T5,'linewidth',3)
    end
else
    if strcmp(color,'red') || strcmp(color,'')
        T6 = errorbar(dataTimeForPlots/24, dataPointsT3Red,errorT3Red,errorT3Red,'. red','markersize',30);
        set(T6,'linewidth',3)
    end
    if strcmp(color,'blue') || strcmp(color,'')
        T2 = errorbar(dataTimeForPlots/24, dataPointsT3Blue,errorT3Blue,errorT3Blue,'. blue','markersize',30);
        set(T2,'linewidth',3)
    end
end
if strcmp(color,'wenzel') || strcmp(color,'')
    plot(wenzelDataTime/24, wenzelDataT3, '. k', 'markersize',15);
end
hold off

%% plot TSH
% figure
% h3 = subplot(1,1,1, 'align');
% plot(time, y(:,7)*TSHconv,'LineWidth', 2); %plot q7 = TSH
% xlabel('hours'); ylabel('mU/l'); legend('TSH'); grid on;
% axis([0,max(time)*1.05,0,max(y(:,7)*TSHconv)*1.25]);