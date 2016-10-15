function PlotExperiment( color )

global T3conv T4conv TSHconv
global time y

%% plot T4 Pools
    subplot(3,2,1*2-1, 'align');
    hold on
    plot(time/24, y(:,1)*T4conv,color,'LineWidth', 4); %plot q1 = T4
    hold off
for i = 2:3
    subplot(3,2,i*2-1, 'align');
    hold on
    plot(time/24, y(:,i)*T4conv,color,'LineWidth', 4); %plot q1 = T4
    hold off
end

%% plot T3 Pools
    subplot(3,2,(4-3)*2, 'align');
    hold on
    plot(time/24, y(:,4)*T3conv,color,'LineWidth', 4); %plot q4 = T3        
    hold off
for i = 5:6
    subplot(3,2,(i-3)*2, 'align');
    hold on
    plot(time/24, y(:,i)*T3conv,color,'LineWidth', 4); %plot q4 = T3        
    hold off
end

% plot1 = p(13) .* y(:,2) ./ (p(14) + y(:,2));
% plot2 = p(15) .* y(:,3) ./ (p(16) + y(:,3));

% if ( color ~= 'k' )
    subplot(3,2,4,'align');
    hold on
    plot(time/24,y(:,20)*400/3.2,'--g','LineWidth',4);
    % plot(time/24,plot1,'g','LineWidth',4);
    hold off
% end
% display(mean(plot1)*1000)
% figure
% plot(time/24,plot2,'g','LineWidth',4);
% display(mean(plot2)*1000)

% figure()
% plot(time/24,y(:,20),color,'LineWidth',2);
% figure()
% global newKprime newKdeg
% plot(time/24,newKprime * y(:,5) - newKdeg * y(:,20),color,'LineWidth',2);

%% plot TSH
% figure
% h3 = subplot(1,1,1, 'align');
% plot(time, y(:,7)*TSHconv,'LineWidth', 2); %plot q7 = TSH
% xlabel('hours'); ylabel('mU/l'); legend('TSH'); grid on;
% axis([0,max(time)*1.05,0,max(y(:,7)*TSHconv)*1.25]);