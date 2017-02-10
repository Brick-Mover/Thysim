function PlotSetup( )

global T3conv T4conv
global time y


poolName = {'Serum T4 Conc'; 'Liver+Kidney T4'; 'Muscle+Brain etc T4'; 'Serum T3 Conc'; 'Liver+Kidney T3 and D1'; 'Muscle+Brain etc T3'};

figure

% plot T4 Pools
i = 1;
    subplot(3,2,i*2-1, 'align');
    
    set(gca,'FontSize',16) %set the font size of everything, including the tick labels
    xlhand = get(gca,'xlabel'); %make a handle for the x axis label
%     xlabel('days') %label the x axis
    set(xlhand,'fontsize',25) %set the font size for the x axis label
    
    set(gca,'FontSize',16) %set the font size of everything, including the tick labels
    ylhand = get(gca,'ylabel'); %make a handle for the x axis label
    ylabel('ug/L') %label the x axis
    set(ylhand,'fontsize',25) %set the font size for the x axis label

    title(poolName(i),'FontSize',30);
    % xlabel('days','FontSize',16);
    % ylabel('nmol','FontSize',16);
    grid on;
    axis([0,max(time)*1.05/24,0,max(y(:,i))*1.25*T4conv]);
for i = 2:3
    subplot(3,2,i*2-1, 'align');
    
    set(gca,'FontSize',16) %set the font size of everything, including the tick labels
    xlhand = get(gca,'xlabel'); %make a handle for the x axis label
    if i == 3 
        xlabel('days') %label the x axis
    end
    set(xlhand,'fontsize',25) %set the font size for the x axis label
    
    set(gca,'FontSize',16) %set the font size of everything, including the tick labels
    ylhand = get(gca,'ylabel'); %make a handle for the x axis label
    ylabel('ug') %label the x axis
    set(ylhand,'fontsize',25) %set the font size for the x axis label
    
    title(poolName(i),'FontSize',30);
%     xlabel('days','FontSize',16);
%     ylabel('nmol','FontSize',16);
    grid on;
    axis([0,max(time)*1.05/24,0,max(y(:,i))*1.25*T4conv]);
end

%plot T3 Pools
i = 4;
    subplot(3,2,(i-3)*2, 'align');
    
    set(gca,'FontSize',16) %set the font size of everything, including the tick labels
    xlhand = get(gca,'xlabel'); %make a handle for the x axis label
%     xlabel('days') %label the x axis
    set(xlhand,'fontsize',25) %set the font size for the x axis label
    
    set(gca,'FontSize',16) %set the font size of everything, including the tick labels
    ylhand = get(gca,'ylabel'); %make a handle for the x axis label
    ylabel('ug/L') %label the x axis
    set(ylhand,'fontsize',25) %set the font size for the x axis label
    
    title(poolName(i),'FontSize',30);
%     xlabel('days','FontSize',16);
%     ylabel('nmol','FontSize',16);
    grid on;
    axis([0,max(time)*1.05/24,0,max(y(:,i))*1.25*T3conv]);
for i = 5:6
    subplot(3,2,(i-3)*2, 'align');
    
    set(gca,'FontSize',16) %set the font size of everything, including the tick labels
    xlhand = get(gca,'xlabel'); %make a handle for the x axis label
    if i == 6
        xlabel('days') %label the x axis
    end
    set(xlhand,'fontsize',25) %set the font size for the x axis label
    
    set(gca,'FontSize',16) %set the font size of everything, including the tick labels
    ylhand = get(gca,'ylabel'); %make a handle for the x axis label
    ylabel('ug') %label the x axis
    set(ylhand,'fontsize',25) %set the font size for the x axis label
    
    title(poolName(i),'FontSize',30);
%     xlabel('days','FontSize',16);
%     ylabel('nmol','FontSize',16);
    grid on;
    axis([0,max(time)*1.05/24,0,max(y(:,i))*1.25*T3conv]);
end