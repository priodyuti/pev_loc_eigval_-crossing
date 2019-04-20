%%
clear  

width = 6; %3;     % Width in inches
height = 6; %3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 10;      % Fontsize
lw = 1.5;      % LineWidth
msz = 6;       % MarkerSize

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

%%
% Each row of contains principal eigenvector entries of transition matrix at iteration i 
load('wheel_random_500_deloc.mat');

%size(data)
N = size(data,2);
t = size(data,1);
C = linspecer(N,'sequential');

for i=1:N 
 x = data(:,i);
 %semilogy(1:t,x,'color',C(i,:),'linewidth',2); 
 %plot(1:t,x,'color',C(i,:),'linewidth',2); 
 loglog(1:t,x,'color',C(i,:),'linewidth',2); 
 hold on 
end

xlabel('\tau_{evolution}','FontSize',16,'FontName', 'Times New Roman','FontWeight','bold','Color','k')
ylabel('{(x_1)_i}','FontSize',16,'FontName','Times New Roman','FontWeight','bold','Color','k');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','FontSize',16)

caxis([1 N]);
oldcmap = colormap(C);
%colormap(flipud(oldcmap));
colormap(oldcmap);
colorbar

%% Print setup
print('all_evec', '-dpng', '-r300');
print -depsc all_evec.eps    


