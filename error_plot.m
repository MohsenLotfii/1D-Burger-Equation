% Errorplot in t=40

% inviscous for 40 & 80 & 160 cells
% Nrm_Err_FOU = [0.3487 0.1507 0.0760];
% Nrm_Err_SOU = [0.1087 0.1636 0.1870];
% Nrm_Err_QUICK = [0.2009 0.1410 0.2227];
% cell_size = [0.513 0.253 0.126];

% viscous for 20 & 160 cells
Nrm_Err_FOU = [0.4371 0.1236];
Nrm_Err_SOU = [0.1402 0.1810];
Nrm_Err_QUICK = [0.2334 0.2227];
cell_size = [1.053 0.126];

figure(7)
loglog(cell_size , Nrm_Err_FOU ,'-b*', 'linewidth' , 3)
hold on
loglog(cell_size , Nrm_Err_SOU ,'-gs', 'linewidth' , 3)
hold on
loglog(cell_size , Nrm_Err_QUICK ,'-r^', 'linewidth' , 3)
legend('FOU Method','SOU Method','QUICK Method','Lax_friedrichs Method','Lax_Wendroff Method')
xlabel('Cell Size')
ylabel('Error Norm')
title('Error Norm vs Cell Sizes')
set(gca,'FontName','Times New Roman','FontSize',10,'fontWeight','bold');
grid on