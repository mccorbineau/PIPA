function [] = my_plot(filename)
load(filename,...
    'nRow','nCol','Xtrue','materials','X',...
    'obj_vec','snr_global_vec','snr_mats_vec','time_vec','method')
     
Xtrue_resh = reshape(Xtrue,size(Xtrue,1),nRow,nCol);
X_resh     = reshape(X,size(X,1),nRow,nCol);

%%%%%%%%%%%%%%%%%%% plot visual results 
for ii = 1:length(materials)
    figure
    set(gcf, 'units','centimeters','outerposition',[0 15 30 15]);

    % ground-truth
    subplot(121)
    imagesc(squeeze(Xtrue_resh(ii,:,:))); axis off
    colorbar('Location','southoutside');
    title('Ground-truth','interpreter','latex','fontsize',18)

    % reconstruction
    subplot(122)
    imagesc(squeeze(X_resh(ii,:,:))); axis off
    colorbar('Location','southoutside');
    title(strcat('Reconstruction:',{' '},num2str(snr_mats_vec(ii,end),'%.2f'),' dB'),'interpreter','latex','fontsize',18)
    
    sgtitle(strcat('\textbf{Visual results for',{' '},materials{ii},...
        {' '},'with',{' '},method,{' '},'after',{' '},...
        num2str(time_vec(end),'%.0f'),' seconds}'),'interpreter','latex','fontsize',20)
end

%%%%%%%%%%%%%%%%%%% plot SNR and objective function 
if strcmp(method,'PIPA')
    style='-r';
elseif strcmp(method,'PIPA-VM')
    style='-b';
end
figure
set(gcf, 'units','centimeters','outerposition',[0 0 30 15]);
        
% objective function
subplot(121)
plot(time_vec,obj_vec,style,'Linewidth',2)
ylabel({'$f(x)+g(x)$'},'Interpreter','latex','fontsize',18)
xlabel({'Time~(s)'},'Interpreter','latex','fontsize',18)
set(gca,'Fontsize',18,'TickLabelInterpreter','latex')
title('Objective function','Interpreter','latex','fontsize',18)
        
% signal-to-noise ratio
subplot(122)
plot(time_vec,snr_global_vec,style,'Linewidth',2)
ylabel({'SNR'},'Interpreter','latex','fontsize',18)
xlabel({'Time~(s)'},'Interpreter','latex','fontsize',18)
set(gca,'Fontsize',18,'TickLabelInterpreter','latex')
xlim([0,time_vec(end)])
ylim([0.9*min(snr_global_vec),1.05*max(snr_global_vec)])
title('Signal-to-noise ratio','Interpreter','latex','fontsize',18)
    
sgtitle(strcat('\textbf{',method,{' '},'results}'),'interpreter','latex','fontsize',20)
end