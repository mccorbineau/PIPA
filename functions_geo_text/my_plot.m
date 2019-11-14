function [] = my_plot(filename)
%==========================================================================
% Plots visual results, the objective function, and the signal-to-noise
% ratio obtained with regards to time using PIPA-VM.
%
% Input -------------------------------------------------------------------
%    filename (str): file name to load the results from
%==========================================================================

load(filename,...
    'sample_name','time_max','x_true','y','chi','H','weightTV','time_init',...
    'obj_vec','snr_vec','time_vec','x','x_true')

[m,n]  = size(x_true);            % image size
N      = m*n;                     % number of pixels
xt     = reshape(x(1:N),m,n);     % texture
xg     = reshape(x(N+1:end),m,n); % geometry
my_snr = @(ut,ug) -10*log10(norm(ut(:)+ug(:)-x_true(:))^2/norm(x_true(:))^2);
       
%%%%%%%%%%%%%%%%%%% plot visual results 
figure
set(gcf, 'units','centimeters','outerposition',[0 15 60 15]);

% texture
subplot(141)
imagesc(xt); colormap(gray); axis off
colorbar('Location','southoutside');
title('Texture $x^t$','interpreter','latex','fontsize',18)
      
% geometry
subplot(142)
imagesc(xg); colormap(gray); axis off
colorbar('Location','southoutside');
title('Geometry $x^g$','interpreter','latex','fontsize',18)
      
% texture+geometry
subplot(143)
imagesc(xt+xg); colormap(gray); axis off
colorbar('Location','southoutside');
title(strcat('Reconstruction $x^{t+g}$:',{' '},num2str(my_snr(xt,xg),'%.2f'),' dB'),'interpreter','latex','fontsize',18)
      
% ground-truth
subplot(144)
imagesc(x_true); colormap(gray); axis off
colorbar('Location','southoutside');
title('Ground-truth','interpreter','latex','fontsize',18)
      
sgtitle(strcat('\textbf{Visual results with PIPA-VM after',{' '},num2str(time_vec(end),'%.0f'),' seconds}'),'interpreter','latex','fontsize',20)
      
%%%%%%%%%%%%%%%%%%% plot SNR and objective function 
figure
set(gcf, 'units','centimeters','outerposition',[0 0 30 15]);
        
% objective function
subplot(121)
plot(time_vec,obj_vec,'-b','Linewidth',2)
ylabel({'$f(x)+g(x)$'},'Interpreter','latex','fontsize',18)
xlabel({'Time~(s)'},'Interpreter','latex','fontsize',18)
set(gca,'Fontsize',18,'TickLabelInterpreter','latex')
title('Objective function','Interpreter','latex','fontsize',18)
        
% signal-to-noise ratio
subplot(122)
plot(time_vec,snr_vec,'-b','Linewidth',2)
ylabel({'SNR'},'Interpreter','latex','fontsize',18)
xlabel({'Time~(s)'},'Interpreter','latex','fontsize',18)
set(gca,'Fontsize',18,'TickLabelInterpreter','latex')
xlim([0,time_vec(end)])
ylim([0.9*min(snr_vec),1.05*max(snr_vec)])
title('Signal-to-noise ratio','Interpreter','latex','fontsize',18)
    
sgtitle('\textbf{PIPA-VM results}','interpreter','latex','fontsize',20)