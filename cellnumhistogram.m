% NumSamples = 1e4;
% MinNumPBMCs = 4e6; %per mL
% MaxNumPBMCs = 6e6; %per mL
% Vp = 1e-3;
% for i = 1:NumSamples
%     [DC(i,1),NK(i,1),BC(i,1),CD4(i,1),CD8(i,1),MC(i,1)] = collectdonorPBMC_analyze(MinNumPBMCs,MaxNumPBMCs,Vp);
%     disp(i)
% end

load EstimatedNumberofTCells.mat
histogram(CD4)
xlabel('Number of CD4+ Cells')
ylabel('Frequency (out of 10,000 Simulations)')
AxFS = 24; %Ax Fontsize
AxLW = 2; %Ax LineWidth
set(gcf,'color','w');
set(gca,'fontsize', AxFS);
set(gca,'LineWidth',AxLW);
title('{\mu}:2.3453 10^6{\pm}3.8443 10^5 Cells/mL')
axis square