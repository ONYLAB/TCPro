close all
load experiment.mat

AxFS = 20; %Ax Fontsize
AxLW = 2; %Ax LineWidth
angle = 90;

list = [2 4 6 7 9 13 14 15 16 17 10 18 19];
tab = 0.1;

IncoAll = 0;
ELISAll = 0;
COMBOAll = 0;
N = 0;

for i = 1:traj-1

    load([num2str(i) '_matlabRcutoff100Prol0597.mat'],'IncoBinary', 'ELISBinary','COMBOBinary');    
    
    IncoAll = IncoAll + sum(IncoBinary(:,1:50)')' + sum(IncoBinary(:,[51:54 56:100])')';
    ELISAll = ELISAll + sum(ELISBinary(:,1:50)')' + sum(ELISBinary(:,[51:54 56:100])')';
    COMBOAll = COMBOAll + sum(COMBOBinary(:,1:50)')' + sum(COMBOBinary(:,[51:54 56:100])')';
    
    N = N + [ones(10,1)*50;ones(9,1)*49];

end

% % % Merge 9 and 12
IncoAll(9) = IncoAll(9) + IncoAll(12);
ELISAll(9) = ELISAll(9) + ELISAll(12);
COMBOAll(9) = COMBOAll(9) + COMBOAll(12);
N(9) = N(9) + N(12);
experiment(9,:) = 0.5*(experiment(9,:) + experiment(12,:));
experiment = experiment(list,:);
% % % % % 

LBIncoAll = 100*(IncoAll(1:19)-   binoinv(0.05,N,IncoAll(1:19)./N)    )./N;
UBIncoAll = 100*(   binoinv(0.95,N,IncoAll(1:19)./N)-IncoAll(1:19)    )./N;

LBELISAll = 100*(ELISAll(1:19)-   binoinv(0.05,N,ELISAll(1:19)./N)    )./N;
UBELISAll = 100*(   binoinv(0.95,N,ELISAll(1:19)./N)-ELISAll(1:19)    )./N;

LBCOMBOAll = 100*(COMBOAll(1:19)-  binoinv(0.05,N,COMBOAll(1:19)./N)   )./N;
UBCOMBOAll = 100*(  binoinv(0.95,N,COMBOAll(1:19)./N)-COMBOAll(1:19)   )./N;

IncoAll = 100*IncoAll(1:19)./N;
ELISAll = 100*ELISAll(1:19)./N;
COMBOAll = 100*COMBOAll(1:19)./N;

IncoAll = IncoAll(list,:);
ELISAll = ELISAll(list,:);
COMBOAll = COMBOAll(list,:);

LBIncoAll = LBIncoAll(list,:);
LBELISAll = LBELISAll(list,:);
LBCOMBOAll = LBCOMBOAll(list,:);

UBIncoAll = UBIncoAll(list,:);
UBELISAll = UBELISAll(list,:);
UBCOMBOAll = UBCOMBOAll(list,:);

% % % % % % % % % % % % 
subplot(1,3,1)
h1 = bar((1:13)-tab,experiment(:,1),0.5);
hold on
h2 = bar((1:13)+tab,IncoAll,0.5);
hold on
errorbar((1:13)+tab,IncoAll,LBIncoAll,UBIncoAll,'s','MarkerFaceColor','auto','LineWidth',1)
title('Proliferation')
ylabel('%Responding Donors')
xticks(1:13);
xticklabels({'Pep2','Pep4', 'Pep6' , 'Pep7', 'Pep9', 'Pep13', 'Pep14', 'Pep15', 'Pep16', ...
    'Pep17', 'HuA33','Exenatide', 'KLH'});
xtickangle(angle)
axis square
legend boxoff
set(gcf,'color','w');
set(gca,'fontsize', AxFS);
set(gca,'LineWidth',AxLW);
ylim([0 100])
box off
plot([10.5 10.5],[0 100],'--')
legend([h1 h2],{'Assay','TCPro'},'Location','NorthWest')

subplot(1,3,2)
bar((1:13)-tab,experiment(:,2),0.5)
hold on
bar((1:13)+tab,ELISAll,0.5)
hold on
errorbar((1:13)+tab,ELISAll,LBELISAll,UBELISAll,'s','MarkerFaceColor','auto','LineWidth',1)
title('ELISpot')
% ylabel('%Responders')
xticks(1:13);
xticklabels({'Pep2','Pep4', 'Pep6' , 'Pep7', 'Pep9', 'Pep13', 'Pep14', 'Pep15', 'Pep16', ...
    'Pep17', 'HuA33','Exenatide', 'KLH'});
xtickangle(angle)
axis square
set(gcf,'color','w');
set(gca,'fontsize', AxFS);
set(gca,'LineWidth',AxLW);
ylim([0 100])
box off
plot([10.5 10.5],[0 100],'--')

subplot(1,3,3)
bar((1:13)-tab,experiment(:,3),0.5)
hold on
bar((1:13)+tab,COMBOAll,0.5)
hold on
errorbar((1:13)+tab,COMBOAll,LBCOMBOAll,UBCOMBOAll,'s','MarkerFaceColor','auto','LineWidth',1)
title('Proliferation+ELISpot')
% ylabel('%Responders')
xticks(1:13);
xticklabels({'Pep2','Pep4', 'Pep6' , 'Pep7', 'Pep9', 'Pep13', 'Pep14', 'Pep15', 'Pep16', ...
    'Pep17', 'HuA33', 'Exenatide', 'KLH'});
xtickangle(angle)
axis square
set(gcf,'color','w');
set(gca,'fontsize', AxFS);
set(gca,'LineWidth',AxLW);
ylim([0 100])
box off
plot([10.5 10.5],[0 100],'--')
