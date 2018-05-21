function pdMS0 = runalldonor()

% SPW = [];
% 
% % ELISpotResponse(n,SimType+1)
% for i = 1:1000
%     disp(i)
%     temp = Main_human(i,1,0.3,0.1); %#ok<AGROW>
%     SPW = [SPW; temp];
% end

load SPWvsMS0.mat;
MS0 = SPW(:,2); 
SPW = round(SPW(:,1));

load onlymediumSPW.mat
histfit([AbzenaImedium;AbzenaIImedium],round(sqrt(597)),'exponential');
pd = fitdist([AbzenaImedium;AbzenaIImedium],'exponential');
% truncate(pd,0,inf)

% text([num2str(1/pd.mean,'%.2f') 'exp^' num2str(1/pd.mean,'%.2f')]);

xlabel('SPW (Medium Only)')
ylabel('Frequency');
AxFS = 24; %Ax Fontsize
AxLW = 2; %Ax LineWidth
set(gcf,'color','w');
set(gca,'fontsize', AxFS);
set(gca,'LineWidth',AxLW);
axis square

figure
scatter(MS0,SPW,'k.')
ylabel('Simulated SPW (Medium Only)')
xlabel('MS_0 (ng/L)');
AxFS = 24; %Ax Fontsize
AxLW = 2; %Ax LineWidth
set(gcf,'color','w');
set(gca,'fontsize', AxFS);
set(gca,'LineWidth',AxLW);
axis square
hold on
plot(0:25,3.113*(0:25),'r-','LineWidth',3)

figure
histfit([AbzenaImedium;AbzenaIImedium]/3.113,round(sqrt(597)),'exponential');
pdMS0 = fitdist([AbzenaImedium;AbzenaIImedium]/3.113,'exponential');
xlabel('Fitted MS_0 (ng/L)');
ylabel('Frequency');
AxFS = 24; %Ax Fontsize
AxLW = 2; %Ax LineWidth
set(gcf,'color','w');
set(gca,'fontsize', AxFS);
set(gca,'LineWidth',AxLW);
axis square
