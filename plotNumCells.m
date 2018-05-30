function plotNumCells(t_record,y_record,N)

NA = 6.0221367e23;
% % Get Initial conditions
% Ag1=y_record(numel(t_record),1);
% MS1=y_record(numel(t_record),2); %No new Maturation signal
% ID1=y_record(numel(t_record),3);
% MD1=y_record(numel(t_record),4);
% AgE1=y_record(numel(t_record),5);
% pE1=y_record(numel(t_record),(5+1):(5+N));
% ME1=y_record(numel(t_record),(6+N):(11+N));
% pME1=y_record(numel(t_record),(12+N):(11+7*N));
% pM1=y_record(numel(t_record),(12+7*N):(11+13*N));
% M1=y_record(numel(t_record),(12+13*N):(17+13*N));
% NT1=y_record(numel(t_record),18+13*N);
% AT_N1=y_record(numel(t_record),(19+13*N):(18+14*N));
% Prolif1=y_record(numel(t_record),(19+14*N):(19+15*N));
% 
% % Initial condition vector
% Prolif1 = 0.0*Prolif1; %This is a place holder for proliferation
% yic2 = [Ag1, MS1, ID1, MD1, AgE1/scaleVOL, pE1/scaleVOL,  ME1/scaleVOL, pME1/scaleVOL,  pM1/scaleVOL, M1/scaleVOL,NT1, AT_N1, Prolif1]';
% yic2 = yic2*scaleVOL;
% 
% % dosing interval
% t_start=DayLimit;
% t_end=DayLimit+IncubationTime;
% tspan2 = linspace(t_start, t_end, numberoftimesamples);
% 
% %pars, Vp changes
% pars(3)=Vp*scaleVOL;
% 
% % Call ODE
% [T2,Y2]=ode15s(@f, tspan2, yic2, options,pars);
% 
% % Record the result into new matrix
% t_record=[t_record;T2]; %#ok<NASGU>
% y_record=[y_record; Y2];
% 
% 
% 
% 
% Output all state variables
% Ag, Ag in the plasma compartment, pmole
Ag_vector=y_record(:,1); %#ok<NASGU>
% MS, Maturation signal for activation of immature dendritic cells
% (Here is LPS, ng)
MS_vector=y_record(:,2); %#ok<NASGU>
% ID,	Immature dendritic cells, number of cells
ID_vector=y_record(:,3); 
% MD,	Mature dendritic cells, number of cells
MD_vector=y_record(:,4); 

%AgE, Ag in the endosomes, pmole
AgE_vector=y_record(:,5); %#ok<NASGU>
% pE,	free epitope peptides from Ag digestion , pmole
pE_vector=y_record(:,(5+1):(5+N)); %#ok<NASGU>
% ME,	free MHC II molecule in endosome , pmole
ME_vector=y_record(:,(6+N):(11+N)); %#ok<NASGU>
% pME,	T-epitope-MHC-II complex in endosome, pmole
pME_vector=y_record(:,(12+N):(11+7*N)); %#ok<NASGU>
% pM,	T-epitope-MHC-II on dendritic cell membrane, pmole
pM_vector=y_record(:,(12+7*N):(11+13*N));
% M, free MHC II molecule on dendritic cell menbrane, pmole
M_vector=y_record(:,(12+13*N):(17+13*N)); %#ok<NASGU>
% NT,	Na�ve helper T cells, number of cells
NT_vector=y_record(:,18+13*N); 
% AT_N,	activated helper T cells derived from naive T cells, number of cells
AT_N_vector=y_record(:,(19+13*N):(18+14*N));
% Prolif,	activated helper T cells derived from memory T cells, number of cells
Prolif_vector=y_record(:,(19+14*N):(19+15*N));

% Post-processing calculation
% Antigen presentation processes
pM_NUMBER_M1=pM_vector(:,(1:N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 1
pM_NUMBER_M2=pM_vector(:,(N+1):(2*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 2
pM_NUMBER_M3=pM_vector(:,(2*N+1):(3*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 3
pM_NUMBER_M4=pM_vector(:,(3*N+1):(4*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 4
pM_NUMBER_M5=pM_vector(:,(4*N+1):(5*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 5
pM_NUMBER_M6=pM_vector(:,(5*N+1):(6*N))*1E-12*NA; % T-epitope-MHC II complex number associated with MHCII molecule 6

Total_pM(:,1:N)=pM_NUMBER_M1+pM_NUMBER_M2+pM_NUMBER_M3+pM_NUMBER_M4+pM_NUMBER_M5+pM_NUMBER_M6;   %#ok<NASGU> % Total number of T-epitope-MHC II complex

% % Save the results
% savefile='results.mat';
% 
% save(savefile, 'koff', 't_record', 'Ag_vector', 'MS_vector', 'ID_vector', 'MD_vector',...
%     'AgE_vector', 'pE_vector', 'ME_vector', 'pME_vector', 'pM_vector', 'M_vector', ...
%     'NT_vector', 'AT_N_vector', 'Prolif_vector', 'Total_pM');

% % % % % % % % % % % % % % % % %
figure
LT = 3; %Line thickness
AxFS = 24; %Ax Fontsize
AxLW = 2; %Ax LineWidth
xlabeltext = 'Time (days)';
ylabeltext = '#cells';
legenddata = {'Immature DC','Mature DC','Naive T','Activated T'};

plot(t_record,ID_vector(:,1),'LineWidth',LT);
hold on
plot(t_record,MD_vector(:,1),'LineWidth',LT);
plot(t_record,NT_vector(:,1),'LineWidth',LT);
plot(t_record,AT_N_vector(:,1:N),'LineWidth',LT);

set(gcf,'color','w');
set(gca,'fontsize', AxFS);
set(gca,'LineWidth',AxLW);
legend(legenddata,'Location','eastoutside');
xlabel(xlabeltext);
ylabel(ylabeltext);
set(gca,'yscale','log')
axis square
