function runalldonor()

cd ..\Input
[antigenname,cohortname,Nrun,MeanFp,StdFp,SampleConcentration,cohort,SIcutoff,sigcutoff] = readinputs();
cd ..\Bin

pdMS0 = makedist('Exponential','mu',5.2285); %Probability distribution for initial Maturation signal concentration in the well

randomseeds = randi(2^32-1,Nrun,height(cohort));
IncoBinary = zeros(Nrun,height(cohort));
ELISBinary = IncoBinary;
similarity = zeros(Nrun,1);

for traj = 1:Nrun
    
    for i = 1:height(cohort)
        
        theseed = randomseeds(traj,i);
        rng(theseed); %So that Fps MS0 cell type distriburions are consistently random per cohort member
        randomFp = random('Normal',MeanFp,StdFp); %always a fixed but random precursor frequency in cohort member i
        
        [Response{traj,i},pval{traj,i}] = Main_human(theseed,cohort(i,:),SampleConcentration,randomFp,pdMS0); 
        temp = Response{traj,i};
        temp2 = pval{traj,i};
        IncoBinary(traj,i) = sign(sum(((temp(1:4))>SIcutoff).*(temp2(1:4)<sigcutoff))); 
        ELISBinary(traj,i) = (temp(5)>SIcutoff).*(temp2(5)<sigcutoff); 
        
    end
    % Percent Similarity ELISpot vs Thymidine Incorporation
    similarity(traj,1) = sum(~xor(IncoBinary(traj,:)',ELISBinary(traj,:)'))/height(cohort); 
    
end

COMBOBinary = IncoBinary.*ELISBinary;

%% Calculate statistics
% Percent Responders
Inco = sum(IncoBinary,2)/height(cohort)*100;
ELIS = sum(ELISBinary,2)/height(cohort)*100;
COMBO = sum(COMBOBinary,2)/height(cohort)*100;

% Mean and 95% CI in percent responders
pd = fitdist(Inco,'Normal');
IncoAllStats = icdf(truncate(pd,0,100),[0.025 0.5 0.975]);
pd = fitdist(ELIS,'Normal');
ELISAllStats = icdf(truncate(pd,0,100),[0.025 0.5 0.975]);
pd = fitdist(COMBO,'Normal');
COMBOAllStats = icdf(truncate(pd,0,100),[0.025 0.5 0.975]);

% ELISpot vs Thymidine incorporation %similarity and 95% CI
pd = fitdist(similarity,'Normal');
SimilarityStats = icdf(truncate(pd,0,1),[0.025 0.5 0.975]);

%% Write output
cd ..\Output
%Write summary output 
fileID = fopen([antigenname '_' cohortname '_Summary.dat'],'w');
header = ['%Responders to ' antigenname ' in the cohort (' cohortname '):']; 
fprintf(fileID,'%s\n', header);
str1 = ['by Thymidine Incroporation (95% CI): %' num2str(IncoAllStats(2)) ' (' num2str(IncoAllStats(1)) ',' num2str(IncoAllStats(3)) ')'];
fprintf(fileID,'%s\n', str1);
str2 = ['by IL-2 ELISpot (95% CI): %' num2str(ELISAllStats(2)) ' (' num2str(ELISAllStats(1)) ',' num2str(ELISAllStats(3)) ')'];
fprintf(fileID,'%s\n', str2);
str3 = ['by Combined Thymidine Incroporation IL-2 ELISpot (95% CI): %' num2str(COMBOAllStats(2)) ' (' num2str(COMBOAllStats(1)) ',' num2str(COMBOAllStats(3)) ')'];
fprintf(fileID,'%s\n', str3);
str4 = '--';
fprintf(fileID,'%s\n', str4);
str5 = ['%Similarity between Thymidine Incroporation IL-2 ELISpot predictions (95% CI): %' num2str(SimilarityStats(2)) ' (' num2str(SimilarityStats(1)) ',' num2str(SimilarityStats(3)) ')'];
fprintf(fileID,'%s\n', str5);
fclose(fileID);

%Write Raw data
dlmwrite([antigenname '_' cohortname '_IncorporationRaw.dat'],IncoBinary,'delimiter',',');
dlmwrite([antigenname '_' cohortname '_ELISpotRaw.dat'],ELISBinary,'delimiter',',');
dlmwrite([antigenname '_' cohortname '_COMBORaw.dat'],COMBOBinary,'delimiter',',');
cd ..\Bin

delete Parameters.mat