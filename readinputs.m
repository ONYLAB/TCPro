function  [antigenname,cohortname,Nrun,MeanFp,StdFp,SampleConcentration,cohort,SIcutoff,sigcutoff] = readinputs()

%Temporary
load KDeff.mat;

%Read options file
fileID = fopen('options.dat');
options = textscan(fileID,'%s %*[^\n]');
fclose(fileID);

antigenname = options{1,1}{1,1};
cohortname = options{1,1}{2,1};
cohorthlafilename = options{1,1}{3,1};
antigenfilename = options{1,1}{4,1};
Nrun = str2num(options{1,1}{5,1});
MeanFp = str2num(options{1,1}{6,1})/1e6;
StdFp = str2num(options{1,1}{7,1})/1e6;
SampleConcentration = (1/15)*str2num(options{1,1}{8,1})/1e6; %Scale the concentration so it reflects the effective concentration of a "15-mer" antigen
SIcutoff = str2num(options{1,1}{9,1});
sigcutoff = str2num(options{1,1}{10,1});

%Extract sequences
antigen = fastaread(antigenfilename);
numchains = length(antigen);

%Read cohort HLAs
cohort = readtable(cohorthlafilename,'Delimiter',',');

%Find unique alleles
uniqHLA = unique([cohort{:,2};cohort{:,3}]);

%Assign allele indexes to the table
for i = 1:length(uniqHLA)    
    
    cohort(strcmp(cohort{:,2},uniqHLA{i}),4) = {i};
    cohort(strcmp(cohort{:,3},uniqHLA{i}),5) = {i};
    KD = [];
    for j = 1:numchains
        peptidesequence = antigen(j).Sequence;
        KD = data(i);%[KD;IEDBScraper(peptidesequence,1,uniqHLA(i))]; %Temporary
    end
    KDeff = 1./sum(1./KD);
    cohort(strcmp(cohort{:,2},uniqHLA{i}),6) = {KDeff};
    cohort(strcmp(cohort{:,3},uniqHLA{i}),7) = {KDeff};
    
end
cohort(find(cohort{:,6}==cohort{:,7}),8)={1}; %If homozygous = 1

cohort.Properties.VariableNames(4) = {'Allele1index'};
cohort.Properties.VariableNames(5) = {'Allele2index'};
cohort.Properties.VariableNames(6) = {'KDeff1'};
cohort.Properties.VariableNames(7) = {'KDeff2'};
cohort.Properties.VariableNames(8) = {'Homozygous'};

cd ..\Output
% Write cohort member Effective KDs into output folder
writetable(cohort,[antigenname '_' cohortname '_cohortKDeff.csv'],'Delimiter',',');
cd ..\Input