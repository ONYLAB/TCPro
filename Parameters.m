function Parameters(SimType,ProteinLength,Va,seed,SampleConcentration,Fp,pdMS0)

%% Therapeutic protein dosage

%Va is the whole volume of the assay (sample + cells)
% SampleConcentration [Molar] This is the actual concentration of the samples with respect to the cell-excluded volume
Dose = SimType*SampleConcentration*Va*1e12;  % pmole

% LE = 0.5; %ng/L Lower bound for Endoxin
% UE = 25; %ng/L Upper bound
rng(seed);
Endotoxin = random(pdMS0); %LE + (UE-LE)*rand; %ng/L

% NA: Avogadro constant
NA = 6.0221367e23;

%% Celltype distribution
%Assay uncertainty
MinNumPBMCs = 4e6; %per ml
MaxNumPBMCs = 6e6; %per ml
VolProliferationCellStock = Va * 0.5; %Initial total cell volume is half of the total well volume
[ID0,NK,BC,NT0,~,MC] = collectdonorPBMC_analyze(MinNumPBMCs,MaxNumPBMCs,VolProliferationCellStock,seed);

%% T-epitope characteristics of therapeutic proteins
% N: the number of T-epitope (protein-specific)
% ME0:  initial amount of MHC-II molecule in a single mature dendritic cell pmole
% kon: on rate for for T-epitope-MHC-II binding
% koff:	off rate for for T-epitope-MHC-II binding

% Number of HLA DRB1 molecules
numHLADR = 1e5;
[kon,koff,N,ME0] = doNETMHCIIpan(SimType,ProteinLength,NA,numHLADR); %N: #Epitopes

%% Dendritic cells
% BetaID: death rate of immature dendritic cells.
BetaID=0.0924; % day-1

% DeltaID:	maximum activation rate of immature dendritic cells
DeltaID=1.66; % day-1

% Km,MS: LPS concentration at which immature DC activation rate is 50% maximum.
KMS= 9.852E3; % ng/L

% BetaMD:	death rate of mature dendritic cells
BetaMD=0.2310; % day-1

%% Antigen presentation

% AlphaAgE:  Ag internalization rate constant by mature dendritic cells
AlphaAgE=14.4; %day-1

%	BetaAgE	: degradation rate for AgE in acidic vesicles
BetaAgE = 17.28; % day-1

% Beta_p: degradation rate for T-epitope peptide (same for all peptides)
Beta_p= 14.4; % day-1

% BetaM: degradation rate for MHC-II
BetaM= 1.663; % day-1.

% Beta_pM: degradation rate for pME
Beta_pM= 0.1663; % day-1

% kext: exocytosis rate for peptide-MHC-II complex from endosome to cell menbrane
kext = 28.8;  % day-1

% kin: internalization rate for peptide-MHC-II complex on DC membrane
kin = 14.4;  % day-1

% KpM_N: number of peptide-MHCII to achieve 50% activation rate of na�ve helper T cells
KpM_N = 400;  % number of complex

% VD: volume of a dendritic cell
VD=2.54e-12; % L

% VE: volume of endosomes in a dendritic cell
VE=1e-14; % L

%% T helper cells
% Fp, Epitope dependent frequency of precursor CD4

% RhoNT: maximum proliferation rate for activated helper T cells
RhoNT=0.19; % day-1

%BetaNT: death rate of naive helper T cells
BetaNT=0.45; % day-1

%DeltaNT: maximum activation rate of naive helper T cells
DeltaNT=1.5; % day-1

% RhoAT: maximum proliferation rate for activated helper T cells
RhoAT=0.5973; % day-1

%BetaAT: death rate of activated helper T cells
BetaAT=0.18; % day-1

%% Initial conditions for state variables in the differential equations:
% Ag0: intitial therapeutic protein in the plasma compartment
Ag0=Dose; %#ok<NASGU> % pmole

% MS0: Maturation signal (MS), particularly, endotoxin, LPS
MS0=Endotoxin*Va; %#ok<NASGU> % ng

% MD0: the initial number of mature dendritic cells
MD0=0; %#ok<NASGU> % cells

% AgE0: initial amount of Ag in endosome
AgE0=0; %#ok<NASGU> % pmole

% pE0:	initial amount of T-epitope peptides from Ag digestion in endosome
pE0=ones(N,1)*0;  %#ok<NASGU> % pmole

% pME0:	initial amount of T-epitope-MHC-II complex in endosome
pME0=ones(6*N,1)*0;  %#ok<NASGU> % pmole

% pM0:	initial amount of T-epitope-MHC-II complex on dendritic cell membrane
pM0=ones(6*N,1)*0;  %#ok<NASGU> % pmole

% M0, free MHC-II molecule on dendritic cell membrane
M0=ones(6,1)*0; %#ok<NASGU> % pmole

% AT_N0: initial number of activated T helper cell derived from naive T cells
AT_N0=0.0; %#ok<NASGU> % cells

% AT_M0: initial number of activated T helper cell derived from memory T cells
Prolif0=ones(4,1)*0.0; %#ok<NASGU> % cells

%% Parameter vector
pars(1)=NA;
% Therapeutic protein PK parameters
pars(2)=Dose;
pars(3)=Va;
% T-epitope characteristics of therapeutic proteins
pars(4)=N;
pars(5:(4+6*N))=reshape(kon,6*N,1);
pars((5+6*N):(4+12*N))= reshape(koff, 6*N,1);
% Dendritic cells
pars(5+12*N)=BetaID;
pars(6+12*N)=DeltaID;
pars(7+12*N)=KMS;
pars(8+12*N)=BetaMD;
% Antigen presentation
pars(9+12*N)=AlphaAgE;
pars(10+12*N)=BetaAgE;
pars(11+12*N)=Beta_p;
pars(12+12*N)=BetaM;
pars(13+12*N)=Beta_pM;
pars(14+12*N)=kext;
pars(15+12*N)=kin;
pars(16+12*N)=KpM_N;
pars(17+12*N)=VD;
pars(18+12*N)=VE;
% T helper cells
pars(19+12*N)=RhoNT;
pars(20+12*N)=BetaNT;
pars(21+12*N)=DeltaNT;
pars(22+12*N)=RhoAT;
pars(23+12*N)=BetaAT;
pars(24+12*N)=Fp; %Epitope Dependent - Size 1
% Initial conditions as parameters in the differential equations:
pars(25+12*N)=ID0;
pars((26+12*N):(31+12*N))=ME0;
pars=pars'; %#ok<NASGU>

save Parameters.mat

function [kon,koff,N,ME0] = doNETMHCIIpan(SimType,ProteinLength,NA,numHLADR)
%% Collect -on, -off rates, number of epitopes and amount of initial MHC

rankcutoff  = 100;

[netmhciipanresult,N,numHLAalleles] = givekon();
Affinity_DR = reshape(netmhciipanresult(1:N*numHLAalleles,1),N,numHLAalleles);
Rank_DR = reshape(netmhciipanresult(1:N*numHLAalleles,2),N,numHLAalleles);

Affinity_DPQ=repmat([4000 4000 4000 4000],N,1); %place holder for other alleles (NOT USED in this version)

if numHLAalleles==1 %homozygot
    disp('Homozygote');
    rankadjustment = [Rank_DR<rankcutoff zeros(N,6-1)]; %Nx6 matrix
    ME0=[numHLADR; 0.0;  34E3/2; 34E3/2; 17.1E3/2; 17.1E3/2]/NA*1E12; % pmole
    Affinity_DR = [Affinity_DR zeros(N,1)];
else
    rankadjustment = [Rank_DR<rankcutoff zeros(N,6-2)]; %Nx6 matrix
    ME0=[numHLADR/2; numHLADR/2;  34E3/2; 34E3/2; 17.1E3/2; 17.1E3/2]/NA*1E12; % pmole
end

EpitopeSurvived = ones(N,6);
% EpitopeLength = 15;
% if ProteinLength>100
%     for i = 1:N
%         surviveP = ((1-(EpitopeLength-1)/(ProteinLength-1))^round(ProteinLength/EpitopeLength))>rand;
%         EpitopeSurvived(i,:) = surviveP*ones(1,6);
%     end
% end

% kon: on rate for for T-epitope-MHC-II binding
kon=SimType*EpitopeSurvived.*rankadjustment*8.64*1E-3; %  pM-1day-1

%	koff:	off rate for for T-epitope-MHC-II binding
koff=8.64*1E-3*[Affinity_DR Affinity_DPQ]*1E3; %  day-1

function [AffinityandRank_DR,N,numHLAalleles] = givekon()
% Read out.dat file
if exist('shortout.dat')==0 %#ok<EXIST>
    getshorty();
end
table = readtable('shortout.dat');
AffinityandRank_DR = table{:,9:10}; 
% Structured as 
% [Epitope1 Allele1; 
%     Epitope2 Allele1; 
%     Epitope1 Allele2;
% Epitope2 Allele2]

N = length(unique(table{:,3})); %Number of epitopes

numHLAalleles = length(unique(table{:,2})); %number of HLAs

%% Post process NetMHCIIPan v3.2 out.dat file
function getshorty()

ind = 0;
fid = fopen('out.dat');
fidshort = fopen('shortout.dat','w');

if contains(pwd,'seqs\10') %exanitide patch
    lists = [57 270]; %use +3 for version 3.1
elseif contains(pwd,'seqs\18') %exanitide patch 
    lists = [32 65]; %use +3 for version 3.1
elseif contains(pwd,'seqs\19') %exanitide patch 
    lists = [14:10:254 263:10:503]; %use +3 for version 3.1
else
    lists = [12 21]; %use +3 for version 3.1
end

tline = fgetl(fid);
while ischar(tline)
    ind = ind + 1;
    if ismember(ind,lists)
        %     disp(tline)
        fprintf(fidshort,'%s\n',tline);
    end
    tline = fgetl(fid);
end

fclose(fid);
fclose(fidshort);
