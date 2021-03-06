function dydt=f(t,y,pars) %#ok<INUSL>

% Parameter vector
NA=pars(1);
% Therapeutic protein PK parameters
Dose=pars(2); %#ok<NASGU>
Vp=pars(3);
% T-epitope characteristics of therapeutic proteins
N=pars(4);
kon=reshape(pars(5:(4+6*N)),N,6);
koff= reshape(pars((5+6*N):(4+12*N)), N,6);
% Dendritic cells
BetaMS=pars(5+12*N);
BetaID=pars(6+12*N);
DeltaID=pars(7+12*N);
KMS=pars(8+12*N);
BetaMD=pars(9+12*N);
% Antigen presentation
konN=pars((10+12*N):(15+12*N));
koffN=pars((16+12*N):(21+12*N));
AlphaAgE=pars(22+12*N);
BetaAgE=pars(23+12*N);
Beta_p=pars(24+12*N);
BetaM=pars(25+12*N);
Beta_pM=pars(26+12*N);
kext=pars(27+12*N);
kin=pars(28+12*N);
KpM_N=pars(29+12*N);
VD=pars(30+12*N);
VE=pars(31+12*N);
% T helper cells
RhoNT=pars(32+12*N);
BetaNT=pars(33+12*N);
DeltaNT=pars(34+12*N);
RhoAT=pars(35+12*N);
BetaAT=pars(36+12*N);
Fp=pars(37+12*N:36+13*N);
% Initial conditions as parameters in the differential equations:
ID0=pars(37+13*N); %#ok<NASGU>
cp0=pars(38+13*N);
ME0=pars((39+13*N):(44+13*N));
NT0=pars((45+13*N):(44+14*N));
MT0=pars((45+14*N):(44+15*N));

% State variables for the differential equations
% Ag, antigenic protein in the plasma compartment
Ag=y(1);
% MS, maturation signal for Immature dendritic cells
MS=y(2);
% ID,	Immature dendritic cells
ID=y(3);
% MD,	mature dendritic cells
MD=y(4);
% cpE, endogenous competing protein in the endosomes
cpE=y(5);
% cptE, endogenous competing peptide in the endosomes
cptE=y(6);
% cptME, competing peptide-MHC complex in the endosomes
cptME=y(7:12);
% cptM, competing peptide-MHC complex on dendritic cell membrane
cptM=y(13:18);
%AgE, Ag in the endosomes
AgE=y(19);
% pE,	free epitope peptides from Ag digestion
pE=y((19+1):(19+N));
% ME(1:6),	free MHC II molecule in endosome
ME=y((20+N):(25+N));
% pME(N,6),	epitope peptide-MHC II complex in endosomes
% pME is in a matrix form, N epitopes against 6 possible MHC alleles
pME=reshape(y((26+N):(25+7*N)),N,6);
% pM(N,6),	epitope peptide-MHC II complex on dendritic cell membrane
pM=reshape(y((26+7*N):(25+13*N)),N,6);
% M, free MHC II molecule on dendritic cell membrane
M=y((26+13*N):(31+13*N));
% NT, na�ve helper T cells
NT=y((32+13*N):(32+13*N));
% AT_N	activated helper T cells
AT_N=y((33+13*N):(32+14*N));
% AT_M	activated helper T cells
AT_M=y((33+14*N):(32+15*N));
% MT, memory helper T cells
MT=y((33+15*N):(32+16*N));
% FT, functional helper T cells
FT=y((33+16*N):(32+17*N)); %#ok<NASGU>

% Calculate functions for helper T cells activation, proliferation, or differentiation

% Adding up all the pM molecules from one epitope against 6 different MHC alleles
pM_NUMBER=(pM+0.0*repmat(cptM',N,1))*ones(6,1)*1E-12*NA;%*(1-sign(kon(1,1)));

% The activation function D for helper T cells
D_N=(MD/(MD+sum(Fp.*NT)+sum(AT_N)+0.0*sum(AT_M)+0.0*sum(MT)))*(((pM_NUMBER)./(pM_NUMBER+KpM_N)));   % for naive T cells
if isnan(D_N)
    D_N=0.0;
end
% The proliferation/differentiation function E for helper T cells
E_N=(MD/(MD+sum(Fp.*NT)+sum(AT_N)+0.0*sum(AT_M)+0.0*sum(MT)))*((pM_NUMBER-KpM_N)./(pM_NUMBER+KpM_N));  % for naive T cells
if isnan(E_N)
    E_N=0.0;
end
% Differential equations
% Ag, y(2), total amount of antigenic protein in the well, pmole
dydt(1,1)=-(ID+MD)*AlphaAgE*VD*(Ag/Vp);

% MS, y(5), maturation signal, particularly, LPS, for immature dendritic cells, ng
dydt(2,1)=-(ID+MD)*AlphaAgE*VD*(MS/Vp);

% ID, y(6), immature dendritic cells, cells
dydt(3,1)=-BetaID*ID-DeltaID*ID*(MS/Vp)/((MS/Vp)+KMS);

% MD,	y(7), mature dendritic cells , cells
dydt(4,1)=DeltaID*ID*(MS/Vp)/((MS/Vp)+KMS)-BetaMD*MD;

% cpE, y(8), endogenous competing protein in endosome, pmole
dydt(5,1)=0.0;%BetaAgE*(cp0-cpE);

% cptE, y(9), endogenous competing peptide in endosome, pmole
dydt(6,1)=0.0;%BetaAgE*cpE -Beta_p*cptE -(konN.*(cptE*(ME/VE)))'*ones(6,1)+(koffN.*cptME)'*ones(6,1);

% cptME, y(10:15), endogenous competing peptide-MHC complex in endosome, pmole
dydt((7:12),1)=0.0;%konN*cptE.*(ME/VE)-koffN.*cptME -Beta_pM*cptME -kext*cptME;

% cptM, y(16:21), endogenous peptide-MHC complex on dendritic cell membrane, pmole
dydt((13:18),1)=0.0;%kext*cptME -koffN.*cptM;

%AgE, y(22), Ag in the endosome, pmole
dydt(19,1)=AlphaAgE*VD*(Ag/Vp)-BetaAgE*AgE;

% pE, y(23:(22+N)), T-epitope peptides from Ag digestion, pmole
dydt((20:(19+N)),1)=BetaAgE*AgE-Beta_p*pE-kon.*(pE*(ME'/VE))*ones(6,1)+koff.*pME*ones(6,1);

% ME, y((23+N):(28+N)),	free MHC-II molecule in endosome, pmole
dydt((20+N):(25+N),1)= BetaM*(ME0-ME)-(ones(1,N)*(kon.*(pE*(ME'/VE))))'+(ones(1,N)*(koff.*pME))'-konN*cptE.*(ME/VE)+koffN.*cptME +kin*M;

% pME, y((29+N):(28+7*N)), T-epitope-MHC-II complex in endosome, pmole
dydt(((26+N):(25+7*N)),1)= reshape(kon.*(pE*(ME'/VE))-koff.*pME -Beta_pM*pME -kext*pME,6*N,1);

% pM,	y((29+7*N):(28+13*N)), T-epitope-MHC-II complex on dendritic cell membrane, pmole
dydt((26+7*N):(25+13*N),1)=reshape(kext*pME -koff.*pM,6*N,1);

% M, y((29+13*N):(34+13*N)), free MHC-II molecule on dendritic cell menbrane, pmole
dydt((26+13*N):(31+13*N),1)=-kin*M +(ones(1,N)*(koff.*pM))'+koffN.*cptM;

% NT, y((35+13*N):(34+14*N)), na�ve helper T cells,  cells
dydt((32+13*N):(32+13*N),1)=MD*RhoNT-BetaNT*NT-sum(DeltaNT*D_N.*Fp.*NT);                                    %PREVIOUSLY: BetaNT*(NT0-NT)-DeltaNT*D_N.*NT;

%AT_N, y((35+14*N):(34+15*N)), activated helper T cells derived from NT, cells
dydt((33+13*N):(32+14*N),1)=DeltaNT*D_N.*Fp.*NT+RhoAT.*E_N.*AT_N-BetaAT*AT_N;

% AT_M, y((35+15*N):(34+16*N)),	activated helper T cells derived from MT, cells
sAt = E_N>=0;
dydt((33+14*N):(32+15*N),1)= MD*RhoNT + sAt.*RhoAT.*E_N.*AT_N;

% NT0Rest, y((35+16*N):(34+17*N)), NT0Rest
dydt((33+15*N):(32+16*N),1)=0.0;%MD*0.1616-0.3048*MT;%This says MT but it's actually NT0Rest

%place holder, y((35+17*N):(34+18*N)) 
dydt((33+16*N):(32+17*N),1)= repmat(0.0,N,1); %#ok<REPMAT>