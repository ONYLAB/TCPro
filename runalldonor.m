function [Response,IncoBinary,ELISBinary,COMBOBinary,Inco,ELIS,COMBO,corrs] = runalldonor(traj)

% IncorporationResponse(n,SimType+1,DayLimit-4)
% ELISpotResponse(n,SimType+1)

load pdMS0.mat;

SIcutoff = 1.9;
sigcutoff = 0.05;

ProteinLength = 40*ones(1,19);
ProteinLength(10) = 1350; 

AbzI = [2 4 6 7 9 10];
A33 = 10;
AbzIList = 1:50;
AbzIIList = 51:100;%[51:54 56:100];
AbzII = [12 13 14 15 16 17 18 19];
Exenatide = 18;
KLH = 19;

% Precursor frequencies
     
Fp = [0.00    0.31    0.00    0.52    0.00    0.39    0.45    0.00    0.91    0.61    0.00    0.91    0.73    0.56    0.78    0.98    0.58    1.78  9.86]/1e6;
FpSTD = [0    0.08    0       0.11    0       0.09    0.10    0       0.11    0.13    0       0.23    0.15    0.12    0.15    0.18    0.12    0.29 1.45]/1e6;

tic
for s = AbzI
    cd(num2str(s))
    for i = AbzIList
        
        cd(num2str(i));
        if s==A33
            SampleConcentration = 0.3e-6;
        else
            SampleConcentration = 5e-6;
        end
        
        rng(i+traj*50); %So that Fps are randomly consistent
        [Response{s,i},pval{s,i}] = Main_human(i+traj*50,ProteinLength(s),SampleConcentration,random('Normal',Fp(s),FpSTD(s)),pdMS0); %#ok<AGROW>
        temp = Response{s,i};
        temp2 = pval{s,i};
        IncoBinary(s,i) = sign(sum(((temp(1:4))>SIcutoff).*(temp2(1:4)<sigcutoff))); %#ok<AGROW>
        ELISBinary(s,i) = (temp(5)>SIcutoff).*(temp2(5)<sigcutoff); %#ok<AGROW>
        cd ..
        
    end
    disp(s)
    cd ..
end
toc
% % % %
for s = AbzII
    cd(num2str(s))
    for i = AbzIIList
        if i~=55
            cd(num2str(i));
            if s==Exenatide
                SampleConcentration = 0.3e-6;
            elseif s== KLH
                SampleConcentration = 20*0.3e-6;                
            else
                SampleConcentration = 5e-6;
            end
            
            rng(i+traj*50); %So that Fps are randomly consistent           
            [Response{s,i},pval{s,i}] = Main_human(i+traj*50,ProteinLength(s),SampleConcentration,random('Normal',Fp(s),FpSTD(s)),pdMS0); %#ok<AGROW>
            temp = Response{s,i};
            temp2 = pval{s,i};
            IncoBinary(s,i) = sign(sum(((temp(1:4))>SIcutoff).*(temp2(1:4)<sigcutoff))); %#ok<AGROW>
            ELISBinary(s,i) = (temp(5)>SIcutoff).*(temp2(5)<sigcutoff); %#ok<AGROW>
            cd ..
        end
    end
    disp(s)
    cd ..
end

COMBOBinary = IncoBinary.*ELISBinary;

Inco = sum(IncoBinary(:,1:50)')'*2 + sum(IncoBinary(:,[51:54 56:100])')'/49*100;
ELIS = sum(ELISBinary(:,1:50)')'*2 + sum(ELISBinary(:,[51:54 56:100])')'/49*100;
COMBO = sum(COMBOBinary(:,1:50)')'*2 + sum(COMBOBinary(:,[51:54 56:100])')'/49*100;

% Correlation
Inco1 = IncoBinary(1:10,1:50); 
Inco2 = IncoBinary(11:end,[51:54 56:100]);
ELIS1 = ELISBinary(1:10,1:50); 
ELIS2 = ELISBinary(11:end,[51:54 56:100]);
corrs = [diag(corr(Inco1',ELIS1'));diag(corr(Inco2',ELIS2'))];

save([num2str(traj) '_matlabRcutoff100Prol0597.mat'])
