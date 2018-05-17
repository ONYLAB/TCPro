function [Response,IncoBinary,ELISBinary,Inco] = runalldonor(traj)

% IncorporationResponse(n,SimType+1,DayLimit-4)
% ELISpotResponse(n,SimType+1)

SIcutoff = 1.9;
sigcutoff = 0.05;

ProteinLength = 40*ones(1,18);
ProteinLength(10) = 1350; 

AbzI = [2 4 6 7 9 10];
A33 = 10;
AbzIList = 1:50;
AbzIIList = 51:100;%[51:54 56:100];
AbzII = [12 13 14 15 16 17 18];
Exenatide = 18;

% Precursor frequencies
Fp = [0.00    0.04    0.00    0.26    0.00    0.18    0.18    0.00    0.51    0.26    0.00    0.51    0.18    0.04    0.22    0.31    0.04    1.23]/1e6;

% Medium Baseline Precursor frequencies
% Fp = [0.00    0.18    0.00    0.48    0.00    0.35    0.35    0.00    0.70    0.48    0.00    0.70    0.18    0.09    0.22    0.44    0.13    1.10]/1e6;

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
        
        [Response{s,i},pval{s,i}] = Main_human(i+traj*50,ProteinLength(i),SampleConcentration,Fp(s)); %#ok<AGROW>
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
            else
                SampleConcentration = 5e-6;
            end
            
            [Response{s,i},pval{s,i}] = Main_human(i+traj*50,ProteinLength(i),SampleConcentration,Fp(s)); %#ok<AGROW>
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

Inco = sum(IncoBinary')'*2;
ELIS = sum(ELISBinary')'*2;
% Inco = Inco([AbzI AbzII]);
% ELIS = ELIS([AbzI AbzII]);

save([num2str(traj) '_matlabRcutoff100Prol0597.mat'])
