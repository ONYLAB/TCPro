function cost = LamberthFig4B(newpars)

% pars(32+12*N)=newpars(1); %RhoNT;
% pars(35+12*N)=newpars(2); %RhoAT;
% pars(36+12*N)=newpars(3); %BetaAT;
dlmwrite('newparamsforoptimization.dat',newpars);

seqtable = readtable('runpeptides.dat');

for i = 1:height(seqtable)
    
    seqid = seqtable{i,'SEQNAME'}{1};
    epitope = seqtable{i,'SEQUENCE'}{1};
    disp(seqid);
    %     isHuman(i) = isAntigenHuman(epitope);
    if mod(i,2)==0 %isHuman(i)==0 %Exogenous
        mkdir(seqid)
        cd(seqid)
        [responsesummary,significantresponsesummary,kon,responsevector,significancevector] = runalldonor(epitope);
        kons{i}=kon;
        allmeanresponses{i} = responsevector;
        allsignificance{i} = significancevector;
        CUMresponses(i,:)=responsesummary;
        CUMsignificance(i,:)=significantresponsesummary;
        cd ..
        %         save
    end
end

for i = 1:height(seqtable)
    if mod(i,2)==0 %isHuman(i)==0 %Exogenous
        for j = 1:51
            
            responses = allmeanresponses{1,i}{1,j}>=2;
            significance = allsignificance{1,i}{1,j}<=0.05;
            
            sigresponse(i,j,:) = responses.*significance;
            CUMsigresponses(i,j) = sign(sum(sigresponse(i,j,1:4))).*sigresponse(i,j,5);
            
        end
    end
end

% save
costfunction = [0 10 0 22 0 18 0 24];%[12 10 10 22 10 18 10 24]; %responders Lamberth et. al.
cost = norm(100*sum(CUMsigresponses')/51 - costfunction);

dlmwrite('results.dat',[newpars(:)' cost],'-append')