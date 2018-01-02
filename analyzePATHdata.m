function analyzePATHdata()

% read .xls
thepwd = pwd;
JoeTable = readtable([thepwd '/GeneVariants/cleaned_path_data.csv']);
numberofpatients = height(JoeTable);

%list of all alleles
uniquelistofalleles = unique([JoeTable{:,9} JoeTable{:,10}]); %#ok<NASGU>
for i = 1:length(uniquelistofalleles)
    if length(uniquelistofalleles{i,1})==4
        uniquelistofalleles{i,1} = ['HLA-DRB1*0' uniquelistofalleles{i,1}];
    else
        uniquelistofalleles{i,1} = ['HLA-DRB1*' uniquelistofalleles{i,1}];
    end
end
save('uniquelistofalleles.mat','uniquelistofalleles');

% FLseq = fastaread([thepwd '/Drugs/FullLength/Advate.fasta']);
% disp('Start FL');
% [~,~,~,cleaverankscoreFL,affinityrankFL] = obtain_cleave_aff(FLseq.Sequence,'FullLength.txt',uniquelistofalleles);
load('FullLength.mat','cleaverankscore','affinityrank');
cleaverankscoreFL = cleaverankscore;
affinityrankFL = affinityrank;
clear('cleaverankscore','affinityrank');

BDDseq = fastaread([thepwd '/Drugs/BDD/Refacto.fasta']);
disp('Start BDD');
[~,~,~,cleaverankscoreBDD,affinityrankBDD] = obtain_cleave_aff(BDDseq.Sequence,'BDD.txt',uniquelistofalleles);

save('DrugsRanks.mat','cleaverankscoreFL','affinityrankFL','cleaverankscoreBDD','affinityrankBDD');

for i = 1:numberofpatients
    disp(['Patient ' num2str(i)]);
    if strcmp(JoeTable{i,'AnyDrug'},'TRUE')
        LocationID = JoeTable{i,'X'}{1}(7:9);
        PatID = num2str(JoeTable{i,'Pat_ID'});
        patiendindex = [LocationID '-' PatID '.txt'];
        patseqfasta = fastaread([thepwd '/GeneVariants/patient_seqs/' patiendindex]);
        patseqfasta = patseqfasta.Sequence;
        seqlength(i) = length(patseqfasta);
        listofalleles{1,1} = JoeTable{i,9}{1};
        listofalleles{2,1} = JoeTable{i,10}{1};
        if strcmp(JoeTable{i,9}{1},JoeTable{i,10}{1})
            listofalleles(2) = [];
            disp('homozygot patient');
        end
        
        for ax = 1:length(listofalleles)
            if length(listofalleles{ax,1})==4
                listofalleles{ax,1} = ['HLA-DRB1*0' listofalleles{ax,1}];
            else
                listofalleles{ax,1} = ['HLA-DRB1*' listofalleles{ax,1}];
            end
        end
        
        [minvectors,~,~,~,~] = obtain_cleave_aff(patseqfasta,patiendindex,listofalleles);
        if JoeTable{i,'BDDDrug'}==1
            cleanBDDaff = cleaner(uniquelistofalleles,listofalleles,affinityrankBDD);
            BDDresult(i,:) = howmany(cleanBDDaff,cleaverankscoreBDD,minvectors);
        end
        if JoeTable{i,'FullLengthDrug'}==1
            cleanFLaff = cleaner(uniquelistofalleles,listofalleles,affinityrankFL);
            FLresult(i,:) = howmany(cleanFLaff,cleaverankscoreFL,minvectors);
        end
    end
    save
end

function [minvectors,pepcontainer,cleavescore,cleaverankscore,affinityrank] = obtain_cleave_aff(fastaseq,patid,listofalleles)

[pepcontainer,cleavescore,cleaverankscore,affinityrank] = extractcleavability(fastaseq,listofalleles);
save(strrep(patid,'txt','mat'),'cleavescore','cleaverankscore','affinityrank','pepcontainer');

[~,ind] = min(min(affinityrank')'.*cleaverankscore);

minvectors = [min(cleaverankscore) min(min(affinityrank)) cleaverankscore(ind) affinityrank(ind)];

function cleanaff = cleaner(uniquelistofalleles,listofalleles,affinityrank)

%find the indices
for i = 1:length(listofalleles)
    ind(i) = find(strcmp(listofalleles{i},uniquelistofalleles)==1);
end

if length(ind)==1
    ind = [ind ind];
end

cleanaff = affinityrank(:,ind);

function result = howmany(cleanaff,cleaverankscore,minvectors)

cleaverankonly = minvectors(1);
affrankonly = minvectors(2);
bothrankandaff = minvectors(3:4);

numbettercleavers = sum(cleaverankscore<cleaverankonly);
numbetteraffiners = sum(sum(cleanaff<affrankonly));

numberboth = sum((cleaverankscore<bothrankandaff(1)).*sign(sum((cleanaff<bothrankandaff(2))'))');

result = [numbettercleavers numbetteraffiners numberboth];