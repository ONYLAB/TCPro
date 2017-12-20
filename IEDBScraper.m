function kon = IEDBScraper(peptidesequence,nallele,allelelist)

url = 'http://tools-cluster-interface.iedb.org/tools_api/mhcii/';
method = 'method';
methodSelect = 'recommended';
sequence_text = 'sequence_text';
sequence_textSelect = peptidesequence;
allele = 'allele';
alleleSelect = '';
for i = 1:nallele
    if i>1
        alleleSelect = [alleleSelect ','];
    end
    al = allelelist{i};
    alleleSelect = [alleleSelect 'HLA-' al(1:4) '*' al(6:7) ':' al(8:9)];
    % alleleSelect = 'HLA-DRB1*10:01,HLA-DRB1*12:01';
end


response = webwrite(url,method,methodSelect,sequence_text,sequence_textSelect,allele,alleleSelect);

fileID = fopen('IEDBresult.dat','w');
nbytes = fprintf(fileID,'%c',response);
fclose(fileID);

res = readtable('IEDBresult.dat');
for i = 1:nallele
    temp = res(strcmp(res.allele,allelelist{i}),'percentile_rank');
    rank(i,1) = temp.mean_percentile_rank;
end

% convert rank to affinity
% kon = f(rank);