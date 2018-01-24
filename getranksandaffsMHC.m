function [Affinity_DR,Rank_DR] = getranksandaffsMHC()

table = readtable('out.dat');
Affinity_DR(1,1) = table{1,'nM'};
Affinity_DR(1,2) = table{1,'nM_1'};
Rank_DR(1,1) = table{1,'Rank'};
Rank_DR(1,2) = table{1,'Rank_1'};