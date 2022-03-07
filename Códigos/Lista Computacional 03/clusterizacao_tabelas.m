% ----------------------------------------------------------------------- %
%   Version: 1.0                                                          %
%   Author:  Kenneth Brenner dos Anjos Benicio                            %
%   Date:    19/02/2022                                                   %
%   E-mail:  kenneth@gtel.ufc.br                                          %
% ----------------------------------------------------------------------- %
%% Geração das tabelas de dados que serão utilizadas no relatório.
% Os links abaixo foram utilizados por mim para desenvolver o código
% a seguir que cria as tabelas para as posições dos protótipos e para as 
% estatísticas de agrupamento. Em seguinda, realizo a conversão das tabelas
% para arquivos .tex para facilitar a visualização no relatório,

% https://www.mathworks.com/matlabcentral/answers/393282-how-can-i-create-a-descriptive-statistics-table-for-columns
% https://www.mathworks.com/matlabcentral/answers/395567-how-can-i-transpose-a-dataset-or-table
% https://www.mathworks.com/matlabcentral/fileexchange/69063-matlab-table-to-latex-conversor

clc;
close all;

% Valor ótimo encontrado segundo a análise dos gráficos gerados pela função
% clusterizacao_kopt.m
Kopt = 10;

[~, Particao, W, ~] = clusterizacao.kmseq(Data, Kopt, 'n');

% Tabela da posição dos protótipos para K-means sequencial.
TabelasKMSEQ = cell(Kopt,1); 
TabelaWKMSEQ = array2table(W);
TabelaWKMSEQ.Properties.RowNames = {'w_{1}' 'w_{2}' 'w_{3}' 'w_{4}' 'w_{5}' 'w_{6}' 'w_{7}' 'w_{8}' 'w_{9}' 'w_{10}'};
TabelaWKMSEQ.Properties.VariableNames = {'Par1' 'Par2' 'Par3' 'Par4' 'Par5' 'Par6'};
table2latex(TabelaWKMSEQ,'KMSEQ_W.tex');

% Tabela das estatísticas de agrupamento para K-means sequencial.
for k = 1:Kopt
    aux = array2table(Particao{k});
    aux = varfun(@(x) [min(x); max(x); median(x); mean(x); std(x)], aux);
    aux = table2cell(aux);
    tabela = cell2table(aux.');
    tabela.Properties.RowNames = {'Par1' 'Par2' 'Par3' 'Par4' 'Par5' 'Par6'};
    tabela.Properties.VariableNames = {'Mínimo' 'Máximo' 'Mediana' 'Média' 'Desvio'};
    TabelasKMSEQ{k,1} = tabela;
end



[~, Particao, W, ~] = clusterizacao.kmbat(Data, Kopt, 'n');

% Tabela da posição dos protótipos para K-means em batch.
TabelasKMBAT = cell(Kopt,1);
TabelaWKMBAT = array2table(W);
TabelaWKMBAT.Properties.RowNames = {'w_{1}' 'w_{2}' 'w_{3}' 'w_{4}' 'w_{5}' 'w_{6}' 'w_{7}' 'w_{8}' 'w_{9}' 'w_{10}'};
TabelaWKMBAT.Properties.VariableNames = {'Par1' 'Par2' 'Par3' 'Par4' 'Par5' 'Par6'};
table2latex(TabelaWKMBAT,'KMBAT_W.tex');

% Tabela das estatísticas de agrupamento para K-means em batch.
for k = 1:Kopt
    clear aux
    aux = array2table(Particao{k});
    aux = varfun(@(x) [min(x); max(x); median(x); mean(x); std(x)], aux);
    aux = table2cell(aux);
    tabela = cell2table(aux.');
    tabela.Properties.RowNames = {'Par1' 'Par2' 'Par3' 'Par4' 'Par5' 'Par6'};
    tabela.Properties.VariableNames = {'Mínimo' 'Máximo' 'Mediana' 'Média' 'Desvio'};
    TabelasKMBAT{k,1} = tabela;
end



% Aqui apenas geros arquivos .tex para serem inseridos no relátorio. Desse
% modo, essa é apenas uma etapa para melhorar a visualização dos dados.
for k = 1:Kopt
    str1 = 'KMSEQ_Agrupamento_';
    str2 = 'KMBAT_Agrupamento_';
    str3 = '.tex';
    aux  = num2str(k);
    filenameKMSEQ = strcat(str1,aux,str3);
    filenameKMBAT = strcat(str2,aux,str3);
    table2latex(TabelasKMSEQ{k},filenameKMSEQ);
    table2latex(TabelasKMBAT{k},filenameKMBAT);
end