% ----------------------------------------------------------------------- %
%   Version: 1.0                                                          %
%   Author:  Kenneth Brenner dos Anjos Benicio                            %
%   Date:    19/02/2022                                                   %
%   E-mail:  kenneth@gtel.ufc.br                                          %
% ----------------------------------------------------------------------- %
%% Primeira etapa desconsiderando a utilização da compressão de dados com a PCA

Nr     = 100;
Ptrain = 80;

% Item 1, Item 2, Item 3, Item 4 e Item 5 sem PCA.
PCA    = 'n';
[STATS1, cfsmtx1, ~, ~, ~, ~, ~, ~] = classificadores.CGQ12(Dados,Nr,Ptrain,PCA);
[STATS2, cfsmtx2, ~, ~, ~, ~, ~, ~] = classificadores.CGQ44(Dados,Nr,Ptrain,PCA);
[STATS3, cfsmtx3, ~, ~, ~, ~, ~, ~] = classificadores.CGQ17(Dados,Nr,Ptrain,PCA);
[STATS4, cfsmtx4, ~, ~, ~, ~, ~, ~] = classificadores.CGQ39(Dados,Nr,Ptrain,PCA);
[STATS5, cfsmtx5, ~, ~, ~]          = classificadores.LMQ(Dados,Nr,Ptrain,PCA);


%% Segunda etapa considerando a utilização da compressão de dados com a PCA

% Item 1, Item 2, Item 3, Item 4 e Item 5 com PCA.
PCA    = 'y';
[STATS6, cfsmtx6, ~, ~, ~, ~, ~, ~] = classificadores.CGQ12(Dados,Nr,Ptrain,PCA);
[STATS7, cfsmtx7, ~, ~, ~, ~, ~, ~] = classificadores.CGQ44(Dados,Nr,Ptrain,PCA);
[STATS8, cfsmtx8, ~, ~, ~, ~, ~, ~] = classificadores.CGQ17(Dados,Nr,Ptrain,PCA);
[STATS9, cfsmtx9, ~, ~, ~, ~, ~, ~] = classificadores.CGQ39(Dados,Nr,Ptrain,PCA);
[STATS10, cfsmtx10, ~, ~, ~]        = classificadores.LMQ(Dados,Nr,Ptrain,PCA);

%% Montagem da tabela desejada

tabela = [STATS1; STATS2; STATS3; STATS4; STATS5; STATS6; STATS7; STATS8; STATS9; STATS10;];
Tabela = array2table(tabela);
Tabela.Properties.RowNames = {'Item 1' 'Item 2' 'Item 3' 'Item 4' 'Item 5' 'Item 1 + PCA' 'Item 2 + PCA' 'Item 3 + PCA' 'Item 4 + PCA' 'Item 5 + PCA'};
Tabela.Properties.VariableNames = {'Média' 'Desvio' 'Mediana' 'Mínimo' 'Máximo'};
table2latex(Tabela,'classificador_tabela.tex');

cfsmtx = {cfsmtx1 cfsmtx2 cfsmtx3 cfsmtx4 cfsmtx5 cfsmtx6 cfsmtx7 cfsmtx8 cfsmtx9 cfsmtx10};
for i = 1:10
    aux = cfsmtx{i}; 
    for j = 1:2
       if j == 1
           close all;
           str1 = 'Item_';
           str2 = num2str(i);
           str3 = '_min.pdf';
           nome = strcat(str1,str2,str3);
           confusionchart(aux{j})
           print(gcf, nome, '-dpdf', '-bestfit', '-loose')
       else
           close all;
           str1 = 'Item_';
           str2 = num2str(i);
           str3 = '_max.pdf';
           nome = strcat(str1,str2,str3);
           confusionchart(aux{j})
           print(gcf, nome, '-dpdf', '-bestfit', '-loose')
       end
    end
end
