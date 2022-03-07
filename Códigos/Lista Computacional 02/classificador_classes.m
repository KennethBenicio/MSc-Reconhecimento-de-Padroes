% ----------------------------------------------------------------------- %
%   Version: 1.0                                                          %
%   Author:  Kenneth Brenner dos Anjos Benicio                            %
%   Date:    19/02/2022                                                   %
%   E-mail:  kenneth@gtel.ufc.br                                          %
% ----------------------------------------------------------------------- %
% Esse código gera uma simples tabela descritiva das classes do conjunto 
% de dados que será utilizado no processo de classificação.
aux = Dados(:,end);
classe1 = (aux == 1);
classe2 = (aux == 2);
classe1 = Dados(classe1,:);
classe2 = Dados(classe2,:);

dados_c1 = [length(classe1) rank(cov(classe1)) rcond(cov(classe1))];
dados_c2 = [length(classe2) rank(cov(classe2)) rcond(cov(classe2))];

tabela = [dados_c1; dados_c2;];
Tabela = array2table(tabela);
Tabela.Properties.RowNames = {'Classe 1' 'Classe 2'};
Tabela.Properties.VariableNames = {'Amostras' 'Posto Mat.' 'N. de Cond.'};
table2latex(Tabela,'classificador_classes.tex');