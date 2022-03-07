% ----------------------------------------------------------------------- %
%   Version: 1.0                                                          %
%   Author:  Kenneth Brenner dos Anjos Benicio                            %
%   Date:    19/02/2022                                                   %
%   E-mail:  kenneth@gtel.ufc.br                                          %
% ----------------------------------------------------------------------- %
%% Obtendo o comportamento médio dos índices
clc;
close all;

% Define o intervalo de valores que será utilizado para obtenção do Kopt.
K = 2:20;
runs = 100;

% index1 = zeros(runs,max(K)-1,3);
% index2 = zeros(runs,max(K)-1,3);

% Processo de obtenção do comportamento médio dos índices para um
% determinado número de realizações considerando o valor de moda.

% for r = 1:runs
%     r
%     parfor k = K
%         [~, ~, ~, aux1] = clusterizacao.kmseq(Data, k,'n');
%         [~, ~, ~, aux2] = clusterizacao.kmbat(Data, k,'n');
%         index1(r,k-1,:) = aux1;
%         index2(r,k-1,:) = aux2;
%     end
% end
% 
% [~,A1] = min(index1(:,:,1), [],2);
% [~,A2] = min(index1(:,:,2), [],2);
% [~,A3] = max(index1(:,:,3), [],2);
% A1 = mode(A1);
% A2 = mode(A2);
% A3 = mode(A3);
% 
% [~,B1] = min(index2(:,:,1), [],2);
% [~,B2] = min(index2(:,:,2), [],2);
% [~,B3] = max(index2(:,:,3), [],2);
% B1 = mode(B1);
% B2 = mode(B2);
% B3 = mode(B3);

% Processo de obtenção do comportamento médio dos índices para um
% determinado número de realizações considerando um processo de
% obtenção de comportamento médio.

index1 = zeros(max(K)-1,3);
index2 = zeros(max(K)-1,3);

for r = 1:runs
    r
    parfor k = K
        [~, ~, ~, aux1] = clusterizacao.kmseq(Data, k,'n');
        [~, ~, ~, aux2] = clusterizacao.kmbat(Data, k,'n');
        index1(k-1,:) = index1(k-1,:) + aux1;
        index2(k-1,:) = index2(k-1,:) + aux2;
    end
end

% Desempenho médio após os experimentos independentes
index1 = index1./runs;
index2 = index2./runs;

% Tabelas para os índices
TabelaIndex1 = array2table(index1);
TabelaIndex1.Properties.RowNames = {'2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30'};
TabelaIndex1.Properties.VariableNames = {'Dunn' 'Davies-Bouldin' 'Calinski-Harabasz'};
table2latex(TabelaIndex1,'KMSEQ_index.tex');

TabelaIndex2 = array2table(index2);
%TabelaIndex2.Properties.RowNames = {'2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15'};
TabelaIndex2.Properties.RowNames = {'2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30'};
TabelaIndex2.Properties.VariableNames = {'Dunn' 'Davies-Bouldin' 'Calinski-Harabasz'};
table2latex(TabelaIndex2,'KMBAT_index.tex');

%% Plotando as curvas de desempenho.

figure;
subplot(3,1,1);
semilogy(K,index1(:,1),'b','linewidth',3);
title('Indice de Dunn x Clusters para K-Means Sequencial');
xlabel('Número de Clusters');
ylabel('Indice de Dunn');
grid on;
xlim tight
set(gca, "fontsize", 8);

subplot(3,1,2);
semilogy(K,index1(:,2),'y','linewidth',3);
title('Indice de Davies-Bouldin x Clusters para K-Means Sequencial');
xlabel('Número de Clusters');
ylabel('Indice DB');
grid on;
xlim tight
set(gca, "fontsize", 8);

subplot(3,1,3);
semilogy(K,index1(:,3),'m','linewidth',3);
title('Indice de Calinski-Harabasz x Clusters para K-Means Sequencial');
xlabel('Número de Clusters');
ylabel('Indice CH');
grid on;
xlim tight
set(gca, "fontsize", 8);



figure;
subplot(3,1,1);
semilogy(K,index2(:,1),'b','linewidth',3);
title('Indice de Dunn x Clusters para K-Means Agrupado');
xlabel('Número de Clusters');
ylabel('Indice Dunn');
grid on;
xlim tight
set(gca, "fontsize", 8);

subplot(3,1,2);
semilogy(K,index2(:,2),'y','linewidth',3);
title('Indice de Davies-Bouldin x Clusters para K-Means Agrupado');
xlabel('Número de Clusters');
ylabel('Indice DB');
grid on;
xlim tight
set(gca, "fontsize", 8);

subplot(3,1,3);
semilogy(K,index2(:,3),'m','linewidth',3);
title('Indice de Calinski-Harabasz x Clusters para K-Means Agrupado');
xlabel('Número de Clusters');
ylabel('Indice CH');
grid on;
xlim tight
set(gca, "fontsize", 8);
