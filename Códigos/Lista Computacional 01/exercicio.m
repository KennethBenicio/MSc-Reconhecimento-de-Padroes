% ----------------------------------------------------------------------- %
%   Version: 1.0                                                          %
%   Author:  Kenneth Brenner dos Anjos Benicio                            %
%   Date:    19/02/2022                                                   %
%   E-mail:  kenneth@gtel.ufc.br                                          %
% ----------------------------------------------------------------------- %
% O conjunto de dados foi carregado inicialmente em um objeto do tipo
% "tabela" e para que seja possível manipulá-lo numericamente é necessário 
% utilizar a função table2array() nativa do MATLAB.


%% ---------- PARTE 1 ---------- %%
disp('%% ---------- PARTE 1 ---------- %% ')

% Transformando a tabela em um conjunto numérico de dados:
ionosphere_array = table2array(ionosphere(:,1:34));

% Matrizes de covariancia considerando os tres metodos sugeridos para a matriz de correlação:
covm_nao_matricial = covariancia.nao_matricial(ionosphere_array.');
covm_matricial = covariancia.matricial(ionosphere_array.');
covm_recursivo = covariancia.recursivo(ionosphere_array.');

% Matriz de covariancia utilizando o comando nativo do MATLAB:
covm_nativo = cov(ionosphere_array);

% Norma das matrizes de diferença/erro:

E_nao_matricial = norm(covm_nao_matricial - covm_nativo);
Z = sprintf('Erro para o método não matricial: %d', E_nao_matricial);
disp(Z)

E_matricial = norm(covm_matricial - covm_nativo);
Z = sprintf('Erro para o método matricial: %d', E_matricial);
disp(Z)

E_recursivo = norm(covm_recursivo - covm_nativo);
Z = sprintf('Erro para o método recursivo: %d', E_recursivo);
disp(Z)


%% ---------- PARTE 2 ---------- %%
disp('%% ---------- PARTE 2 ---------- %% ')

tempo_nao_matricial = 0;
tempo_matricial = 0;
tempo_recursivo = 0;
tempo_nativo = 0;

% Esse fluxo foi criado com a intenção de obter um comportamento médio
% considerando 10000 realizações do experimento pedido. Dessa forma, será
% possível apontar com maior segurança o método de estimação da matriz de
% covariancia que de fato demanda menor tempo computacional.
observacoes = 10000;
for i = 1:observacoes
    
    tic()
    covm_nao_matricial = covariancia.nao_matricial(ionosphere_array.');
    aux = toc();
    tempo_nao_matricial = tempo_nao_matricial + aux;

    tic()
    covm_matricial = covariancia.matricial(ionosphere_array.');
    aux = toc();
    tempo_matricial = tempo_matricial + aux;
    
    tic()
    covm_recursivo = covariancia.recursivo(ionosphere_array.');
    aux = toc();
    tempo_recursivo = tempo_recursivo + aux;
    
    tic()
    covm_nativo = cov(ionosphere_array);
    aux = toc();
    tempo_nativo = tempo_nativo + aux;
    
end

tempo_nao_matricial = tempo_nao_matricial/observacoes
tempo_matricial = tempo_matricial/observacoes
tempo_recursivo = tempo_recursivo/observacoes
tempo_nativo = tempo_nativo/observacoes



%% ---------- PARTE 3 ---------- %% 
disp('%% ---------- PARTE 3 e 4 ---------- %% ')

% Separando o conjunto de dados nas duas classes disponíveis:
ionosphere_class_g = ionosphere(ionosphere.g=='g',:);
ionosphere_class_b = ionosphere(ionosphere.g=='b',:);

% Transformando as tabelas em conjuntos de dados numéricos:
class_g = table2array(ionosphere_class_g(:,1:34));
class_b = table2array(ionosphere_class_b(:,1:34));

% Escolhendo o método matricial para o calculo das matrizes de covariancia
% locais, visto que este apresente o melhor desempenho:
covmtx_g = covariancia.matricial(class_g.');
covmtx_b = covariancia.matricial(class_b.');
covmtx_global = covariancia.matricial(ionosphere_array.');


%% ---------- PARTE 4 ---------- %% 

% Analisando a invertibilidade das matrizes de covariancia locais e global:

rank_class_g = rank(covmtx_g);
condicionamento_class_g = rcond(covmtx_g);
Z = sprintf('O rank da matriz de covariancia local da classe g é %d e o numero de condicionamento é %d', rank_class_g, condicionamento_class_g);
disp(Z)

rank_class_b = rank(covmtx_b);
condicionamento_class_b = rcond(covmtx_b);
Z = sprintf('O rank da matriz de covariancia local da classe b é %d e o numero de condicionamento é %d', rank_class_b, condicionamento_class_b);
disp(Z)

rank_global = rank(covmtx_global);
condicionamento_global = rcond(covmtx_global);
Z = sprintf('O rank da matriz de covariancia global é %d e o numero de condicionamento é %d', rank_global, condicionamento_global);
disp(Z)

