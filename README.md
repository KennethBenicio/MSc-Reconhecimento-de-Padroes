Em cada uma das subpastas existem muitos arquivos que foram gerados durante as simulações, então para evitar maiores confusões explico brevemente aqui a organização dos códigos.

# ~/Lista Computacional 01: 
covariancia.m >>> Módulo de funções para cálculo da covariância.

exercicio.m   >>> Função que gera todos os dados apresentados na LC 01.



# ~/Lista Computacional 02:
classificadores.m >>> Módulo de funções para os diferentes algoritmos de classificação.

classificadores_tabela >>> Código que gera diversas tabelas internamente ao matlab e as converte para latex.

classificadores_classes >>> Código para criação de tabelas de estatísticas descritivas das classes do conjunto de dados (Antes da etapa de classificação).



# ~/Lista Computacional 03:
clusterizacao.m >>> Módulo de funções com algoritmos de clusterização de dados.

clusterizacao_kopt.m >>> Código para a obtenção do valor ótimo para o número de clusters. Apresento aqui dois procedimento: O primeiro, que não foi utilizado no relátorio mas apresenta-se funcional, é baseado no cálculo da moda apos um determinado número de repetições dos algoritmos. O segundo, recomendado pelo professor, é um procedimento simples de Monte Carlo para a obtenção do comportamento médio considerando um determinado número de rodadas.

clusterizacao_tabelas.m >>> Código que gera diversas tabelas internamente ao matlab e as converte para latex.

obs1: Todos os outros arquivos que não foram citados são gerados automaticamente pelas rotinas criadas. 

obs2: Foram utilizados parcialmente os códigos disponibilizados na disciplina.
