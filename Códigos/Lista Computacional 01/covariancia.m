% ----------------------------------------------------------------------- %
%   Version: 1.0                                                          %
%   Author:  Kenneth Brenner dos Anjos Benicio                            %
%   Date:    19/02/2022                                                   %
%   E-mail:  kenneth@gtel.ufc.br                                          %
% ----------------------------------------------------------------------- %
classdef covariancia
    methods(Static)
        % CONVENÇÃO DE USO IMPORTANTE: p representa o número de classes, enquanto N o número de observações.
        % Dessa forma, espera-se a entrada de uma matriz de dados pxN pertencente aos reais.
        
        %% Estimação de Rx utilizando a equação não matricial.
        function  Cx = nao_matricial(X)
            [p, N] = size(X);
            Rx = zeros(p);
            m = mean(X,2); %% Média.
            
            for n = 1:N
                Rx = Rx + X(:,n)*(X(:,n).'); %% Matriz de correlação.
            end
            
            Cx = Rx/N - m*(m.'); %% Matriz de covariancia
        end
        
        
        
        %% Estimação de Rx utilizando a equação matricial.
        function Cx = matricial(X)
           [~, N] = size(X);
           m = mean(X,2); %% Média.
           Rx = (X*(X.'))/N; %% Matriz de correlação.
           Cx = Rx - m*(m.'); %% Matriz de covariancia.
        end
        
        
        
        %% Estimação de Rx utilizando a equação recursiva.
        function Cx = recursivo(X)
            [p, N] = size(X);
            m = zeros(p,1);
            Rx = eye(p);
            
            for n = 1:N
                a = (n-1)/n;
                x = X(:,n);
                m = a*m + (1-a)*x;
                Rx = a*Rx + (1-a)*(x*(x.')); %% Matriz de correlação.
                Cx = Rx - m*(m.'); %% Matriz de covariancia.
            end
            
        end
        
        
        
    end
    
end
