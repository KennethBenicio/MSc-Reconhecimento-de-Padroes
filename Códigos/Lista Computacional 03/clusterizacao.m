% ----------------------------------------------------------------------- %
%   Version: 1.0                                                          %
%   Author:  Kenneth Brenner dos Anjos Benicio                            %
%   Date:    19/02/2022                                                   %
%   E-mail:  kenneth@gtel.ufc.br                                          %
% ----------------------------------------------------------------------- %
classdef clusterizacao
    methods(Static)
        %% Clusterização via K-means Sequential
        function [Icluster, Particao, W, index] = kmseq(Data, K, PCA)        
            % Normalização dos dados
            Data = (Data - mean(Data))./std(Data);
            
            % Principal Component Analysis se desejado
            if PCA == 'y'  
                data_only = Data(:,1:end);
                Cx = cov(data_only);
                [V, ~, Contribution]= pcacov(Cx);
                tol = 0.65;
                aux = cumsum(Contribution)/sum(Contribution);
                components = length(find(aux<tol));
                V = V(:,1:components);
                data_only = (V.')*(data_only.');
                Data = [data_only.'];  
            end
            
            [N, P] = size(Data);

            % Posicao inicial dos prototipos
            I=randperm(N,K); 
            W=Data(I,:);   
            
            Ne=1000;   % Numero de iteracoes
            lr=0.01;  % Passo de aprendizagem

            for r=1:Ne
              % Embaralha conjunto de dados a cada epoca de treinamento
              I=randperm(N); Data=Data(I,:);

              % Busca pelo prototipo mais proximo do vetor de atributos atual
              soma=0;
              for t=1:N
                for i=1:K
                  Dist2(i)=norm(Data(t,:)-W(i,:));
                end
                [~, Istar]=min(Dist2);   % Indice do prototipo mais proximo

                aux=Data(t,:)-W(Istar,:);
                soma=soma+aux*aux';     % Usado no calculo do SSD por epoca

                % Atualiza posicao do prototipo mais proximo
                W(Istar,:) = W(Istar,:) + lr*aux;
              end

              SSD(r)=soma;

            end

            % Realiza particao pos-treinamento do conjunto original em K subconjuntos
            for t=1:N
                for i=1:K
                  Dist2(i)=norm(Data(t,:)-W(i,:));
                end
                [~, Istar]=min(Dist2);   % Indice do prototipo mais proximo

                Icluster(t)=Istar; 
            end

            for k=1:K
                I=find(Icluster==k); % Indice de todos os exemplos mais proximos do prototipo W_k
                Particao{k}=Data(I,:);  % Separa todos os exemplos da k-esima particao
            end
            
            % Indice de Dunn
            % Dispersão interna ao grupo
            interno = [];
            for k= 1:K
                cluster_k     = find(Icluster==k);
                aux = Data(cluster_k,:); %#ok<*FNDSB>
                [I, ~] = size(aux);
                auxiliar2 = zeros(I,1);
                if (I == 1) || (I == 0) || isempty(I)
                    interno(k,1) = 1e-10;
                else
                    for i = 1:I
                        J = 1:I;
                        J(i) = [];
                        z = 0;
                        auxiliar1 = zeros(I-1,1);
                        for j = J
                            z = z + 1;
                            auxiliar1(z,1) = sqrt(sum((aux(i,:) - aux(j,:)).^2,2)); %#ok<*AGROW>
                        end
                        auxiliar2(i) = max(auxiliar1);
                    end
                    interno(k,1) = max(auxiliar2);
                end
            end

            % Dispersão externa ao grupo
            externo = [];
            for k= 1:K
                cluster_k     = find(Icluster==k);
                cluster_not_k = find(Icluster~=k);
                aux1 = Data(cluster_k,:); %#ok<*FNDSB>
                aux2 = Data(cluster_not_k,:); 
                [I, ~] = size(aux1);
                auxiliar = {};
                if (I == 0) || isempty(I)
                    externo(k,1) = 1e+10;
                else
                    for i = 1:I
                        auxiliar{i} = sqrt(sum((aux1(i,:) - aux2).^2,2));
                    end
                    externo(k,1) = min(min(cell2mat(auxiliar)));
                end
            end
            dunn = min(externo)./max(interno);

            % Calculo do indice Davies-Bouldin
            % Dispersão interna ao grupo              
            interno  = []; 
            externo  = []; 
            auxiliar = {};

            for k= 1:K
                soma = 0;
                cluster_k = find(Icluster==k);
                aux = Data(cluster_k,:); %#ok<*FNDSB>
                [I, ~] = size(aux);
                if (I == 0) || isempty(I)
                    interno{k} = zeros(1,6);
                else
                    for i = 1:I
                        soma = soma + (aux(i,:) - W(k,:)).^2;
                    end
                    interno{k} = sqrt(soma./I);
                end
            end

            % Dispersão externa ao grupo
            for k= 1:K    
                AUX = 1:K;
                AUX(k) = [];
                for aux = AUX
                    externo(k,aux) = sqrt(sum((W(k,:) - W(aux,:)).^2));
                end
            end

            % Calculo da razão
            for k = 1:K
                z = 0; 
                AUX = 1:K;
                AUX(k) = [];
                %auxiliar = zeros(K-1,1);
                for aux = AUX
                    z = z + 1;
                    if (externo(k,aux) == 0)
                        auxiliar{z} = 0;
                        continue;
                    else
                        auxiliar{z} = (interno{k}+interno{aux})./externo(k,aux);
                    end
                end
                razao(k) = max(cell2mat(auxiliar));
            end
            db = sum(razao)./K;

            % Calculo do indice Calinski-Harabasz
            % Matriz de dispersão entregrupos Sb
            Sb = zeros(1,1);
            for k = 1:K
                Nk = sum(Icluster == k);
                if (Nk == 0) || isempty(Nk)
                    continue;
                else
                    Sb = Sb + Nk*((W(k,:) - mean(Data)).')*(W(k,:) - mean(Data));  
                end
            end

            % Matriz de dispersão intragrupo Sw
            Sw = zeros(1,1);
            for k = 1:K
                cluster_k_index = find(Icluster == k);
                cluster_k = Data(cluster_k_index,:);
                [Nk, ~] = size(cluster_k);
                if (Nk == 0) || isempty(Nk)
                    continue;
                else
                    Sw = Sw + Nk*cov(cluster_k);
                end
            end
            ch = ((N - K)./(K - 1))*(trace(Sb)./trace(Sw));
            index = [dunn db ch];
 
%             % Mostra evolucao do SSD ao longo das epocas de treinamento
%             figure; 
%             plot(SSD,'linewidth',3); 
%             title('K-Means Sequencial: SSD x Iteração');
%             xlabel('Iteração');
%             ylabel('SSD'); 
%             grid on;
%             ylim tight;
%             set(gca, "fontsize", 8);
             
%             if P == 2
%                 figure;
%                 aux = gscatter(Data(:,1),Data(:,2),Icluster,'bg','o',8);        
%                 % fill markers 
%                 for g = 1:prod(size(aux)) %#ok<PSIZE>
%                     aux(g).MarkerFaceColor = aux(g).Color;
%                 end  
%                 grid on;
%             end     
        end
        
        %% Clusterização via K-means em Batch
        function [Icluster, Particao, W, index] = kmbat(Data, K, PCA)
            % Normalização dos dados 
            Data = (Data - mean(Data))./std(Data);
            
            % Principal Component Analysis se desejado
            if PCA == 'y'  
                data_only = Data(:,1:end);
                Cx = cov(data_only);
                [V, ~, Contribution]= pcacov(Cx);
                tol = 0.65;
                aux = cumsum(Contribution)/sum(Contribution);
                components = length(find(aux<tol));
                V = V(:,1:components);
                data_only = (V.')*(data_only.');
                Data = [data_only.'];  
            end
            
            [N, P]=size(Data);

            % Posicao inicial dos prototipos
            I=randperm(N,K); 
            W=Data(I,:);  
            
            Ne=100;   % Numero de iteracoes

            for r=1:Ne
              % Busca pelo prototipo mais proximo do vetor de atributos
              for t=1:N
                  for k=1:K
                      Dist2(k)=norm(Data(t,:)-W(k,:),2);
                  end
                  [Dmin(t), Icluster(t)]=min(Dist2);   % Indice do prototipo mais proximo
              end  
              SSD(r)=sum(Dmin.^2);     % Soma das distancias quadraticas por iteracao

              % Particiona dados em K subconjuntos e atualiza prototipo correspondente
              for k=1:K
                I=find(Icluster==k); % Indice de todos os exemplos mais proximos do prototipo W_k
                Particao{k}=Data(I,:);  % Separa todos os exemplos da k-esima particao
                if isempty(I)
                    W(k,:) = zeros(1,P);
                else
                    W(k,:) = mean(Particao{k});     % Atualiza posicoes dos prototipos
                end
              end

              % Calcula SSD 
              for t=1:N
                  for k=1:K
                      Dist2(k)=norm(Data(t,:)-W(k,:),2);
                  end
                  [~, Icluster(t)]=min(Dist2);   % Indice do prototipo mais proximo
              end 

            end
            
            % Indice de Dunn
            % Dispersão interna ao grupo
            interno = [];
            for k= 1:K
                cluster_k     = find(Icluster==k);
                aux = Data(cluster_k,:); %#ok<*FNDSB>
                [I, ~] = size(aux);
                auxiliar2 = zeros(I,1);
                if (I == 1) || (I == 0) || isempty(I)
                    interno(k,1) = 1e-10;
                else
                    for i = 1:I
                        J = 1:I;
                        J(i) = [];
                        z = 0;
                        auxiliar1 = zeros(I-1,1);
                        for j = J
                            z = z + 1;
                            auxiliar1(z,1) = sqrt(sum((aux(i,:) - aux(j,:)).^2,2)); %#ok<*AGROW>
                        end
                        auxiliar2(i) = max(auxiliar1);
                    end
                    interno(k,1) = max(auxiliar2);
                end
            end

            % Dispersão externa ao grupo
            externo = [];
            for k= 1:K
                cluster_k     = find(Icluster==k);
                cluster_not_k = find(Icluster~=k);
                aux1 = Data(cluster_k,:); %#ok<*FNDSB>
                aux2 = Data(cluster_not_k,:); 
                [I, ~] = size(aux1);
                auxiliar = {};
                if (I == 0) || isempty(I)
                    externo(k,1) = 1e+10;
                else
                    for i = 1:I
                        auxiliar{i} = sqrt(sum((aux1(i,:) - aux2).^2,2));
                    end
                    externo(k,1) = min(min(cell2mat(auxiliar)));
                end
            end
            dunn = min(externo)./max(interno);

            % Calculo do indice Davies-Bouldin
            % Dispersão interna ao grupo              
            interno  = []; 
            externo  = []; 
            auxiliar = {};

            for k= 1:K
                soma = 0;
                cluster_k = find(Icluster==k);
                aux = Data(cluster_k,:); %#ok<*FNDSB>
                [I, ~] = size(aux);
                if (I == 0) || isempty(I)
                    interno{k} = zeros(1,6);
                else
                    for i = 1:I
                        soma = soma + (aux(i,:) - W(k,:)).^2;
                    end
                    interno{k} = sqrt(soma./I);
                end
            end

            % Dispersão externa ao grupo
            for k= 1:K    
                AUX = 1:K;
                AUX(k) = [];
                for aux = AUX
                    externo(k,aux) = sqrt(sum((W(k,:) - W(aux,:)).^2));
                end
            end

            % Calculo da razão
            for k = 1:K
                z = 0; 
                AUX = 1:K;
                AUX(k) = [];
                %auxiliar = zeros(K-1,1);
                for aux = AUX
                    z = z + 1;
                    if (externo(k,aux) == 0)
                        auxiliar{z} = 0;
                        continue;
                    else
                        auxiliar{z} = (interno{k}+interno{aux})./externo(k,aux);
                    end
                end
                razao(k) = max(cell2mat(auxiliar));
            end
            db = sum(razao)./K;

            % Calculo do indice Calinski-Harabasz
            % Matriz de dispersão entregrupos Sb
            Sb = zeros(1,1);
            for k = 1:K
                Nk = sum(Icluster == k);
                if (Nk == 0) || isempty(Nk)
                    continue;
                else
                    Sb = Sb + Nk*((W(k,:) - mean(Data)).')*(W(k,:) - mean(Data));  
                end
            end

            % Matriz de dispersão intragrupo Sw
            Sw = zeros(1,1);
            for k = 1:K
                cluster_k_index = find(Icluster == k);
                cluster_k = Data(cluster_k_index,:);
                [Nk, ~] = size(cluster_k);
                if (Nk == 0) || isempty(Nk)
                    continue;
                else
                    Sw = Sw + Nk*cov(cluster_k);
                end
            end
            ch = ((N - K)./(K - 1))*(trace(Sb)./trace(Sw));
            index = [dunn db ch];
            
%             % Mostra evolucao do SSD ao longo das epocas de treinamento
%             figure; 
%             plot(SSD,'linewidth',3); 
%             title('K-Means em Batch: SSD x Iteração');
%             xlabel('Iteração');
%             ylabel('SSD'); 
%             grid on;
%             ylim tight
%             set(gca, "fontsize", 12);
             
%             if P == 2
%                 figure;
%                 aux = gscatter(Data(:,1),Data(:,2),Icluster,'bg','o',8);        
%                 % fill markers 
%                 for g = 1:prod(size(aux)) %#ok<PSIZE>
%                     aux(g).MarkerFaceColor = aux(g).Color;
%                 end  
%                 grid on;
%             end
        end
        
    end
end
        