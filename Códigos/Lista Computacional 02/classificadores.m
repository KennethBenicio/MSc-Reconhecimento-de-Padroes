% ----------------------------------------------------------------------- %
%   Version: 1.0                                                          %
%   Author:  Kenneth Brenner dos Anjos Benicio                            %
%   Date:    19/02/2022                                                   %
%   E-mail:  kenneth@gtel.ufc.br                                          %
% ----------------------------------------------------------------------- %
classdef classificadores
    methods(Static)        
        %% Classificador Gaussiano de Log-Verossimilhança Variante 1 (EQ.12)
        %gi (xn) = − Qi (xn) − ln |Σi| + ln p(ωi).
        function [STATS, cfsmtx, acr, TX_OK, X, m, S, posto] = CGQ12(Dados,Nr,Ptrain,PCA)
            %
            % Mahalanobis with one COV matrix per class.  
            %
            % INPUTS: * data (matrix): dataset matrix (N x (p+1))
            %	  	OBS1: feature vectors along the rows of data matrix
            %	  	OBS2: last column must contain numerical class labels
            %	  * Nr (scalar): Number of runs (Nr>=1)
            %	  * Ptrain (scalar): Percentage of training data (0 < Ptrain < 100)
            %
            % OUTPUTS: X (struct) - the data samples separated per class
            %          m (struct) - the classes centroids
            %          S (struct) - the COV matrices per class
            %          STATS (vector) - Statistics of test data (mean, median, min/max, sd)
            %
            % Author: Guilherme Barreto
            % Date: 21/10/2018    
            
            % Principal Component Analysis
            if PCA == 'y'  
                labels = Dados(:,end);
                data_only = Dados(:,1:(end-1));

                Cx = cov(data_only);
                [V, ~, Contribution]= pcacov(Cx);
                %pcacov already sorts out the columns of V;
                %[V, L]= eig(Cx);
                %[L, Positions] = sort(L, 'descend');
                %V = V(:,Positions);

                tol = 0.5;
                aux = cumsum(Contribution)/sum(Contribution);
                components = length(find(aux<tol));
                V = V(:,1:components);

                data_only = (V.')*(data_only.');
                Dados = [data_only.' labels];  
            end
            
            [N, p]=size(Dados);  % Get dataset size (N)

            Ntrn=round(Ptrain*N/100);  % Number of training samples
            Ntst=N-Ntrn; % Number of testing samples

            K=max(Dados(:,end)); % Get the number of classes
            TX_OK = zeros(Nr,1);
            for r = 1:Nr  % Loop of independent runs

              I=randperm(N);
              Dados=Dados(I,:); % Shuffle rows of the data matrix

              % Separate into training and testing subsets
              Dtrn=Dados(1:Ntrn,:);  % Training data
              Dtst=Dados(Ntrn+1:N,:); % Testing data

              % Partition of training data into K subsets
              lamb=0.05;  % regularization parameter
              for k = 1:K
                I=find(Dtrn(:,end)==k);  % Find rows with samples from k-th class
                X{k}=Dtrn(I,1:end-1); % Data samples from k-th class
                m{k}=mean(X{k})';   % Centroid of the k-th class
                S{k}=cov(X{k})+lamb*eye(p-1); % Regularized cov matrix of the k-th class
                posto{k}=rank(S{k}); % Check invertibility of covariance matrix by its rank
                iS{k}=pinv(S{k});    % Inverse covariance matrix of the k-th class
                %iS{k}=inv(S{k});    % Inverse covariance matrix of the k-th class
              end
              
              % À priori probabilities
              aux = Dados(:,end);
              for k = 1:K
                apriori(k) = (sum(aux == k)-1)/N;
              end
              
              % Testing phase
              dist = zeros(K,1);
              correct=0;  % number correct classifications
              for i = 1:Ntst
                Xtst=Dtst(i,1:end-1)';   % test sample to be classified
                Label_Xtst=Dtst(i,end);   % Actual label of the test sample
                for k = 1:K
                  % Mahalanobis distance to k-th class
                  dist(k) = -0.5*((Xtst-m{k}).'*iS{k}*(Xtst-m{k})) - 0.5*(log(det(S{k}))) + log(apriori(k));
                end
                [~, Pred_class]=max(dist);  % index of the minimum distance class
                Label_Hat(i,1) = Pred_class; %#ok<*AGROW>
                
                if Pred_class == Label_Xtst
                    correct=correct+1;
                end
              end
              Label_True = Dtst(:,end); 
              cfsmtx{r}  = confusionmat(Label_True, Label_Hat);
              TX_OK(r)   = 100*correct/Ntst;   % Recognition rate of r-th run
            end
            
            %Accuracy
            for r = 1:Nr
                mtx = cfsmtx{r};
                acr(r) = sum(diag(mtx))/(sum(sum(mtx)));
            end
            acr = mean(acr);
            
            [~, k1] = min(TX_OK);
            [~, k2] = max(TX_OK);
            cfsmtx = {cfsmtx{k1} cfsmtx{k2}};
            STATS  = [mean(TX_OK) std(TX_OK) median(TX_OK) min(TX_OK) max(TX_OK)];
        end
           
        %% Classificador Gaussiano de Log-Verossimilhança Variante 2 (EQ.44)
        function [STATS, cfsmtx, acr, TX_OK, X, m, S, posto] = CGQ44(Dados,Nr,Ptrain,PCA)
            %
            % Mahalanobis with one COV matrix per class.
            %
            % INPUTS: * data (matrix): dataset matrix (N x (p+1))
            %	  	OBS1: feature vectors along the rows of data matrix
            %	  	OBS2: last column must contain numerical class labels
            %	  * Nr (scalar): Number of runs (Nr>=1)
            %	  * Ptrain (scalar): Percentage of training data (0 < Ptrain < 100)
            %
            % OUTPUTS: X (struct) - the data samples separated per class
            %          m (struct) - the classes centroids
            %          S (struct) - the COV matrices per class
            %          STATS (vector) - Statistics of test data (mean, median, min/max, sd)
            %
            % Author: Guilherme Barreto
            % Date: 21/10/2018

            if PCA == 'y'  
                labels = Dados(:,end);
                data_only = Dados(:,1:(end-1));

                Cx = cov(data_only);
                [V, ~, Contribution]= pcacov(Cx);
                %pcacov already sorts out the columns of V;
                %[V, L]= eig(Cx);
                %[L, Positions] = sort(L, 'descend');
                %V = V(:,Positions);

                tol = 0.5;
                aux = cumsum(Contribution)/sum(Contribution);
                components = length(find(aux<tol));
                V = V(:,1:components);

                data_only = (V.')*(data_only.');
                Dados = [data_only.' labels];  
            end
            
            [N, p]=size(Dados);  % Get dataset size (N)

            Ntrn=round(Ptrain*N/100);  % Number of training samples
            Ntst=N-Ntrn; % Number of testing samples

            K=max(Dados(:,end)); % Get the number of classes

            %ZZ=sprintf('The problem has %d classes',K);
            %disp(ZZ);
            TX_OK = zeros(Nr,1);
            for r = 1:Nr  % Loop of independent runs

              I=randperm(N);
              Dados=Dados(I,:); % Shuffle rows of the data matrix

              % Separate into training and testing subsets
              Dtrn=Dados(1:Ntrn,:);  % Training data
              Dtst=Dados(Ntrn+1:N,:); % Testing data

              % Partition of training data into K subsets
              lamb=0.05;  % regularization parameter
              for k = 1:K
                I=find(Dtrn(:,end)==k);  % Find rows with samples from k-th class
                X{k}=Dtrn(I,1:end-1); % Data samples from k-th class
                m{k}=mean(X{k})';   % Centroid of the k-th class
                S{k}=cov(X{k})+lamb*eye(p-1); % Regularized cov matrix of the k-th class
                posto{k}=rank(S{k}); % Check invertibility of covariance matrix by its rank
                iS{k}=pinv(S{k});    % Inverse covariance matrix of the k-th class
                %iS{k}=inv(S{k});    % Inverse covariance matrix of the k-th class
              end

              % Testing phase
              dist = zeros(K,1);
              correct=0;  % number correct classifications
              for i = 1:Ntst
                Xtst=Dtst(i,1:end-1)';   % test sample to be classified
                Label_Xtst=Dtst(i,end);   % Actual label of the test sample
                for k = 1:K 
                  % Mahalanobis distance to k-th class
                  %dist(k) = (Xtst-m{k}).'*iS{k}*(Xtst-m{k}) - 0.5*(log(det(S{k})));
                  dist(k) = -0.5*(Xtst.'*iS{k}*Xtst)+(m{k}.'*iS{k}*Xtst)-0.5*(m{k}.'*iS{k}*m{k})-0.5*(log(det(S{k})));
                end
                [~, Pred_class]=max(dist);  % index of the minimum distance class
                Label_Hat(i,1) = Pred_class; %#ok<*AGROW>
                
                if Pred_class == Label_Xtst
                    correct=correct+1;
                end
              end
              Label_True = Dtst(:,end); 
              cfsmtx{r}  = confusionmat(Label_True, Label_Hat);
              TX_OK(r)   = 100*correct/Ntst;   % Recognition rate of r-th run
            end
            
            %Accuracy
            for r = 1:Nr
                mtx = cfsmtx{r};
                acr(r) = sum(diag(mtx))/(sum(sum(mtx)));
            end
            acr = mean(acr);
            
            [~, k1] = min(TX_OK);
            [~, k2] = max(TX_OK);
            cfsmtx = {cfsmtx{k1} cfsmtx{k2}};
            STATS  = [mean(TX_OK) std(TX_OK) median(TX_OK) min(TX_OK) max(TX_OK)];
        end
        
        %% Classificador Gausiano Quadrático com Matriz de Covariância Agregada Variante 1 (EQ.17)   
        function [STATS, cfsmtx, acr, TX_OK, X, m, S, posto] = CGQ17(Dados,Nr,Ptrain,PCA)
            %
            % Mahalanobis with a common pooled covariance matrix.
            %
            % INPUTS: * data (matrix): dataset matrix (N x (p+1))
            %	  	OBS1: feature vectors along the rows of data matrix
            %	  	OBS2: last column must contain numerical class labels
            %	  * Nr (scalar): Number of runs (Nr>=1)
            %	  * Ptrain (scalar): Percentage of training data (0 < Ptrain < 100)
            %
            % OUTPUTS: X (struct) - the data samples separated per class
            %          m (struct) - the classes centroids
            %          S (struct) - the COV matrices per class
            %          STATS (vector) - Statistics of test data (mean, median, min/max, sd)
            %
            % Author: Guilherme Barreto
            % Date: 21/10/2018
            
            if PCA == 'y'  
                labels = Dados(:,end);
                data_only = Dados(:,1:(end-1));

                Cx = cov(data_only);
                [V, ~, Contribution]= pcacov(Cx);
                %pcacov already sorts out the columns of V;
                %[V, L]= eig(Cx);
                %[L, Positions] = sort(L, 'descend');
                %V = V(:,Positions);

                tol = 0.5;
                aux = cumsum(Contribution)/sum(Contribution);
                components = length(find(aux<tol));
                V = V(:,1:components);

                data_only = (V.')*(data_only.');
                Dados = [data_only.' labels];  
            end
            
            [N, p]=size(Dados);  % Get dataset size (N)

            Ntrn=round(Ptrain*N/100);  % Number of training samples
            Ntst=N-Ntrn; % Number of testing samples

            K=max(Dados(:,end)); % Get the number of classes
            
            %ZZ=sprintf('The problem has %d classes',K);
            %disp(ZZ);

            for r = 1:Nr  % Loop of independent runs

              I=randperm(N);
              Dados=Dados(I,:); % Shuffle rows of the data matrix

              % Separate into training and testing subsets
              Dtrn=Dados(1:Ntrn,:);  % Training data
              Dtst=Dados(Ntrn+1:N,:); % Testing data

              % Partition of training data into K subsets
              Spool=zeros(p-1);
              for k = 1:K
                I=find(Dtrn(:,end)==k);  % Find rows with samples from k-th class
                n{k}=length(I);          % number of samples of the k-th class
                X{k}=Dtrn(I,1:end-1);    % Data samples from k-th class
                m{k}=mean(X{k})';        % Centroid of the k-th class
                S{k}=(n{k}-1)*cov(X{k});          % Covariance matrix of the k-th class
                posto{k}=rank(S{k});    % Check invertibility of covariance matrix by its rank
                Spool = Spool + S{k};
              end
              
              lamb   = 0.05;
              Spool  = Spool/(Ntrn-K);  % Pooled covariance matrix
              Spool  = Spool + lamb*eye(p-1);
              iSpool = pinv(Spool);  % inverse of Spool
              
%               Spool  = Spool/(Ntrn-K);  % Pooled covariance matrix
%               iSpool = pinv(Spool);  % inverse of Spool

              % Testing phase
              correct=0;  % number correct classifications
              for i = 1:Ntst
                Xtst=Dtst(i,1:end-1)';   % test sample to be classified
                Label_Xtst=Dtst(i,end);   % Actual label of the test sample
                for k = 1:K
                  % Mahalanobis distance to k-th class
                  dist(k) = -0.5*(Xtst - m{k}).'*iSpool*(Xtst - m{k});
                end
                [~, Pred_class]=max(dist);  % index of the minimum distance class
                Label_Hat(i,1) = Pred_class; %#ok<*AGROW>
                
                if Pred_class == Label_Xtst
                    correct=correct+1;
                end
              end
              Label_True = Dtst(:,end); 
              cfsmtx{r}  = confusionmat(Label_True, Label_Hat);
              TX_OK(r)   = 100*correct/Ntst;   % Recognition rate of r-th run
            end
            
            %Accuracy
            for r = 1:Nr
                mtx = cfsmtx{r};
                acr(r) = sum(diag(mtx))/(sum(sum(mtx)));
            end
            acr = mean(acr);
            
            [~, k1] = min(TX_OK);
            [~, k2] = max(TX_OK);
            cfsmtx = {cfsmtx{k1} cfsmtx{k2}};
            STATS  = [mean(TX_OK) std(TX_OK) median(TX_OK) min(TX_OK) max(TX_OK)];
        end
        
        %% Classificador Gausiano Quadrático com Matriz de Covariância Agregada Variante 2 (EQ.39)
        %gi (xn) = mi^T Σ⁻¹ xn − mi^T Σ⁻¹ mi .   
        function [STATS, cfsmtx, acr, TX_OK, X, m, S, posto] = CGQ39(Dados,Nr,Ptrain,PCA)
            %
            % Mahalanobis with a common pooled covariance matrix.
            %
            % INPUTS: * data (matrix): dataset matrix (N x (p+1))
            %	  	OBS1: feature vectors along the rows of data matrix
            %	  	OBS2: last column must contain numerical class labels
            %	  * Nr (scalar): Number of runs (Nr>=1)
            %	  * Ptrain (scalar): Percentage of training data (0 < Ptrain < 100)
            %
            % OUTPUTS: X (struct) - the data samples separated per class
            %          m (struct) - the classes centroids
            %          S (struct) - the COV matrices per class
            %          STATS (vector) - Statistics of test data (mean, median, min/max, sd)
            %
            % Author: Guilherme Barreto
            % Date: 21/10/2018

            if PCA == 'y'  
                labels = Dados(:,end);
                data_only = Dados(:,1:(end-1));

                Cx = cov(data_only);
                [V, ~, Contribution]= pcacov(Cx);
                %pcacov already sorts out the columns of V;
                %[V, L]= eig(Cx);
                %[L, Positions] = sort(L, 'descend');
                %V = V(:,Positions);

                tol = 0.5;
                aux = cumsum(Contribution)/sum(Contribution);
                components = length(find(aux<tol));
                V = V(:,1:components);

                data_only = (V.')*(data_only.');
                Dados = [data_only.' labels];  
            end
            
            [N, p]=size(Dados);  % Get dataset size (N)

            Ntrn=round(Ptrain*N/100);  % Number of training samples
            Ntst=N-Ntrn; % Number of testing samples

            K=max(Dados(:,end)); % Get the number of classes
            
            %ZZ=sprintf('The problem has %d classes',K);
            %disp(ZZ);

            for r = 1:Nr  % Loop of independent runs

              I=randperm(N);
              Dados=Dados(I,:); % Shuffle rows of the data matrix

              % Separate into training and testing subsets
              Dtrn=Dados(1:Ntrn,:);  % Training data
              Dtst=Dados(Ntrn+1:N,:); % Testing data

              % Partition of training data into K subsets
              Spool=zeros(p-1);
              for k = 1:K
                I=find(Dtrn(:,end)==k);  % Find rows with samples from k-th class
                n{k}=length(I);          % number of samples of the k-th class
                X{k}=Dtrn(I,1:end-1);    % Data samples from k-th class
                m{k}=mean(X{k})';        % Centroid of the k-th class
                S{k}=(n{k}-1)*cov(X{k});          % Covariance matrix of the k-th class
                posto{k}=rank(S{k});    % Check invertibility of covariance matrix by its rank
                Spool = Spool + S{k};
              end
              
              lamb   = 0.05;
              Spool  = Spool/(Ntrn-K);  % Pooled covariance matrix
              Spool  = Spool + lamb*eye(p-1);
              iSpool = inv(Spool);  % inverse of Spool
              
%               Spool  = Spool/(Ntrn-K);  % Pooled covariance matrix
%               iSpool = pinv(Spool);  % inverse of Spool

              % Testing phase
              correct=0;  % number correct classifications
              for i = 1:Ntst
                Xtst=Dtst(i,1:end-1)';   % test sample to be classified
                Label_Xtst=Dtst(i,end);   % Actual label of the test sample
                for k = 1:K
                  % Mahalanobis distance to k-th class
                  dist(k) = (m{k}.'*iSpool*Xtst) - 0.5*(m{k}.'*iSpool*m{k});
                end
                [~, Pred_class]=max(dist);  % index of the minimum distance class
                Label_Hat(i,1) = Pred_class; %#ok<*AGROW>
                
                if Pred_class == Label_Xtst
                    correct=correct+1;
                end
              end
              Label_True = Dtst(:,end); 
              cfsmtx{r}  = confusionmat(Label_True, Label_Hat);
              TX_OK(r)   = 100*correct/Ntst;   % Recognition rate of r-th run
            end
            
            %Accuracy
            for r = 1:Nr
                mtx = cfsmtx{r};
                acr(r) = sum(diag(mtx))/(sum(sum(mtx)));
            end
            acr = mean(acr);
            
            [~, k1] = min(TX_OK);
            [~, k2] = max(TX_OK);
            cfsmtx = {cfsmtx{k1} cfsmtx{k2}};
            STATS  = [mean(TX_OK) std(TX_OK) median(TX_OK) min(TX_OK) max(TX_OK)];
        end
        
        %% Classificador Gaussiano Linear de Mínimos Quadrados
        function [STATS, cfsmtx, acr, TX_OK, W] = LMQ(Dados,Nr,Ptrain,PCA)
            %
            % Euclidean distance to the centroid (mean) of each class
            %
            % INPUTS: * data (matrix): dataset matrix (N x (p+1))
            %	  	OBS1: feature vectors along the rows of data matrix
            %	  	OBS2: last column must contain numerical class labels
            %	  * Nr (scalar): Number of runs (Nr>=1)
            %	  * Ptrain (scalar): Percentage of training data (0 < Ptrain < 100)
            %
            % OUTPUTS: X (struct) - the data samples separated per class
            %          m (struct) - the classes centroids
            %          STATS (vector) - Statistics of test data (mean, median, min/max, sd)
            %
            % Author: Guilherme Barreto
            % Date: 21/10/2018
            
            if PCA == 'y'  
                labels = Dados(:,end);
                data_only = Dados(:,1:(end-1));

                Cx = cov(data_only);
                [V, ~, Contribution]= pcacov(Cx);
                %pcacov already sorts out the columns of V;
                %[V, L]= eig(Cx);
                %[L, Positions] = sort(L, 'descend');
                %V = V(:,Positions);

                tol = 0.5;
                aux = cumsum(Contribution)/sum(Contribution);
                components = length(find(aux<tol));
                V = V(:,1:components);

                data_only = (V.')*(data_only.');
                Dados = [data_only.' labels];  
            end
            
            [N, p]=size(Dados);  % Get dataset size (N)

            Ntrn=round(Ptrain*N/100);  % Number of training samples
            Ntst=N-Ntrn; % Number of testing samples

            K=max(Dados(:,end)); % Get the number of classes
            
            %ZZ=sprintf('The problem has %d classes',K);
            %disp(ZZ);

            for r = 1:Nr  % Loop of independent runs

              I=randperm(N);
              Dados=Dados(I,:); % Shuffle rows of the data matrix

              % Separate into training and testing subsets
              Dtrn=Dados(1:Ntrn,:);  % Training data
              Dtst=Dados(Ntrn+1:N,:); % Testing data

              % Partition of training data into K subsets
              Xtrn=Dtrn(:,1:end-1); % Input data matrix
              Ytrn=Dtrn(:,end); % Target (numerical labels) data vector

              % Routine to convert numerical labels into binary (1-out-of-K) labels
              Ltrn=zeros(Ntrn,K);  % initialization of label matrix
              Labels=eye(K);
              for k = 1:K
                Ik=find(Ytrn==k);
                nk=length(Ik);
                Ltrn(Ik,:)=repmat(Labels(k,:),nk,1);
              end

              W=pinv(Xtrn)*Ltrn;  % Compute the weight matrix

              % Testing phase
              correct=0;  % number correct classifications
              for i = 1:Ntst
                Xtst=Dtst(i,1:end-1);   % test sample to be classified
                Label_Xtst=Dtst(i,end);   % Actual label of the test sample

                Ypred=Xtst*W;   % Predict the output label (binary)

                [~, Pred_class] = max(Ypred); % Find numerical label
                Label_Hat(i,1) = Pred_class; %#ok<*AGROW>

                if Pred_class == Label_Xtst
                    correct=correct+1;
                end
              end
              Label_True = Dtst(:,end); 
              cfsmtx{r}  = confusionmat(Label_True, Label_Hat);
              TX_OK(r)   = 100*correct/Ntst;   % Recognition rate of r-th run
            end
            
            %Accuracy
            for r = 1:Nr
                mtx = cfsmtx{r};
                acr(r) = sum(diag(mtx))/(sum(sum(mtx)));
            end
            acr = mean(acr);
            
            [~, k1] = min(TX_OK);
            [~, k2] = max(TX_OK);
            cfsmtx = {cfsmtx{k1} cfsmtx{k2}};
                      
            STATS  = [mean(TX_OK) std(TX_OK) median(TX_OK) min(TX_OK) max(TX_OK)];
        end
    end
end
