function [L,C] = cluster_kmeans(X,k,varargin)

  my_select = 2;
  
  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end
  
  switch my_select
    case 1
      % function [label,m] = litekmeans(X, k)
      % Perform k-means clustering.
      %   X: d x n data matrix
      %   k: number of seeds
      % Written by Michael Chen (sth4nth@gmail.com).
      n = size(X,2);
      last = 0;
      label = ceil(k*rand(1,n));  % random initialization
      while any(label ~= last)
          [u,~,label] = unique(label);   % remove empty clusters
          label = label(:)';
          k = length(u);
          E = sparse(1:n,label,1,n,k,n);  % transform label into indicator matrix
          m = X*(E*spdiags(1./sum(E,1)',0,k,k));    % compute m of each cluster
          last = label;
          [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1); % assign samples to the nearest centers
      end
      [~,~,label] = unique(label);

    case 2
  %     function [L,C] = cluster_kmeans(X,k)
      % KMEANS Cluster multivariate data using the k-means++ algorithm.
      %   [L,C] = kmeans(X,k) produces a 1-by-size(X,2) vector L with one class
      %   label per column in X and a size(X,1)-by-k matrix C containing the
      %   centers corresponding to each class.
      % 
      %   Version: 2013-02-08
      %   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
      % 
      %   References:
      %   [1] J. B. MacQueen, "Some Methods for Classification and Analysis of 
      %       MultiVariate Observations", in Proc. of the fifth Berkeley
      %       Symposium on Mathematical Statistics and Probability, L. M. L. Cam
      %       and J. Neyman, eds., vol. 1, UC Press, 1967, pp. 281-297.
      %   [2] D. Arthur and S. Vassilvitskii, "k-means++: The Advantages of
      %       Careful Seeding", Technical Report 2006-13, Stanford InfoLab, 2006.

      L = [];
      L1 = 0;

      while length(unique(L)) ~= k

      %     The k-means++ initialization.
          C = X(:,1+round(rand*(size(X,2)-1)));
          L = ones(1,size(X,2));
          for i = 2:k
              D = X-C(:,L);
              D = cumsum(sqrt(sum(D.^2,1)));
%               D = cumsum(sqrt(dot(D,D,1)));
%               fprintf('%g\n',norm(D2-D));
              if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
              C(:,i) = X(:,find(rand < D/D(end),1));
%               [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'));
              [~,L] = max(bsxfun(@minus,2*(C'*X),sum(C.^2,1).'));
          end

      %     The k-means algorithm.
          while any(L ~= L1)
              L1 = L;
%               lM = bsxfun(@eq,[1:k]',L);
%               lMsum = sum(lM,2);
              for i = 1:k
                l = L==i;
                C(:,i) = sum(X(:,l),2)/sum(l);
%                 C(:,i) = sum(X(:,lM(i,:)),2)/lMsum(i);
              end
              [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'),[],1);
          end

      end
    case 3
  %     function [label, centroid, dis] = fkmeans(X, k, options)
      % FKMEANS Fast K-means with optional weighting and careful initialization.
      % [L, C, D] = FKMEANS(X, k) partitions the vectors in the n-by-p matrix X
      % into k (or, rarely, fewer) clusters by applying the well known batch
      % K-means algorithm. Rows of X correspond to points, columns correspond to
      % variables. The output k-by-p matrix C contains the cluster centroids. The
      % n-element output column vector L contains the cluster label of each
      % point. The k-element output column vector D contains the residual cluster
      % distortions as measured by total squared distance of cluster members from
      % the centroid.
      %
      % FKMEANS(X, C0) where C0 is a k-by-p matrix uses the rows of C0 as the
      % initial centroids instead of choosing them randomly from X.
      %
      % FKMEANS(X, k, options) allows optional parameter name/value pairs to 
      % be specified. Parameters are:
      %
      %   'weight' - n-by-1 weight vector used to adjust centroid and distortion
      %              calculations. Weights should be positive.
      %   'careful' - binary option that determines whether "careful seeding"
      %               as recommended by Arthur and Vassilvitskii is used when
      %               choosing initial centroids. This option should be used
      %               with care because numerical experiments suggest it may
      %               be counter-productive when the data is noisy.
      %
      % Notes
      % (1) The careful seeding procedure chooses the first centroid at random
      % from X, and each successive centroid from the remaining points according
      % to the categorical distribution with selection probabilities proportional
      % to the point's minimum squared Euclidean distance from the already chosen
      % centroids. This tends to spread the points out more evenly, and, if the
      % data is made of k well separated clusters, is likely to choose an initial
      % centroid from each cluster. This can speed convergence and reduce the
      % likelihood of getting a bad solution [1]. However, in experiments where
      % 5% uniformly distributed noise data was added to such naturally clustered
      % data the results were frequently worse then when centroids were chosen at
      % random.
      % (2) If, as is possible, a cluster is empty at the end of an iteration,
      % then there may be fewer than k clusters returned. In practice this seems
      % to happen very rarely.
      % (3) Unlike the Mathworks KMEANS this implementation does not perform a
      % final, slow, phase of incremental K-means ('onlinephase') that guarantees
      % convergence to a local minimum. 
      %
      % References
      % [1] "k-means++: The Advantages of Careful Seeding", by David Arthur and
      % Sergei Vassilvitskii, SODA 2007.

      n = size(X,1);
      options = [];

      % option defaults
      weight = 0; % uniform unit weighting
      careful = 0;% random initialization

      if nargin == 3
          if isfield(options, 'weight')
              weight = options.weight;
          end
          if isfield(options,'careful')
              careful = options.careful;
          end
      end

      % If initial centroids not supplied, choose them
      if isscalar(k)
          % centroids not specified
          if careful
              k = spreadseeds(X, k);
          else
              k = X(randsample(size(X,1),k),:);
          end
      end

      % generate initial labeling of points
      [~,label] = max(bsxfun(@minus,k*X',0.5*sum(k.^2,2)));
      k = size(k,1);

      last = 0;

      if ~weight
          % code defactoring for speed
          while any(label ~= last)
              % remove empty clusters
              [~,~,label] = unique(label);
              % transform label into indicator matrix
              ind = sparse(label,1:n,1,k,n,n);
              % compute centroid of each cluster
              centroid = (spdiags(1./sum(ind,2),0,k,k)*ind)*X;
              % compute distance of every point to each centroid
              distances = bsxfun(@minus,centroid*X',0.5*sum(centroid.^2,2));
              % assign points to their nearest centroid
              last = label;
              [~,label] = max(distances);
          end
          dis = ind*(sum(X.^2,2) - 2*max(distances)');
      else
          while any(label ~= last)
              % remove empty clusters
              [~,~,label] = unique(label);
              % transform label into indicator matrix
              ind = sparse(label,1:n,weight,k,n,n);
              % compute centroid of each cluster
              centroid = (spdiags(1./sum(ind,2),0,k,k)*ind)*X;
              % compute distance of every point to each centroid
              distances = bsxfun(@minus,centroid*X',0.5*sum(centroid.^2,2));
              % assign points to their nearest centroid
              last = label;
              [~,label] = max(distances);
          end
          dis = ind*(sum(X.^2,2) - 2*max(distances)');
      end
      label = label';

      % Code below this line reused from the file exchange submission K-means++
      % (http://www.mathworks.com/matlabcentral/fileexchange/28901-k-means) in
      % accordance with the license:
      % vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      % Copyright (c) 2010, Michael Chen
      % All rights reserved.
      % 
      % Redistribution and use in source and binary forms, with or without
      % modification, are permitted provided that the following conditions are
      % met:
      % 
      %     * Redistributions of source code must retain the above copyright 
      %       notice, this list of conditions and the following disclaimer.
      %     * Redistributions in binary form must reproduce the above copyright
      %       notice, this list of conditions and the following disclaimer in the
      %       documentation and/or other materials provided with the distribution
      %       
      % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
      % IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
      % THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
      % PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
      % CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
      % EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
      % PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
      % PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
      % LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
      % NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
      % SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
      % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      % 
  end
  function D = sqrdistance(A, B)
    % Square Euclidean distances between all sample pairs
    % A:  n1 x d data matrix
    % B:  n2 x d data matrix
    % WB: n2 x 1 weights for matrix B
    % D: n2 x n1 pairwise square distance matrix
    %    D(i,j) is the squared distance between A(i,:) and B(j,:)
    % Written by Michael Chen (sth4nth@gmail.com). July 2009.
    n1 = size(A,1); n2 = size(B,2);
    m = (sum(A,1)+sum(B,1))/(n1+n2);
    A = bsxfun(@minus,A,m);
    B = bsxfun(@minus,B,m);
    D = full((-2)*(A*B'));
    D = bsxfun(@plus,D,full(sum(B.^2,2))');
    D = bsxfun(@plus,D,full(sum(A.^2,2)))';
  end

  function [S, idx] = spreadseeds(X, k)
    % X: n x d data matrix
    % k: number of seeds
    % reference: k-means++: the advantages of careful seeding.
    % by David Arthur and Sergei Vassilvitskii
    % Adapted from softseeds written by Mo Chen (mochen@ie.cuhk.edu.hk), 
    % March 2009.
    [n,d] = size(X);
    idx = zeros(k,1);
    S = zeros(k,d);
    D = inf(n,1);
    idx(1) = ceil(n.*rand);
    S(1,:) = X(idx(1),:);
    for i = 2:k
      D = min(D,sqrdistance(S(i-1,:),X));
      idx(i) = find(cumsum(D)/sum(D)>rand,1);
      S(i,:) = X(idx(i),:);
    end
  end
end

function y = randsample(n,k)

  y = randperm(n);
  y = y(1:k);

end