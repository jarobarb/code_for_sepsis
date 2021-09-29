function [L] = my_clustering(A)

  L = [];

  %  First make a distance matrix
  dmatrix = ones(size(A,2),size(A,2));
  my_bump = 1e-6;
  for dmrc = 1:size(A,2)
    for dmcc = 1:size(A,2)
      [~,dmatrix(dmrc,dmcc)] = fminbnd(...
        @(t) sum((abs(exp(1i*t)*A(:,dmrc)-A(:,dmcc))).^2),...
        -pi-my_bump,pi+my_bump);
      if dmatrix(dmrc,dmcc) < 0
        keyboard
      end
    end
    fprintf('%g/%g\n',dmrc,size(A,2));
  end
  
  keyboard

end