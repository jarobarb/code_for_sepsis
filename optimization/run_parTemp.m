close all

%parameters
a = 1.5;
b = 1.5;
sigma = 0.5;

%generate random x and find length of x
x = rand(100,1);
x = sort(x);
m = length(x);

% function
y = a*x+b+sigma*rand(100,1);

%use Matlab polyfit to find ideal function
p_polyfit = polyfit(x,y,1);
y_i = polyval(p_polyfit,x);
p = p_polyfit;


%% Call the Parallel Tempering Matlab Code
[cm,acceptm,accattemptm,expected_valuem,nrpmcm,swattemptm,swapm] =...
 parTemp(p,x,y);

%% plots

for mcc = 1:size(cm,1)

  close(figure(mcc)); figure(mcc);

  c = {cm{mcc,:}};
  accept = acceptm(mcc,:);
  accattempt = accattemptm(mcc,:);
  expected_value = expected_valuem{mcc};
  swattempt = swattemptm(mcc,:);
  swap = swapm(mcc,:);

  %calculate acceptance 
  acceptance = accept./accattempt;
  expected_value

  %calculate swap acceptance
  swapacceptance_ratio = swap./swattempt

  %f_mean = mean(f)
  p_polyfit

  % Plot commands:
  %parameter a in cell 1
  subplot(4,1,1)
  plot(c{1}(:,1))
  %parameter a in cell 2
  subplot(4,1,2)
  plot(c{2}(:,1))
  %parameter a in cell 3 
  subplot(4,1,3)
  plot(c{3}(:,1))
  %histogram of parameter a
  subplot(4,1,4)
  histogram(c{1},100)

  %3D histogram of cloud function 
  close(figure(10+mcc)); figure(10+mcc);
  hist3(c{1},[20,20])

  %cloud function, shows relationship between parameters a&b
  close(figure(10+mcc)); figure(10+mcc);
  plot(c{3}(:,1),c{3}(:,2),'.')
  xlabel('a')
  ylabel('b')
  hold on
  plot(c{2}(:,1),c{2}(:,2), '.')
  hold on
  plot(c{1}(:,1),c{1}(:,2),'.')

end