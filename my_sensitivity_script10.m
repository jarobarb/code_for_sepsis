%%  This folder completely ignores the 1.94e? data
outlier_matrix_194 = [4,1,1;4,1,2;4,1,3;4,1,4;4,1,5;4,1,6;...
  4,2,1;4,2,2;4,2,3;4,2,4; 4,3,1;4,3,2;4,4,1];
td = 24;
plot_all_data(outlier_matrix_194,'td',td);

%%  Location of paths (as of 2019/05/09)
addpath ../optimization/
addpath ../Finalized_Code/
%  This is for comparing structures...otherwise, not needed.
addpath ../../../../../Matlab_tools/

%%  %  List of many different fits.  Go to the last one.
%%  Original values from V1
% Healthy < Bsource = 128.443 < Aseptic < Bsource = 248.891 < Septic.
% R^2 = -0.701999 (normal space ave); 
% R^2 = -0.00842937 (log-transformed space ave); 
% R^2 = -0.300235 (normal space all); 
% R^2 = 0.173534 (log-transformed space all);
%  Normal 4/10
klabelsdef = {'B_infy','f','kD1','Bmf','nu4'};
kcurr = [5.0046602934939192764
         70.999380840544219495
       0.030009711419666851295
         80.838762887398843304
         1468.9689068779709942]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(end) = 1;
reghig(end) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,600,0,0],...
  'fixedklabels',{'use_log_Bc','clot_capacity','kD3','kD2'}};
usegenalg = false; fdebug_script;

%%  Original values from V1 push td back to 24
% Healthy < Bsource = 128.443 < Aseptic < Bsource = 248.891 < Septic.
% R^2 = -0.701999 (normal space ave); 
% R^2 = -0.00842937 (log-transformed space ave); 
% R^2 = -0.300235 (normal space all); 
% R^2 = 0.173534 (log-transformed space all);
%  Normal 4/10
klabelsdef = {'B_infy','f','kD1','Bmf','nu4'};
kcurr = [4.98880494497202
          65.8337684515266
        0.0151293905016911
          46.3653249934213
          2587.35329191681]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(5) = 1;
reghig(5) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,600,0,0,24],...
  'fixedklabels',{'use_log_Bc','clot_capacity','kD3','kD2','td'}};
usegenalg = false; fdebug_script;

%%  Original values from V1+k1 push td back to 24
% Healthy < Bsource = 128.443 < Aseptic < Bsource = 248.891 < Septic.
% R^2 = -0.701999 (normal space ave); 
% R^2 = -0.00842937 (log-transformed space ave); 
% R^2 = -0.300235 (normal space all); 
% R^2 = 0.173534 (log-transformed space all);
%  Normal 4/10
klabelsdef = {'B_infy','f','kD1','Bmf','nu4','k1'};
kcurr = [0.750292172030321
          187.075386297795
        0.0132536162213115
          555.133399950308
          4000.66030867705
          1.26182291667486]';
kcurr = [0.681953653737825
          199.947235257306
        0.0141469258709592
          820.033777707272
          4001.52197195351
           1.3127433029186]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(5) = 1;
reghig(5) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,600,0,0,24],...
  'fixedklabels',{'use_log_Bc','clot_capacity','kD3','kD2','td'}};
usegenalg = false; fdebug_script;

%%  Original values from V1+k1+k2/kBl push td back to 24
% Healthy < Bsource = 128.443 < Aseptic < Bsource = 248.891 < Septic.
% R^2 = -0.701999 (normal space ave); 
% R^2 = -0.00842937 (log-transformed space ave); 
% R^2 = -0.300235 (normal space all); 
% R^2 = 0.173534 (log-transformed space all);
%  Normal 4/10
klabelsdef = {'B_infy','f','kD1','Bmf','nu4','k1','kBl'};
kcurr = [0.750292172030321
          187.075386297795
        0.0132536162213115
          555.133399950308
          4000.66030867705
          1.26182291667486
          0.6]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(5) = 1;
reghig(5) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,600,0,0,24],...
  'fixedklabels',{'use_log_Bc','clot_capacity','kD3','kD2','td'}};
usegenalg = false; fdebug_script;

%%  V1 plus clot_capacity fit
% Healthy < Bsource = 128.443 < Aseptic < Bsource = 248.891 < Septic.
% R^2 = -0.701999 (normal space ave); 
% R^2 = -0.00842937 (log-transformed space ave); 
% R^2 = -0.300235 (normal space all); 
% R^2 = 0.173534 (log-transformed space all);
%  Normal 4/10
klabelsdef = {'B_infy','f','kD1','Bmf','clot_capacity','nu4'};
kcurr = [4.99749753688595
          68.5727998150189
        0.0236375115255366
          66.9862191611089
           379.28607637041
          1569.42348460266]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(5) = 1;
reghig(5) = 505;
reglow(end) = 1;
reghig(end) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,0,0],...
  'fixedklabels',{'use_log_Bc','kD3','kD2'}};
usegenalg = false; fdebug_script;

%%  V1 plus clot_capacity fit and k1
% Healthy < Bsource = 128.443 < Aseptic < Bsource = 248.891 < Septic.
% R^2 = -0.701999 (normal space ave); 
% R^2 = -0.00842937 (log-transformed space ave); 
% R^2 = -0.300235 (normal space all); 
% R^2 = 0.173534 (log-transformed space all);
%  Normal 4/10
klabelsdef = {'B_infy','f','kD1','Bmf','clot_capacity','nu4','k1'};
kcurr = [0.736964699290478
          213.335614287955
        0.0344697636772473
          1207.32981781487
          407.214725448303
          3640.09871963404
          1.27176887084432]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(5) = 1;
reghig(5) = 505;
reglow(end-1) = 1;
reghig(end-1) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,0,0],...
  'fixedklabels',{'use_log_Bc','kD3','kD2'}};
usegenalg = false; fdebug_script;

%%  V1 plus clot_capacity fit and k1 and td = 24
% Healthy < Bsource = 128.443 < Aseptic < Bsource = 248.891 < Septic.
% R^2 = -0.701999 (normal space ave); 
% R^2 = -0.00842937 (log-transformed space ave); 
% R^2 = -0.300235 (normal space all); 
% R^2 = 0.173534 (log-transformed space all);
%  Normal 4/10
klabelsdef = {'B_infy','f','kD1','Bmf','clot_capacity','nu4','k1'};
kcurr = [0.750292172030321
          187.075386297795
        0.0132536162213115
          555.133399950308
          505
          4000.66030867705
          1.26182291667486]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(5) = 1;
reghig(5) = 505;
reglow(end-1) = 1;
reghig(end-1) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,0,0,24],...
  'fixedklabels',{'use_log_Bc','kD3','kD2','td'}};
usegenalg = false; fdebug_script;

%%  V1 plus clot_capacity fit and k1
% Healthy < Bsource = 128.443 < Aseptic < Bsource = 248.891 < Septic.
% R^2 = -0.701999 (normal space ave); 
% R^2 = -0.00842937 (log-transformed space ave); 
% R^2 = -0.300235 (normal space all); 
% R^2 = 0.173534 (log-transformed space all);
%  Normal 4/10
klabelsdef = {'B_infy','f','kD1','Bmf','clot_capacity','nu4','k1'};
kcurr = [0.756294060199746
          189.165154499672
        0.0139245140081729
          568.636183759993
          242.171418816312
          4001.04015915501
          1.25794062220083]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(5) = 1;
reghig(5) = 505;
reglow(6) = 1;
reghig(6) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,0,0,24],...
  'fixedklabels',{'use_log_Bc','kD3','kD2','td'}};
usegenalg = false; fdebug_script;

%%  V1 plus optimal non-collinear parameters (kind of)
% Healthy < Bsource = 128.443 < Aseptic < Bsource = 248.891 < Septic.
% R^2 = -0.701999 (normal space ave); 
% R^2 = -0.00842937 (log-transformed space ave); 
% R^2 = -0.300235 (normal space all); 
% R^2 = 0.173534 (log-transformed space all);
%  Normal 4/10
klabelsdef = {'kD1','Bmf','clot_capacity','k1','kBl','mul','sA'};
kcurr = [0.0345047683493963
          1209.05580432659
          407.026653687292
          1.27189591880854
         0.600000005554719
       0.00199995116883313
        0.0124137227816369]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(3) = 1;
reghig(3) = 505;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,0,0,0.736962304795435,213.182954867838,3640.20598805358],...
  'fixedklabels',{'use_log_Bc','kD3','kD2','B_infy','f','nu4'}};
usegenalg = false; fdebug_script;

%%  V1 plus clot_capacity fit and k1 and mul
% Healthy < Bsource = 128.443 < Aseptic < Bsource = 248.891 < Septic.
% R^2 = -0.701999 (normal space ave); 
% R^2 = -0.00842937 (log-transformed space ave); 
% R^2 = -0.300235 (normal space all); 
% R^2 = 0.173534 (log-transformed space all);
%  Normal 4/10
klabelsdef = {'B_infy','f','kD1','Bmf','clot_capacity','mul','k1'};
kcurr = [0.728530455732297
          214.397633476279
        0.0327149818566835
          1217.57723740552
          422.227256349954
       0.00200167871174815
          1.27747099703467]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(5) = 1;
reghig(5) = 505;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',false,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,0,0,3640.20598805358],...
  'fixedklabels',{'use_log_Bc','kD3','kD2','nu4'}};
usegenalg = false; fdebug_script;

%%  Allow k1 free, use all data
%  Results in k1 doubling and other features
%  nu4 = 1.  With kD3 = kD2 = 0 (and clot_capacity = 600), we have the 
%  original dosing function.
% Healthy < Bsource = 128.031 < Aseptic < Bsource = 249.28 < Septic.
% R^2 = -0.673486 (normal space ave); 
% R^2 = 0.257739 (log-transformed space ave); 
% R^2 = -0.289693 (normal space all); 
% R^2 = 0.429842 (log-transformed space all);
%  Normal 4/10
klabelsdef = {'B_infy','f','kD1','Bmf','nu4','k1'};
kcurr = [1.22732111531223
          149.346518822541
         0.030170304223655
            337.0898084753
          3773.47098609286
          1.00419190709545]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(end-1) = 1;
reghig(end-1) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,600,0,0],...
  'fixedklabels',{'use_log_Bc','clot_capacity','kD3','kD2'}};
usegenalg = false; fdebug_script;

%%  Allow k1 to be free, use all data, also fit k5
%  Story is k5 wants to get smaller...data says it doesn't want to rely on
%  acute inflammation killing off bacteria.  Bacteria equation effectively
%  becomes decoupled with the rest.
%  nu4 = 1.  With kD3 = kD2 = 0 (and clot_capacity = 600), we have the 
%  original dosing function.
% Healthy < Bsource = 177.945 < Aseptic < Bsource = 177.945 < Septic.
% R^2 = -1.37485 (normal space ave); 
% R^2 = 0.655544 (log-transformed space ave); 
% R^2 = -0.453283 (normal space all); 
% R^2 = 0.711635 (log-transformed space all);
%  Normal 0/10
klabelsdef = {'B_infy','f','kD1','Bmf','nu4','k1','k5'};
kcurr = [0.381807052187948
          585.923541905356
        0.0743111274599988
           9096.4936320112
          194.727187451662
          1.41081087349689
         0.215182602395603]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(end-2) = 1;
reghig(end-2) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,600,0,0,inf],...
  'fixedklabels',{'use_log_Bc','clot_capacity','kD3','kD2','fint'}};
usegenalg = false; fdebug_script;

%%  Allow k1 free, use average data
%  nu4 = 1.  With kD3 = kD2 = 0 (and clot_capacity = 600), we have the 
%  original dosing function.
% Healthy < Bsource = 177.945 < Aseptic < Bsource = 177.945 < Septic.
% R^2 = -1.37485 (normal space ave); 
% R^2 = 0.655544 (log-transformed space ave); 
% R^2 = -0.453283 (normal space all); 
% R^2 = 0.711635 (log-transformed space all);
%  Normal 0/10
klabelsdef = {'B_infy','f','kD1','Bmf','nu4','k1'};
kcurr = [0.573129214034351
          853.528346215018
        0.0681520697763593
          12792.1704539541
          9.50996053260375
          1.41098624059482]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(end-1) = 1;
reghig(end-1) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',false,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,600,0,0,inf],...
  'fixedklabels',{'use_log_Bc','clot_capacity','kD3','kD2','fint'}};
usegenalg = false; fdebug_script;

%%  Allow k1 to be free, use average data, also fit k5
%  Story is k5 wants to get smaller...data says it doesn't want to rely on
%  acute inflammation killing off bacteria.  Bacteria equation effectively
%  becomes decoupled with the rest.
%  nu4 = 1.  With kD3 = kD2 = 0 (and clot_capacity = 600), we have the 
%  original dosing function.
% Healthy < Bsource = 177.945 < Aseptic < Bsource = 177.945 < Septic.
% R^2 = -1.37485 (normal space ave); 
% R^2 = 0.655544 (log-transformed space ave); 
% R^2 = -0.453283 (normal space all); 
% R^2 = 0.711635 (log-transformed space all);
%  Normal 0/10
klabelsdef = {'B_infy','f','kD1','Bmf','nu4','k1','k5'};
kcurr = [0.573129214034351
          853.528346215018
        0.0681520697763593
          12792.1704539541
          9.50996053260375
          1.41098624059482
          0.40776411686464]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(end-1) = 1;
reghig(end-1) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',false,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,600,0,0,inf],...
  'fixedklabels',{'use_log_Bc','clot_capacity','kD3','kD2','fint'}};
usegenalg = false; fdebug_script;

%%  No 194? data. logall. Set kD3 and kD2 to 0. Let nu4 roam.  Linux
%  nu4 = 1.  With kD3 = kD2 = 0 (and clot_capacity = 600), we have the 
%  original dosing function.
% Healthy < Bsource = 128.443 < Aseptic < Bsource = 248.891 < Septic.
% R^2 = -0.701999 (normal space ave); 
% R^2 = -0.00842937 (log-transformed space ave); 
% R^2 = -0.300235 (normal space all); 
% R^2 = 0.173534 (log-transformed space all);
%  Normal 4/10
klabelsdef = {'B_infy','f','kD1','Bmf','clot_capacity','nu4','k1'};
kcurr = [0.870930508025787/2
          517.673444134444
        0.0609104463275865
          11927.5081556635
          250
          1.00088553214953
          2*1.41237598051023]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(end-1) = 1;
reghig(end-1) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',false,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,0,0],...
  'fixedklabels',{'use_log_Bc','kD3','kD2'}};
usegenalg = false; fdebug_script;

%%  No 194? data. logall. Set kD3 and kD2 to 0. Let nu4 roam.  Linux
%  nu4 = 1.  With kD3 = kD2 = 0 (and clot_capacity = 600), we have the 
%  original dosing function.
% Healthy < Bsource = 128.443 < Aseptic < Bsource = 248.891 < Septic.
% R^2 = -0.701999 (normal space ave); 
% R^2 = -0.00842937 (log-transformed space ave); 
% R^2 = -0.300235 (normal space all); 
% R^2 = 0.173534 (log-transformed space all);
%  Normal 4/10
klabelsdef = {'B_infy','f','kD1','Bmf','clot_capacity','nu4','k1'};
kcurr = [0.80236446649422
          197.501694491328
         0.029360877138451
          837.595933409809
          420.655632564131
          3697.72081770904
          1.22628683854467]';
kcurr = [7773.98944133487
          685.508762885341
        0.0522599395302567
          9804.35696664821
          419.546731826337
          1.00807736881186
          1.41041362684352]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
regcoefs(1) = 1e-1;
reglow(1) = 1;
reghig(1) = 200;
regcoefs(5) = 1;
reglow(5) = 10;
reghig(5) = 400;
reglow(end-1) = 1;
reghig(end-1) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',false,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,0,0,true],...
  'fixedklabels',{'use_log_Bc','kD3','kD2','effective_carrying_capacity'}};
usegenalg = false; fdebug_script;

%%  V1 plus clot_capacity fit and k1
% Healthy < Bsource = 128.443 < Aseptic < Bsource = 248.891 < Septic.
% R^2 = -0.701999 (normal space ave); 
% R^2 = -0.00842937 (log-transformed space ave); 
% R^2 = -0.300235 (normal space all); 
% R^2 = 0.173534 (log-transformed space all);
%  Normal 4/10
klabelsdef = {'B_infy','f','kD1','Bmf','clot_capacity','nu4','k1'};
kcurr = [0.736964699290478
          213.335614287955
        0.0344697636772473
          1207.32981781487
          407.214725448303
          3640.09871963404
          1.27176887084432]';
uvals = ones(size(kcurr)+[0 1]);
regcoefs = 10*ones(size(kcurr));
reglow = eps*ones(size(kcurr));
reghig = (1/eps)*ones(size(kcurr));
reglow(5) = 1;
reghig(5) = 505;
reglow(end-1) = 1;
reghig(end-1) = 4000;
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'regularization_coefs',regcoefs,...
  'regularization_lower',reglow,'regularization_upper',reghig,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,0,0],...
  'fixedklabels',{'use_log_Bc','kD3','kD2'}};
usegenalg = false; fdebug_script;

%%  Actual optimization procedures
%  INITIALIZE KLABELSDEF, KCURR, ETC IN SECTION ABOVE FIRST!
usegenalg = false; fdebug_script;
useglosea = false;
val = 10;
mytol = 0.01;
if exist('vals') ~= 1
  vals = [];
  kcurrs = [];
end
numsubplotsrc = ceil(sqrt(numel(kcurr)+1));

optionsc = optimset('MaxIter',50,'TolFun',1e-16,...
  'MaxFunEvals',500,'Display','iter');

dudk = ones(size(kcurr)); kmat = {};
close(figure(1)); figure(1);

fdebug_script;

close(figure(1)); 

while (val > mytol)

  %  Structured sampling plus sensitivity estimates
  for sens_frac_nexp = 1:8
    [dudk(sens_frac_nexp,:),~,~,kmat{sens_frac_nexp},...
      uvals(sens_frac_nexp,:)] = my_sensitivity(kcurr,klabelsdef,...
      addlargs,'util',reg_util,'frac',(-0.1)^sens_frac_nexp,...
      'parcomp',true);
    fprintf('Sensitivity #%d/#%d done.\n',sens_frac_nexp,8);
  end
  [val,a2] = min(uvals(:));
  [b1,b2] = ind2sub(size(uvals),a2);
  knext = kmat{b1}(:,b2)';

  %  Since I have to think too hard to figure out the best dudk if we
  %  have Monte Carlo info, I just take the best one from structured
  %  sampling and use that to perform "steepest descent"
  bestdudk = dudk(b1,:)
  maxmult = 0.1^b1;
  lmult = -maxmult; rmult = maxmult;

  kdir = @(mult) kcurr+mult*bestdudk.*min(abs(kcurr./bestdudk));

  f = @(mult) reg_util(kdir(mult),klabelsdef,addlargs);

  [mult,valsd] = par_opt1d(lmult,0,rmult,f,'xtol',maxmult*1e-6,...
    'ftol',maxmult*1e-6,'abstol',maxmult*1e-1);

  if valsd < val
    knext = kdir(mult);
  end
  
  %  We also try out the nelder-mead simplex method:
  fforsimp = @(expk) reg_util(knext'.*exp(expk-1),klabelsdef,addlargs);
  [~,vertices,fvals,~,~] = par_simplex2((kmat{1}./kmat{1}(:,1))',...
    fforsimp,'tol',maxmult*1e-6);
  [val,a2] = min(fvals);
  expk = vertices(a2,:)';
  knext = (knext'.*exp(expk-1))';
  
  %  We also try Matlab's built-in genetic algorithm
  if usegenalg
    gaopts = optimoptions(@ga,'Display','iter','UseParallel',true,...
      'InitialPopulationMatrix',knext,'FunctionTolerance',maxmult*1e-6,...
      'PlotFcn',{@gaplotbestf,@gaplotstopping});
    [X,fval] = ga(fga,numel(knext),[],[],[],[],zeros(size(knext)),...
      inf*ones(size(knext)),[],gaopts);
    if fval < val
      knext = X;
    end
  end
  
  %  We also try to use Matlab's global search algorithm
  if useglosea
    %  Similar to the simplex method above, exponential spacing.
    fforgs = @(expk) reg_util(knext.*exp(expk),klabelsdef,addlargs);
    gs = GlobalSearch('FunctionTolerance',maxmult*1e-6,'Display','iter',...
      'PlotFcn',{@gsplotbestf,@gsplotfunccount});
    expkini = ones(size(knext));
    problem = createOptimProblem('fmincon','x0',expkini,...
      'objective',fforgs,'lb',log(10^-3*expkini),'ub',log(10^3*expkini));
    [X,fval] = run(gs,problem);
  end

  %  Evaluating each point    
  parcomp = true;
  show_mons = true;
  clear my_jobs;
  mc_samples_per = 8;
  if parcomp
    p = gcp();
    mc_samples_per = p.NumWorkers*max(1,round(mc_samples_per/p.NumWorkers));
  end
    
  %  Monte carlo sampling
  %  Setting up the points so they occupy concentric "circular shells" in
  %  parameter space
  for mc_frac_nexp = 1:14
    mc_pts = [];
    while true
      temp = knext.*(1+0.1^mc_frac_nexp/numel(knext)*...
        randn(size(knext)));
      dev = norm((temp-knext)./knext);
      if (dev < 0.1^mc_frac_nexp) && (dev > 0.1^(mc_frac_nexp+1))
        mc_pts = [mc_pts;temp];
      end
      if size(mc_pts,1) == mc_samples_per
        break
      end
    end
    if parcomp
      for mcc = 1:size(mc_pts,1)
        my_jobs(mcc) = parfeval(p,reg_util,1,mc_pts(mcc,:),klabelsdef,...
          addlargs);
      end
      for mcc = 1:size(mc_pts,1)
        [job_index,valnew] = fetchNext(my_jobs,5);
        if ~isempty(valnew)
          if valnew < val
            val = valnew;
            knext = mc_pts(job_index,:);
            disp('MC FOUND SOMETHING!');
          end
        end
        if show_mons
          fprintf('%d/%d for MC\n',...
            (mc_frac_nexp-1)*mc_samples_per+mcc,...
            14*mc_samples_per);
        end
      end
    else
      for mcc = 1:size(mc_pts,1)
        valnew = reg_util(mc_pts(mcc,:),klabelsdef,addlargs);
        if valnew > val
          val = valnew;
          knext = mc_pts(mcc,:);
          disp('MC FOUND SOMETHING!');
        end
        if show_mons, fprintf('%d/%d for MC\n',mcc,size(mc_pts,1)); end
      end
    end
  end
%   if parcomp, delete(p); pause(1); end
        
  kcurr = knext
  try 
    kcurrs(end+1,:) = kcurr;
  catch
    kcurrs = kcurr;
  end
%   kcurr(1) = min(kcurr(1),10000);
  val
  vals(end+1) = val;
  if ~ishandle(1), figure(1); end
  for spc = 1:numel(kcurr)
    subplot(numsubplotsrc,numsubplotsrc,spc);
    plot(abs(kcurrs(:,spc)/mean(kcurrs(:,spc))-1));
    title(klabelsdef{spc}); set(gca,'YScale','log');
  end
  subplot(numsubplotsrc,numsubplotsrc,numel(kcurr)+1);
  plot(abs(vals/mean(vals)-1));
  
  title('vals'); set(gca,'YScale','log');
  
  if ispc
    save currwin vals kcurrs kcurr klabelsdef uvals addlargs usegenalg;
  else
    save currlin vals kcurrs kcurr klabelsdef uvals addlargs usegenalg;
  end
end