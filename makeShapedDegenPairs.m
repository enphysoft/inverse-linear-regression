%
% AUTHOR: Albert S. Kim
% Email: albertsk@hawaii.edu
% Affiliation: Civil and Environmental Engineering, Unviersity of Hawaii at Manoa
% DATE: 2021-07-24T22:30:48 HST 
% PURPOSE: As a part of paper entitled "Inverse Sampling of Degenerate Datasets from a Linear Regression Line"
%          submitted to Journal of American Statistical Association 
% FILENAME: makeShapedDegenPairs.m
% 
%% revised at Sun Jun 27 06:41:32
%% 
clc; 
clear all; close all; 
% close all; 
format long
space=' '; 
comma=', '; 
disp("============================================== ")
disp("== Generation of Anscombe's data sets ")

%% Case parameter of a shape function: -1, 0, or +1 
Ssign           = 1 ; 
%% Values of the six main constraint 
N       = 11    ;
Beta1   = 0.5   ;
Xavg    = 9	;
Xvar    = 11	; 
Yavg    = 7.5	;
Yvar    = 4.125	; 
%% 
Beta0   = 3.0           ; % not a main constraint parameter, but to be calculated using Beta1.
m	= (N+mod(N,2))/2; % mid-point index, regardless of N is even or odd. 
Ysum    = Yavg * N      ;
SYY     = Yvar * (N-1)  ;
SXX     = (N-1) * Xvar  ;
Cbias   = Yvar/4        ; 

%% To determin  x-vector with given  with Xavg= Mean(x) and  Xvar = Var(x). 
% For example, Mean(x) = 9 and  Var(x) = 11 
% 
disp(['== For Xvec and Xdel = Xvec - Xavg given Mean(Xvec) = ', num2str(Xavg), ' and  SXX = ', num2str(SXX) ])
a       = sqrt(Xvar)*sqrt(6/N/m)        ; 
X0      = Xavg - a * m                  ;
XvecN   = 1:N                           ;
Xplot   = 1:0.1:N                       ; % fine x-grid 
Xvec    = X0 + a* XvecN                 ;
Xdel    = Xvec - mean(Xvec)             ;

%  Shape function   
disp(['== Generating a shape function ... ' ])
Yvec_fun        = Xvec          ; % to make the same dimension 
Smag            = sqrt(2.e-4)   ; 
Yvec_fun        = Ssign * Smag* (Xvec -4.15).*(Xvec -7.48).*(Xvec -10.71).*(Xvec -13.85)                        ;
Yvec_fun_plot   = Ssign * Smag* (Xplot-4.15).*(Xplot-7.48).*(Xplot-10.71).*(Xplot-13.85)                        ;
Xcur            = min(Xvec) : (max(Xvec)-min(Xvec))/100 : max(Xvec)                                             ; 
Yvec_cur        = Ssign * Smag* (Xcur-4.15).*(Xcur-7.48).*(Xcur-10.71).*(Xcur-13.85) + Beta0 + Beta1* Xcur      ;

if     (Ssign < 0 )
    fileNameOut = 'outputTableNegative.txt'     ;
    txtSign     = 'negative'            ;
elseif (Ssign > 0 )
    fileNameOut = 'outputTablePositive.txt'     ;
    txtSign     = 'positive'            ;
elseif (Ssign == 0 )
    fileNameOut = 'outputTableZero.txt'     ;
    txtSign     = 'zero'                ;
end
fileNameLRG = 'outputLinearReg.txt'         ;

txtOut1a= ['   The ', txtSign, ' sign of a shape function is ' num2str(Ssign),'.'];
txtOut1b= ['   Possible values of variable "Ssign" are: -1, 0, and +1.'];
txtOut2 = ['   Calculated data values are stored in "',fileNameOut,'".'];
txtOut3 = ['   Linear regression results are stored in "',fileNameLRG,'".']; 
disp(txtOut1a)
disp(txtOut1b)
disp(txtOut2)
disp(txtOut3)

%
disp("== Checking Xvec constraints, calcuated using proposed Xvec. ")
Xavg    = mean(Xvec)            ;
Xvar    = var(Xvec)             ;
Yext    = Beta0 + Beta1 * Xvec	;

%% To determin  y-vector with given  with Yavg= Mean(y) and  Yvar = Var(y).
% For example, Mean(y) = 7.5 and  Var(y) = 4.125

%% Main program starts running.
disp('== Generation of an initial y adn Y = y - yavg.')
disp(['== A random vector with components of ', num2str(N), ' is to be created.'])
RootTest = -1           ;
nRootTestCount = 0      ; 
nRootTestMax   = 1000   ;

while (RootTest < 0 )
  %
  rvec		= rand(1,N)	;
  rvec(m)	= 0.5		; 
  SYYcond	= 0		;
  Ybias_init    = Cbias * ( 2 * rvec(1:N) - 1.0 )       ; 
  Yvec_try_base = Beta0 + Beta1 * Xvec + Yvec_fun	; 
  Yvec_try_init = Yvec_try_base + Ybias_init            ;
  %
  Ybias		= Ybias_init                            ;
  Yvec_try	= Yvec_try_init                         ;
  Yavg_try	= mean(Yvec_try)			; 
  Beta1_try	= var (Yvec_try) / var(Xvec)		;
  Beta0_try	= Yavg_try  - Beta1_try * Yavg_try	; 
  Yfit		= Beta0_try + Beta1_try * Xvec		; 
  Ydel_try	= Yvec_try  - Yavg			; 
  beta1_test	= var(Yvec_try)/var(Xvec)		; 
  %
  a10 = Beta1*SXX + sum ((Xdel(N) - Xdel(2  :m-1) ).* Ydel_try(2  :m-1)) ...
                  + sum ((Xdel(N) - Xdel(m+1:N-1) ).* Ydel_try(m+1:N-1))              ;
  a1  = a10  / ( Xdel(1) - Xdel(N) )			;
  b1  = - ( Xdel(N) - Xdel(m) ) / ( Xdel(N) - Xdel(1) )  ;    
  an0 = Beta1*SXX + sum ((Xdel(1) - Xdel(2  :m-1) ).* Ydel_try(2  :m-1)) ...
                  + sum ((Xdel(1) - Xdel(m+1:N-1) ).* Ydel_try(m+1:N-1))              ; 
  an = an0  / ( Xdel(N) - Xdel(1) )			; 
  bn = - ( Xdel(m) - Xdel(1) ) / ( Xdel(N) - Xdel(1) )  ;    
  SYY_tmp = SYY - sum (Ydel_try(2  :m-1).*Ydel_try(2  :m-1)) ...
                - sum (Ydel_try(m+1:N-1).*Ydel_try(m+1:N-1)) ; 
  B		=   ( a1*b1 + an*bn ) / (1 + b1^2 + bn^2)		;
  C		=   ( a1^2  + an^2  ) / (1 + b1^2 + bn^2)		;
  SYY_tmp	=        SYY_tmp      / (1 + b1^2 + bn^2)		;
  Ym0		= -B			;
  Ym1sq		=  SYY_tmp + B^2 - C	;
  Ym1		=  sqrt(Ym1sq)		;
  RootTest	=  sign(Ym1sq)          ; 
  nRootTestCount=  nRootTestCount + 1   ; 

  if (nRootTestCount > nRootTestMax)
      disp("the iteration limit exceeded. Terminating ... ")
      break 
  end 
end

if (nRootTestCount < nRootTestMax)      %
    disp(['== The RootTest passed after ', num2str(nRootTestCount),' trials.'] )
    %% Storing two solutions of Ym 
    YmPlus  = Ym0 + Ym1         ; 
    YmMinus = Ym0 - Ym1		; 
    % for YmPlus
    Ym      = YmPlus            ; 
    Y1	= a1 + b1 * Ym          ;
    Yn	= an + bn * Ym          ;
    Ydel    = Ydel_try          ;
    Ydel(1) = Y1		;
    Ydel(N) = Yn		;
    Ydel(m) = Ym		;
    YvecPlus= Ydel + Yavg	; 
    % for YmMinus
    Ym      = YmMinus           ;
    Y1	= a1 + b1 * Ym          ;
    Yn	= an + bn * Ym          ;
    Ydel    = Ydel_try          ;
    Ydel(1) = Y1		;
    Ydel(N) = Yn		;
    Ydel(m) = Ym		;
    YvecMinus = Ydel + Yavg	; 
    
    %% 
    Yvec_fun_plot       = Ssign * Smag* (Xplot-4.15).*(Xplot-7.48).*(Xplot-10.71).*(Xplot-13.85)    ;
    xplot               = 0:20                    ;
    yplot               = Beta0 + Beta1 * xplot   ;
    figure
    plot(xplot,yplot            ,'-b'  ); hold on
    plot(Xcur ,Yvec_cur         ,"-r"  ); hold on
    plot(Xvec ,Yvec_try_base    ,"or" ,'MarkerSize',2  ); hold on
    plot(Xvec ,Yvec_try         ,"sr" ,'MarkerSize',12  ); hold on
    plot(Xvec ,YvecPlus         ,"ob"  ); hold on
    plot(Xvec ,YvecMinus        ,"^k"  ); hold on
    % 
    legend({'Given linear line ','Shape func.','shaped data','Yvec tried',...
            'Yvec Sol. w/ positive root', 'Yvec Sol. w/ negative root'} ,'Location','southeast')
    
    %% Screen output 
    avgYvecPlus     = mean(YvecPlus)            ;
    avgYvecMinus    = mean(YvecMinus)           ;
    varYvecPlus     = var(YvecPlus)             ;
    varYvecMinus    = var(YvecMinus)            ;
    avgString=[ '== mean(Y): given, plus, minus = ', num2str(Yavg),comma,num2str(avgYvecPlus),comma,num2str(avgYvecMinus)] ; 
    varString=[ '== var (Y): given, plus, minus = ', num2str(Yvar),comma,num2str(varYvecPlus),comma,num2str(varYvecMinus)] ; 
    disp(avgString)
    disp(varString)
    
    %% File output 
    XYpm=[Xvec; Yvec_try_base ; Yvec_try_init ;  YvecPlus ; YvecMinus];
    fileIDout = fopen(fileNameOut,'w');
    headerNames = ["Xvecm", "YvecBase", "YvecInit",  "YvecPlus", "YvecMinus"]; 
    fprintf(fileIDout, '%12s   %12s    %12s    %12s    %12s\n', headerNames );
    fprintf(fileIDout,'%12f   %12f    %12f    %12f    %12f\n',XYpm);
    fclose(fileIDout);
    
    XYplot=[xplot ; yplot];
    fileIDlgr = fopen(fileNameLRG,'w');
    fprintf(fileIDlgr, '       X             Y \n');
    fprintf(fileIDlgr,'%12f   %12f \n',XYplot);
    fclose(fileIDlgr); 
else
    message= ["RootTest did not passed, trying again ..., nRootTestCount = ", num2str(nRootTestCount)]; 
    disp(message)
    disp("Please, re-adjust your parameters of statistical constraints and/or a shape function.")
end
disp("============================================== ")
% End of File

