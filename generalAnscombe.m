%
% AUTHOR: Albert S. Kim
% Email: albertsk@hawaii.edu
% Affiliation: Civil and Environmental Engineering, Unviersity of Hawaii at Manoa
% DATE: 2021-07-24T22:30:48 HST 
% PURPOSE: As a part of paper entitled "Inverse Sampling of Degenerate Datasets from a Linear Regression Line"
%          submitted to Journal of American Statistical Association 
% FILENAME: generalAnscombe.m
%
clc; clear all; close all; format long
disp("========================================= ")
disp("== Generation of Anscombe's data sets === ")
disp("========================================= ")
% The number of samples 
N       =  11           ;
Beta0   =  3.0          ;
Beta1   =  0.5          ;
%% 
m       =  (N+1)/2      ;
%%
Xavg 	=  9            ;
Xvar 	=  11           ; 
%%
Yavg 	=  7.5          ;
Yvar 	=  4.125	;

%% 
% 
% For variable x with  Mean(x) = 9 and  Var(x) = 11
% 
disp("== For Xvec and Xdel = Xvec - Xavg given Mean(Xvec) = 9 and  SXX = 110.")
a    	=  sqrt(Xvar)*sqrt(6/N/m)       ;
X0   	=  Xavg - a * m                 ;
Xvec	=  1:N                          ;
Xvec 	=  X0 + a* Xvec                 ;
Xdel 	=  Xvec - mean(Xvec)            ;
disp("   Checking Xvec constraints, calcuated using proposed Xvec.")
Xavg 	=  mean(Xvec)                   ;
Xvar 	=  var(Xvec)                    ;
SXX  	=  (N-1) * Xvar                 ;
Yext 	=  Beta0 + Beta1 * Xvec         ;

% 
% For variable y with  Mean(y) = 7.5 and  Var(y) = 4.125
% 
disp("== For Yvec and Ydel = Yvec - Yavg with")
SYY             =  Yvar * (N-1)         ;
Cbias           =  Yvar/2               ; 
Smag            =  1                    ; 
Yvec_fun	=  Smag*sin(Xvec/3)     ;
txtOut          = ['   given Yvec constraints: SYY = ', num2str(SYY) ]     ;
disp(txtOut)

%% Yvec_fun= 0.0;
%  
disp("== Generation of an initial y adn Y = y - yavg using a random vector of 11 components.")
disp("== Adjusting Y1 and Yn by solving for Ym until 'RootTest > 0' is satisfied.")
RootTest= -1    ;
nIter   =  1    ;
while (RootTest < 0 )
  % rvec=[0.756 0.324 0.568 0.715 0.531 0.282 0.226 0.136 0.513 0.942 0.643] ;
  rvec		= rand(1,N)	;
  rvec(m)	= 0.5		; 
  SYYcond	= 0		;
  Ybias_init    = Cbias * ( 2 * rvec(1:N) - 1.0 )		; 
  Yvec_try_init = Beta0 + Beta1 * Xvec + Yvec_fun+ Ybias_init	;
  Ybias		= Ybias_init                                    ;
  Yvec_try	= Yvec_try_init                                 ;
  Yavg_try	= mean(Yvec_try)			; 
  Beta1_try	= var (Yvec_try) / var(Xvec)		;
  Beta0_try	= Yavg_try  - Beta1_try * Yavg_try	; 
  Yfit		= Beta0_try + Beta1_try * Xvec		; 
  Ydel_try	= Yvec_try  - Yavg			; 
  beta1_test	= var(Yvec_try)/var(Xvec)		; 
  a10 = Beta1*SXX + sum ((Xdel(N) - Xdel(2  :m-1) ).* Ydel_try(2  :m-1)) ...
                  + sum ((Xdel(N) - Xdel(m+1:N-1) ).* Ydel_try(m+1:N-1))              ;
  a1 = a10  / ( Xdel(1) - Xdel(N) )			;
  b1 = - ( Xdel(N) - Xdel(m) ) / ( Xdel(N) - Xdel(1) )  ;    
  an0 = Beta1*SXX + sum ((Xdel(1) - Xdel(2  :m-1) ).* Ydel_try(2  :m-1)) ...
                  + sum ((Xdel(1) - Xdel(m+1:N-1) ).* Ydel_try(m+1:N-1))              ; 
  an = an0  / ( Xdel(N) - Xdel(1) )			; 
  bn = - ( Xdel(m) - Xdel(1) ) / ( Xdel(N) - Xdel(1) )  ;    
  SYY_tmp = SYY - sum (Ydel_try(2  :m-1).*Ydel_try(2  :m-1)) ...
                - sum (Ydel_try(m+1:N-1).*Ydel_try(m+1:N-1)) ; 
  %% 
  B		=   ( a1*b1 + an*bn ) / (1 + b1^2 + bn^2)		;
  C		=   ( a1^2  + an^2  ) / (1 + b1^2 + bn^2)		;
  SYY_tmp	=        SYY_tmp      / (1 + b1^2 + bn^2)		;
  Ym0		= -B			;
  Ym1sq		=  SYY_tmp + B^2 - C	;
  Ym1		=  sqrt(Ym1sq)		;
  RootTest	=  sign(Ym1sq)          ;

  if (RootTest < 0)
      nIter     = nIter + 1 ;
      nRem      = mod(nIter,10);
      ordinal='th'  ;
      switch nRem
        case 1
          if (nIter ~= 11)
              ordinal='st'  ;
          end
        case 2
          if (nIter ~= 12)
              ordinal='nd'  ;
          end
        case 3 
          if (nIter ~= 13)
              ordinal='rd'  ;
          end
        otherwise
          ordinal='th'  ;
      end
      txtOut= ['   RootTest did not passed, i.e., RootTest < 0. Trying again ... The ',...
               num2str(nIter),ordinal,' iteration starts.' ];  
      disp(txtOut)
  end
end

if (nIter > 1 )
    txtOut= ['   RootTest successfully passed, i.e., RootTest > 0, with ', num2str(nIter),' iteration(s).' ]; 
    disp(txtOut)
end

%
if (rand(1,1) > 0.5 )
  Ym = Ym0 + Ym1 		; 
else
  Ym = Ym0 - Ym1		; 
end
    
Y1	= a1 + b1 * Ym	;
Yn	= an + bn * Ym  ;
Ydel    = Ydel_try	;
Ydel(1) = Y1		;
Ydel(N) = Yn		;
Ydel(m) = Ym		;
Yvec	= Ydel + Yavg	; 

%% 
xplot	= 0:20			;
yplot	= Beta0+Beta1*xplot     ;
figure(1)
plot(xplot ,yplot  ,'-b' ); hold on
plot(Xvec ,Yvec_try,"xr" ); hold on
plot(Xvec ,Yfit    ,"-.r"); hold on
plot(Xvec ,Yvec    ,"ob" ); hold on
axis([0 16 2 14])         ; hold off 

%%
avgYvec	 	= mean(Yvec)            ;
varYvec  	= var(Yvec)             ;
sumYYsq 	= (N-1) * varYvec       ;
txtOut1 = ['   Calculation results    : Y-mean = ', num2str(avgYvec), ...
          '; Y-var = ', num2str(varYvec),...
          '; sum-YYsq = ', num2str(sumYYsq) ] ;

txtOut2 = ['   Exact values (Anscombe): Y-mean = ', num2str(7.5), ...
          '; Y-var = ', num2str(4.125),...
          '; sum-YYsq = ', num2str(41.25) ] ; 
disp(txtOut1)
disp(txtOut2)
% 
%
