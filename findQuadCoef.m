%
% AUTHOR: Albert S. Kim
% Email: albertsk@hawaii.edu
% Affiliation: Civil and Environmental Engineering, Unviersity of Hawaii at Manoa
% DATE: 2021-07-24T22:30:48 HST 
% PURPOSE: As a part of paper entitled "Inverse Sampling of Degenerate Datasets from a Linear Regression Line"
%          submitted to Journal of American Statistical Association 
% FILENAME: findQuadCoef.m
% 
clear all       ; close all     ; format long 
N	= 11  ;
Beta0	= 3.0 ;
Beta1	= 0.5 ;

%%
m	= (N+1)/2       ;

%%
Xavg    = 9             ;
Xvar    = 11            ; 
Xstd    = sqrt(Xvar)    ; 

%%
Yavg    = 7.5009                ;
Yvar    = 4.125                 ; 
Yvar    = 4.12762909090909      ;
Ystd    = sqrt(Yvar)            ;

%%
sumx_square     = 0 ;
sumx_cube       = 0 ;
sumx_fourth     = 0 ;

%%
X0 	=  3            ; 
Xvec	=  1:N          ;
Xvec 	=  Xvec + X0    ; 

%%
smin	=  4                    ;
smax 	=  14                   ;
Ndivs  	=  1000000              ;
Ns     	=  Ndivs + 1            ;
ds  	=  (smax - smin)/ Ndivs ; 

%%
Xstrvec (1:Ns) 		=   0           ;
Yvarvec (1:Ns) 		=   0           ;
Yvarvec_const (1:Ns)    =   Yvar        ;
Alphavec (1:Ns)         =   0           ;
Q0vec (1:Ns) 		=   0           ;

%%
for j = 1:Ns
    sumx_square 	=  0                    ;
    sumx_cube   	=  0                    ;
    sumx_fourth 	=  0                    ;
    xstar               =  smin + ds * (j-1)    ; 
    for i = 1:N
        sumx_square 	=  sumx_square + (Xvec(i)-xstar)^2                      ; 
        sumx_cube   	=  sumx_cube   + (Xvec(i)-xstar)^2 * (Xvec(i)-Xavg )    ; 
        sumx_fourth 	=  sumx_fourth + (Xvec(i)-xstar)^4                      ; 
    end
    Alpha 	=  Beta1 *(N-1) * Xvar / sumx_cube                      ;
    Q0    	=  Yavg - Alpha/N*sumx_square                           ;
    Yvar_test 	=  Alpha^2/(N-1)*sumx_fourth - N/(N-1)*(Q0-Yavg)^2      ;
    
    Xstrvec (j)	=   xstar               ;
    Alphavec(j)	=   Alpha               ;
    Q0vec   (j)	=   Q0                  ;
    Yvarvec (j)	=   Yvar_test           ;
end

%%
figure(1); plot(Xstrvec,Yvarvec-Yvarvec_const, 'LineWidth',3)
title("Yvarvec-Yvarvec const")
xlim([ 4 14 ])
ylim([ -2 2 ])
line([9,9], ylim, 'Color', 'k', 'LineWidth', 2);
line(xlim, [0,0], 'Color', 'k', 'LineWidth', 2);
grid on 

%%
figure(2); plot(Xstrvec,Alphavec, 'LineWidth',3)
title("Alphavec")
xlim([  8.95 9.05 ])
ylim([ -250  250  ])
line([9,9], ylim, 'Color', 'k', 'LineWidth', 2); 
line(xlim, [0,0], 'Color', 'k', 'LineWidth', 2); 
grid on 

%% printing outputs to a file. 
x 	=  0:.1:1                       ;
A 	=  [x; exp(x)]                  ;
xydata	=  [Xstrvec; Yvarvec ]          ;
outFile =  'outputDiffVarY.txt'               ;
fileID 	=  fopen(outFile,'w')    ;
fprintf(fileID,'%6.2f %12.8f\n',xydata) ;
fclose(fileID)                          ;

for i = 1: Ns-1
    if ((Yvarvec(i) - Yvar) * (Yvarvec(i+1) - Yvar) <0 )
        xval = Xstrvec  (i) ; 
        aval = Alphavec (i) ; 
        qval = Q0vec    (i) ;
        txtOut=['   [Result] Arrary index =',  num2str(i), ...
                '; Approx. Solution = ', num2str(xval),...
                '; Approx. Alpha = ', num2str(aval), ...
                '; Q value = ', num2str(qval) ]; 
        disp(txtOut)
    end
end
disp(['   "',outFile,'" is generated. '])
%% Final Answers of two cases

% ans =     7
% ans =   0.125000000000000
% ans =   5.750900000000000

% ans =    11
% ans =  -0.125000000000000
% ans =   9.250900000000000
