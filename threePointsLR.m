% 
% AUTHOR: Albert S. Kim
% Email: albertsk@hawaii.edu
% Affiliation: Civil and Environmental Engineering, Unviersity of Hawaii at Manoa
% DATE: 2021-07-24T22:30:48 HST 
% PURPOSE: As a part of paper entitled "Inverse Sampling of Degenerate Datasets from a Linear Regression Line"
%          submitted to Journal of American Statistical Association 
% FILENAME: threePointsLR.m
% 
clear all; close all 
% 
xbar    =  9               ;
xvar    = 11               ;
dx1     = - sqrt(xvar)     ; 
dx2     = 0                ; 
dx3     = - dx1            ; 

x1 	= xbar +  dx1   ;
x2 	= xbar +  dx2   ;
x3 	= xbar +  dx3   ;

beta1   = 0.5                   ;
yvar    = 4.125                 ;
B1      = beta1 * xvar / dx3 ;
X       = [x1 x2 x3]            ;
ysign   = -1                    ; 
figure(1)
for i = -1:2:1
    ysign = i                                           ; 
    dy2 = ysign * 2/sqrt(3) * sqrt(yvar - B1^2)         ;
    dy1 = - dy2/2.0 - B1                                ;
    dy3 = - dy2/2.0 + B1                                ;
    ybar = 7.5                  ;
    y1 = dy1 + ybar             ;
    y2 = dy2 + ybar             ;
    y3 = dy3 + ybar             ;
    Y = [ y1 y2 y3]             ;
    if ( i < 0 )
        plot(X,Y,"sb",'MarkerSize',12,'LineWidth',2)    ; hold on 
        plot(X,Y,"-.b",'LineWidth',1)                   ; hold on 
    else
        plot(X,Y,"or",'MarkerSize',12,'LineWidth',2)    ; hold on
        plot(X,Y,"-.r",'LineWidth',1)                   ; hold on
    end 
end

xfit    = 0:20                          ;
yfit    = 3 + 0.5 * xfit                ;
plot(xfit,yfit,"-k",'LineWidth',2)      ;

axis ([0 20 0 15])      ;
xlabel("x")             ;
ylabel("y")             ;
xticks([0:2:20])        ;
yticks([1:2:15])        ; 
grid on;





