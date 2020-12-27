%% Generalised Smooth of Melanoma data
clear;
clc;

%% Load the Melanoma data
tempmat = load('melanoma.dat');
    
x  = tempmat(:,2);

y  = tempmat(:,3);

n=length(y);

scatter(x,y)