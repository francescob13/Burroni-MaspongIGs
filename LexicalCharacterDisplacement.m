function LexicalCharacterDisplacement(~)
%% Simulate Lexical Character Displacement 
%
% This script was used to generate the simulation of degemination and
% degemination inhibition presented in Burroni and Maspong (2020). The
% exemplar model is inspired by Blevins and Wedel 2009 and Wedel 2004.
%
% Coded by Francesco Burroni 9/9/2020

% Clear workspace, close all figures
clearvars global; close all;

% Set seed
rng(60)

%% Model parameters

muDurIG= 125;
muDurSingleton= 75;
sigma= 15; % standard deviation
nExemplars= 100; % n exemplars in each population
nEntranchment= 3 ; % n exemplars used to calculate production target by taking
% their mean
nIterations= 10000;
noise= 25; % noise added to production is uniformly distributed in
% (-noise,+ noise)

%% Initialize Populations

% Singleton Population N(muDurSingleton, sigma)
NSingleton= randn(1,nExemplars)*sigma+muDurSingleton;

% IGs Population N(muDurIG, sigma)
NGeminate= randn(1,nExemplars)*sigma+muDurIG;

%% Plot initial position of singleton and IG in closure duration space

bw= 8; % bandwidth for ksd
C= colororder; % get Matlab default colors

% Figure
f = figure('WindowState','fullscreen');
subplot(6,6,[7 8 9 13 14 15 19 20 21 25 26 27]);  hold on;

% Plot filled ksdensities
filledKsdensity(NSingleton,bw,C(1,:));
filledKsdensity(NGeminate,bw,C(2,:));

% Change appearance of subplot
xL= xline(99.7,"LineStyle","--","label","category boundary",...
    "HandleVisibility","off","FontName","Fira Sans","FontSize",30);
ylim([0 .05])
xlim([10 190])
yticks([])
xlabel("Closure Duration [ms]")
legend("Singleton","Initial Geminate");
set(gca,"FontSize",20,"Fontname","Fira Sans")
t= title("a)","FontSize",20,"Fontname","Fira Sans");
t.Position(1)= 20;
t.Position(2)= 0.048;

%% Simulate degemination

% Update distributions allowing external disambiguation to block
% variant trading
d1= updateDistributionExternal(NSingleton,nEntranchment,nIterations,noise);
d2= updateDistributionExternal(NGeminate,nEntranchment,nIterations,noise);

% Subplot
subplot(6,6,[4 5 6 10 11 12 16 17 18]);  hold on;

% Plot filled ksdensities
filledKsdensity(d1,bw,C(1,:));
filledKsdensity(d2,bw,C(2,:));

% Change appearance of subplot
ylim([0 .04])
xlim([10 190])
yticks([])
set(gca,"FontSize",20,"Fontname","Fira Sans")
t= title("b)","FontSize",20,"Fontname","Fira Sans");
t.Position(1)= 20;
t.Position(2)= 0.038;

%% Simulate degemination inhibition

% Update distributions allowing variant trading on the basis of distance
% from mean of each distribution
[d1, d2]= updateDistribution(NSingleton,NGeminate,nEntranchment,nIterations,...
    noise);

% Change appearance of subplot
subplot(6,6,[22 23 24 28 29 30 34 35 36]);  hold on;
filledKsdensity(d1,bw,C(1,:));
filledKsdensity(d2,bw,C(2,:));
ylim([0 .04])
xlim([10 190])
yticks([])
xlabel("Closure Duration [ms]")
set(gca,"FontSize",20,"Fontname","Fira Sans")
t= title("c)","FontSize",20,"Fontname","Fira Sans");
t.Position(1)= 20;
t.Position(2)= 0.038;
hold off;

saveas(f,"LexicalCharacterDisplacement.png")
close all

fprintf("The Lexical Character Displacement \n model has been printed\n")
end

%% Local Functions

function h= filledKsdensity(v,bandwidth,C)

% Plot ksdensity
[f,xi]= ksdensity(v,"Bandwidth",bandwidth);

% Fill ksdensity
h= fill([xi,flip(xi)],[zeros(1,length(xi)),flip(f)],C(1,:),...
    "Linestyle","none","FaceAlpha",.4);

% Change appearance of contour line
l= plot(xi,f,"Color",C(1,:),"LineWidth",2);
l.HandleVisibility = "off";
end

function d= updateDistributionExternal(d,nExemplars,nIterations,noise)

% Update a distribution one exemplar at a time with no variant trading
% due to external disambiguation (Blevins and Wedel 2009)

for k= 1 : nIterations
    
    % Get ixs of n random exemplars
    exemplarsIx= randi(length(d),1,nExemplars);
    
    % Store values of exemplars in a vector
    exemplarsVector= d(exemplarsIx);
    
    % Add noise uniformly distributed in (-noise, noise)
    noiseVal= rand*noise*2-noise;
    
    % Compute value of new exemplar using mean of exemplar vector plus
    % noise, Blevins&Wedel (2009) explain why we want to use the mean,
    % similar to Pierrehumbert's (2001) entrenchment
    newexemplar= mean(exemplarsVector) + noiseVal;
    
    % Get ix of an exemplar to be replaced by new one
    replaceExemplarIx= randi(length(d));
    
    % Replace old exemplar with new one
    d(replaceExemplarIx)= newexemplar;
    
end

end

function [d1,d2] = updateDistribution(d1,d2,nExemplars,nIterations,noise)

% Update two distributions one exemplar at a time, allow variant trading
% based on minimum distance between new exemplar value and mean of d1 and d2
% Wedel 2004

for k= 1 : nIterations
    
    if rem(k,2) == 0 % At even iterations draw from d1...
        
        % Get ixs of n random exemplars
        exemplarsIx= randi(length(d1),1,nExemplars);
        
        % Store values of exemplars in a vector
        exemplarsVector= d1(exemplarsIx);
        
        % Add noise uniformly distributed in (-noise, noise)
        noiseVal= rand*noise*2-noise;
        
        % Compute value of new exemplar
        newexemplar= mean(exemplarsVector) + noiseVal;
        
    else % ...at odd iterations draw from d2...
        
        exemplarsIx= randi(length(d2),1,nExemplars);
        
        exemplarsVector= d2(exemplarsIx);
        
        noiseVal= rand*noise*2-noise;
        
        newexemplar= mean(exemplarsVector) + noiseVal;
    end
    
    % Calculate distance between new exemplar and mean of each distribution
    distanced1= abs(mean(d1) - newexemplar);
    distanced2= abs(mean(d2) - newexemplar);
    
    if distanced1 < distanced2 % if closer to d1
        
        % Get ix of exemplar to be replaced by new one in d1
        newexemplarIx= randi(length(d1));
        
        % Replace old exemplar with new one in d1
        d1(newexemplarIx)= newexemplar;
        
    else  % if closer to d2
        
        % Get ix exemplar to be replaced by new one in d2
        newexemplarIx= randi(length(d2));
        
        % Replace old exemplar with new one in d2
        d2(newexemplarIx)= newexemplar;
    end
    
end

end
