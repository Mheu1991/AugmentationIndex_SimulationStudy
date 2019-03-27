% Plot relation between left ventricular stroke work (W stroke ) and wave
% reflection indices (i.e. augmentation index (AIx) and backward pressure amplitude (p_bw)
% Maarten Heusinkveld 16-03-2019
clear all;
close all;
clc;

load('Result2.mat')
close all

ok = ~isnan(Result.aix) & ~isnan(Result.wstroke);

%%
vsI = scatteredInterpolant(Result.aix(ok), Result.wstroke(ok), Result.vs(ok), 'natural', 'none');
kI  = scatteredInterpolant(Result.aix(ok), Result.wstroke(ok), Result.k(ok), 'natural', 'none');
[X,Y] = meshgrid(16:0.1:62,0.90:0.001:1.22);

figure
surf(X,Y,vsI(X,Y), kI(X,Y), 'LineStyle', 'none')
xlabel('AIX')
ylabel('W_{stroke}')
zlabel('V_{s}')
ca = colorbar;
ylabel(ca, 'k')
hold on
plot3(Result.aix(ok), Result.wstroke(ok), Result.vs(ok), 'ok')

%%
vs = Result.vs(ok);
nVs = 9;
vsScaledI = (Result.vs(ok) - min(Result.vs(ok)))./ ...
  (max(Result.vs(ok))-min(Result.vs(ok)));
vsScaledI = vsScaledI * (nVs-1)/nVs;

kScaledI = (Result.k(ok) - min(Result.k(ok)))./ ...
  (max(Result.k(ok))-min(Result.k(ok)));
%aix = Result.aix(ok)
pb  = Result.pb(ok);
wstroke = Result.wstroke(ok);


rgbM = hsv2rgb([vsScaledI,kScaledI,ones(size(vsScaledI))]);
%rgbM = hsv2rgb([vsScaledI,kScaledI,kScaledI]);

figure
markersC = {'o','s','d','^','v','<','>','p','h'}
vstemp = 2
legendHV = [1];
for i = 1:length(vsScaledI)
  %h = plot(aix(i), wstroke(i), 'MarkerFaceColor', rgbM(i,:), 'Marker', markersC{vs(i)-1}, 'LineStyle', 'none', 'MarkerEdgeColor', 'k');
  h = plot(pb(i), wstroke(i), 'MarkerFaceColor', rgbM(i,:), 'Marker', markersC{vs(i)-1}, 'LineStyle', 'none', 'MarkerEdgeColor', 'k');

  hold on
  if vs(i) > vstemp
    legendHV = [legendHV;oldH];
    vstemp = vs(i);
  else
    oldH = h;
  end
end
set(gca,'Color',0.9*[1,1,1]);
legendHV = legendHV(2:end);
legendHV = [legendHV;h]; %Don't forget last one!
lh = legend(legendHV, arrayfun(@(x) num2str(x),2:10, 'UniformOutput', false));
set(lh, 'Color', 'w');

figure
[X,Y] = meshgrid(0:0.01:1, 0:0.01:1);

CC = arrayfun(@(x,y) hsv2rgb(x,y,1), X,Y, 'UniformOutput', false);

C = zeros([size(X), 3]);
C(:,:,1) = cellfun(@(x) x(1),CC);
C(:,:,2) = cellfun(@(x) x(2),CC);
C(:,:,3) = cellfun(@(x) x(3),CC);

surf(X,Y,zeros(size(X)),C, 'LineStyle','none')
view(0,90)


