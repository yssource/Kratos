function [] = MainCode()

for i = [1:3]
    figure(i)
    hold off;
end

global M;
global ShapeN;
global SpacingR;

M = 0.98;
ShapeN = 10;
SpacingR = 12;

ThisExists = true;
for i = [0:100]
%     ThisExists = MakeThisFile([num2str(i), 'drained_triaxial.csv']);
    ThisExists = MakeThisFile([num2str(i),'undrained_triaxial.csv']);
    if (ThisExists == false )
        break
    end
    pause(0.01)
    %      MakeThisFile([num2str(i),'oedometer.csv']);
    %      MakeThisFile([num2str(i),'isotropic.csv']);
    
end


figure(2)
set(gca, 'FontSize', 12)
xlabel('$\| \epsilon^d \|$', 'interpreter', 'latex')
ylabel('$q$ (kPa)', 'interpreter', 'latex')
yy = ylim();

ylim([0, yy(2)]);

figure(1)
set(gca, 'FontSize', 12)
xlabel('$p^{\prime}$ (kPa)', 'interpreter', 'latex')
ylabel('$q$ (kPa)', 'interpreter', 'latex')
xlim([0, 80])
ylim([0, 60])

figure(1)
print('EffectiveStressPath-Casm', '-dpng')
figure(2)
print('DevDef-DevStress-Casm', '-dpng')

function [ThisExists] = MakeThisFile(XFILE)
ThisExists = true;

if ( isfile(XFILE) == false)
    ThisExists = false;
    return
end
rawData = csvread(XFILE);

global ShapeN;
global SpacingR;

ShapeN = rawData(1,1);
SpacingR = rawData(1,2);
rawData = rawData(2:end,:);

Stress = rawData(:,2:7);
Strain = rawData(:,8:13);

ps = rawData(:,14);
ps(2)

[p, J] = ComputeStressInvariants( Stress);
[eVol, eDev] = ComputeStrainInvariants( Strain);
%
% Slope(1) = 0;
% for k = 2:length(eVol)
%     Slope(k) = log(p(k)/p(k-1));
%     Slope(k) = Slope(k)/(eVol(k)-eVol(k-1));
% end
% figure(121)
% semilogy(p, Slope)
% hold on


figure(1)
PlotYieldSurface( ps(1),  'k-.');
% PlotYieldSurface( ps(end), 'g-.');

% Triaxial plane
plot(p, J, 'linewidth', 1.4);
xlabel('p')
ylabel('q')
axis equal

figure(2)
plot( eDev-eDev(1), J, 'linewidth', 1.4);
xlabel('eDev')
ylabel('sDev')
hold on



% figure(3)
% semilogx( p, -eVol);
% xlabel('p')
% ylabel('eVol')
% hold on

max(J)/2
J(end)/2

function PlotYieldSurface( p0, SPEC)

global M;
global ShapeN;
global SpacingR;
n = ShapeN;
r = SpacingR;

p0 = -p0;


pp = linspace(0, p0, 1000);



plot([0, 1.5*p0], M*[0, 1.5*p0], 'r-.')
hold on

qq = M*pp .* ( -log(pp/p0)/log(r)).^(1/n);
% plot(pp, qq, SPEC, 'linewidth', 1.0)


axis equal




function [p, J] = ComputeStressInvariants( Stress)

p = sum(Stress(:,1:3)')/3;

for i = 1:size(Stress,1)
    J(i) = 0;
    for e = 1:3
        J(i) = J(i) + (Stress(i,e)-p(i))^2;
    end
    for e = 4:6
        J(i) = J(i) + 2*(Stress(i,e))^2;
    end
    J(i) = sqrt(0.5*J(i))*sqrt(3);
end
p = - p;


function [eVol, eDev] = ComputeStrainInvariants( Strain)

eVol = sum(Strain(:,1:3)')/3;

for i = 1:size(Strain,1)
    eDev(i) = 0;
    for e = 1:3
        eDev(i) = eDev(i) + (Strain(i,e)-eVol(i))^2;
    end
    for e = 4:6
        eDev(i) = eDev(i) + 2*(Strain(i,e)/2.0)^2;
    end
    eDev(i) = sqrt(0.5*eDev(i))*sqrt(3);
end
eVol = - eVol*3;
