function P = Pressure_cal_truncation(gap)

global delta_L threshold  randomFAI Parameter Penalty

gamma = 1.5;
M = 10;
L = 0.001;
NL = 512 *1;
delta_L = L/(NL+1);

x = 1e-8:L/NL:1e-8+L;   x = x';
y = 1e-8:L/NL:1e-8+L;
z = zeros(length(x),length(y));

for m = 1:M
    for n = Parameter.n_min:Parameter.n_max
        z = z + gamma^((Parameter.Dim-3)*n).*(cos(randomFAI(m,n+1-Parameter.n_min)) - ...
            cos((2*pi*gamma^n.*(x.^2+y.^2).^0.5./L).*cos(atan(y./x)-pi*m/M) + randomFAI(m,n+1-Parameter.n_min)) );
    end
end

z1 = zeros(length(x), 1);
m = 1;
for n = Parameter.n_min:Parameter.n_max
    z1 = z1 + gamma^((Parameter.Dim-3)*n) .* (cos(randomFAI(m,n+1-Parameter.n_min)) - cos((2*pi*gamma^n.*(x)/L - randomFAI(m,n+1-Parameter.n_min)))  );
end

z1 = z1 * L*(Parameter.G/L)^(Parameter.Dim-2)*(log(gamma))^0.5;
z1 = z1 - mean(z1);
[sigma, ~, k, ~] = Surface_parameter(x, z1);


z = z  *  L*(Parameter.G/L)^(Parameter.Dim-2)*(log(gamma)/M)^0.5;
z = z - mean(z(:));


Rt = (max(z(:))-min(z(:)));
ht = (Rt-gap)/Rt;

if ht < 0
    P = 0; return;
elseif ht > 0.3
%     Continue to run, or terminate, or change to  penalty method
%     P = 0.001*Penalty*ht*Rt; return;
end

threshold = max(z(:))  -  ht*Rt;

z_fig = z;
z_fig(z_fig < threshold) = NaN;


global E_roughness H kZ 
E1 = Parameter.E_rou;   E2 = Parameter.E_rou;
miu1 = Parameter.miu;   miu2 = Parameter.miu;
E_roughness = 1/((1-miu1^2)/E1 + (1-miu2^2)/E2);
kZ = 0.454+0.43*miu1;

H = 4*E_roughness/(3*pi*kZ*Parameter.fai)*sqrt(sigma/k);


regionStats = analyzeSurface(z_fig, threshold);
rowsToDelete = cellfun(@(x) x == 0, regionStats(:,5));  
regionStats(rowsToDelete, :) = [];  

asperity_w = zeros(length(regionStats(:,5)),1);  
asperity_R = zeros(length(regionStats(:,5)),1);
asperity_f = zeros(length(regionStats(:,5)),1);
asperity_k = zeros(length(regionStats(:,5)),1);
asperity_a = zeros(length(regionStats(:,5)),2);

for i = 1:length(regionStats(:,1))
    asperity_w(i,1) = regionStats{i,3} - threshold;
    asperity_R(i,1) = asperity_R_cal(asperity_w(i,1), regionStats{i,2} );
    asperity_f(i,1) = asperity_f_cal(asperity_R(i,1),  asperity_w(i,1));
%     asperity_k(i,1) = asperity_k_cal(asperity_R(i,1),  asperity_w(i,1));
%     asperity_a(i,:) = asperity_a_cal(asperity_R(i,1),  asperity_w(i,1));
%     r_real(i) = sqrt(asperity_a(i,1)/pi);
%     r_t(i) = sqrt(regionStats{i,2}/pi);
end




P = sum(asperity_f)/L^2;




end