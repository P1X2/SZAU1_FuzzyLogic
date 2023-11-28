clear;

n = 400;
C1 = 0.95;
C2 = 0.95;
a1 = 17;
a2 = 15;
tau = 40;
start = 100;

% punkt pracy
F10 = 52;
FD0 = 12;
h10 = 14.1730;
h20 = 18.2044;

V10 = C1*h10^2;
V20 = C2*h20^2;

figure
title('Przebiegi wyjscia dla skoku wartosci sterowania w chwili 100.')
xlabel('czas[t]')
ylabel('h2[cm]')
hold on
for i = 82:-10:22
    F1in = F10 * ones(1,n);
    F1in(start:n) = i;
    F1 = F10 * ones(1,n);
    FD = FD0 * ones(1,n);
    h1 = h10 * ones(1,n);
    h2 = h20 * ones(1,n);
    V1 = C1*h10*h10 * ones(1,n);
    V2 = C2*h20*h20 * ones(1,n);
    for t = tau+1 : n
        F1(t) = F1in(t-tau);
        V1(1,t) = V1(1,t-1) + F1(t-1)+ FD(t-1) - a1*h1(1,t-1)^0.5;
        V2(1,t) = V2(1,t-1) + a1*h1(1,t-1)^0.5 - a2*h2(1,t-1)^0.5;
        h1(1,t) = (V1(1,t)/C1)^0.5;
        h2(1,t) = (V2(1,t)/C2)^0.5;
    end
    plot(start:n,h2(1,start:end))
    legend('82', '72', '62', '52', '42','32', '22')
end


figure
title('Przebiegi wyjscia dla skoku wartosci zaklocenia w chwili 100.')
xlabel('czas[t]')
ylabel('h2[cm]')
hold on
for i = 18:-2:6
    F1in = F10 * ones(1,n);
    F1 = F10 * ones(1,n);
    FD = FD0 * ones(1,n);
    FD(start:n) = i;
    h1 = h10 * ones(1,n);
    h2 = h20 * ones(1,n);
    V1 = C1*h10*h10 * ones(1,n);
    V2 = C2*h20*h20 * ones(1,n);
    for t = tau+1 : n
        F1(t) = F1in(t-tau);
        V1(1,t) = V1(1,t-1) + F1(t-1)+ FD(t-1) - a1*h1(1,t-1)^0.5;
        V2(1,t) = V1(1,t-1) + a1*h1(1,t-1)^0.5 - a2*h2(1,t-1)^0.5;
        h1(1,t) = (V1(1,t)/C1)^0.5;
        h2(1,t) = (V2(1,t)/C2)^0.5;
    end
    plot(start:n,h2(1,start:end))
end
legend('18','16', '14', '12', '10', '8','6')