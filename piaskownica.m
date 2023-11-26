clear variables
close all

n = 400;
C1 = 0.95;
C2 = 0.95;
a1 = 16;
a2 = 16;
tau = 50;
start = 60;

F10 = 54;
FD0 = 10;
h10 = 17;
h20 = 16;
V10 = C1*h10*h10;
V20 = C2*h20*h20;

figure
title('Przebiegi wyjœcia dla skoku wartoœci sterowania w chwili 60s.')
xlabel('czas[t]')
ylabel('h2[cm]')
hold on
for i = 84:-10:24
    F1in = F10 * ones(1,n);
    F1in(start:n) = i;
    F1 = F10 * ones(1,n);
    FD = FD0 * ones(1,n);
    h1 = h10 * ones(2,n);
    h2 = h20 * ones(2,n);
    V1 = C1*h10*h10 * ones(2,n);
    V2 = C2*h20*h20 * ones(2,n);
    for t = tau+1 : n
        F1(t) = F1in(t-tau);
        V1(1,t) = V1(1,t-1) + F1(t-1)+ FD(t-1) - a1*h1(1,t-1)^0.5;
        V2(1,t) = V2(1,t-1) + a1*h1(1,t-1)^0.5 - a2*h2(1,t-1)^0.5;
        h1(1,t) = (V1(1,t)/C1)^0.5;
        h2(1,t) = (V2(1,t)/C2)^0.5;
        
        V1(2,t) = V1(2,t-1) + (F1(t-1) - F10) + (FD(t-1) - FD0) - a1/2*h10^-0.5 * (h1(2,t-1) - h10);
        V2(2,t) = V2(2,t-1) + a1/2*h10^-0.5 * (h1(2,t-1) - h10) - a1/2*h20^-0.5 * (h2(2,t-1) - h20);
        h1(2,t) = h10 + 1/2*(C1*V10)^-0.5 * (V1(2,t) - V10);
        h2(2,t) = h20 + 1/2*(C2*V20)^-0.5 * (V2(2,t) - V20);
    end
    plot(start:n,h2(1,start:end),'b')
    plot(start:n,h2(2,start:end),'m')
end

figure
title('Przebiegi wyjœcia dla skoku wartoœci zak³ocenia w chwili 60s.')
xlabel('czas[t]')
ylabel('h2[cm]')
hold on
for i = 15:-2.5:5
    F1in = F10 * ones(1,n);
    F1 = F10 * ones(1,n);
    FD = FD0 * ones(1,n);
    FD(start:n) = i;
    h1 = h10 * ones(2,n);
    h2 = h20 * ones(2,n);
    V1 = C1*h10*h10 * ones(2,n);
    V2 = C2*h20*h20 * ones(2,n);
    for t = tau+1 : n
        F1(t) = F1in(t-tau);
        V1(1,t) = V1(1,t-1) + F1(t-1)+ FD(t-1) - a1*h1(1,t-1)^0.5;
        V2(1,t) = V1(1,t-1) + a1*h1(1,t-1)^0.5 - a2*h2(1,t-1)^0.5;
        h1(1,t) = (V1(1,t)/C1)^0.5;
        h2(1,t) = (V2(1,t)/C2)^0.5;
        
        V1(2,t) = V1(2,t-1) + (F1(t-1) - F10) + (FD(t-1) - FD0) - a1/2*h10^-0.5 * (h1(2,t-1) - h10);
        V2(2,t) = V2(2,t-1) + a1/2*h10^-0.5 * (h1(2,t-1) - h10) - a1/2*h20^-0.5 * (h2(2,t-1) - h20);
        h1(2,t) = h10 + 1/2*(C1*V10)^-0.5 * (V1(2,t) - V10);
        h2(2,t) = h20 + 1/2*(C2*V20)^-0.5 * (V2(2,t) - V20);
    end
    plot(start:n,h2(1,start:end),'b')
    plot(start:n,h2(2,start:end),'m')
end