clear variables
close all

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
title('Przebiegi wyjscia dla skoku wartoœci sterowania w chwili 100s.')
xlabel('czas[t]')
ylabel('h2[cm]')
hold on
for i = 82:-10:22
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
        %nlin
        V1(1,t) = V1(1,t-1) + F1(t-1)+ FD(t-1) - a1*h1(1,t-1)^0.5;
        V2(1,t) = V2(1,t-1) + a1*h1(1,t-1)^0.5 - a2*h2(1,t-1)^0.5;
        h1(1,t) = (V1(1,t)/C1)^0.5;
        h2(1,t) = (V2(1,t)/C2)^0.5;
        %lin
        V1(2,t) = V1(2,t-1) + (F1(t-1) - F10) + (FD(t-1) - FD0) - a1/2*h10^-0.5 * (h1(2,t-1) - h10);
        V2(2,t) = V2(2,t-1) + a1/2*h10^-0.5 * (h1(2,t-1) - h10) - a1/2*h20^-0.5 * (h2(2,t-1) - h20);
        h1(2,t) = h10 + 1/2*(C1*V10)^-0.5 * (V1(2,t) - V10);
        h2(2,t) = h20 + 1/2*(C2*V20)^-0.5 * (V2(2,t) - V20);
    end
    plot(1:n,h2(1,1:end),'b')
    plot(1:n,h2(2,1:end),'m')
    legend('m.nlin', 'm.lin',Location='northwest')
end


figure
title('Przebiegi wyjscia dla skoku wartosci zaklocenia w chwili 100s.')
xlabel('czas[t]')
ylabel('h2[cm]')
hold on
for i = 18:-2:6
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
        %nlin
        V1(1,t) = V1(1,t-1) + F1(t-1)+ FD(t-1) - a1*h1(1,t-1)^0.5;
        V2(1,t) = V2(1,t-1) + a1*h1(1,t-1)^0.5 - a2*h2(1,t-1)^0.5;
        h1(1,t) = (V1(1,t)/C1)^0.5;
        h2(1,t) = (V2(1,t)/C2)^0.5;
        %lin
        V1(2,t) = V1(2,t-1) + (F1(t-1) - F10) + (FD(t-1) - FD0) - a1/2*h10^-0.5 * (h1(2,t-1) - h10);
        V2(2,t) = V2(2,t-1) + a1/2*h10^-0.5 * (h1(2,t-1) - h10) - a1/2*h20^-0.5 * (h2(2,t-1) - h20);
        h1(2,t) = h10 + 1/2*(C1*V10)^-0.5 * (V1(2,t) - V10);
        h2(2,t) = (h20 + 1/2*(C2*V20)^-0.5 * (V2(2,t) - V20)) ;
    end
    plot(1:n,h2(1,1:end),'b')
    plot(1:n,h2(2,1:end),'m')
    legend('m.nlin', 'm.lin', Location='northwest')
end