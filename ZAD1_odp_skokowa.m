function [s,z] = ZAD1_odp_skokowa(Y0, draw)

    n = 1000;
    tau = 40;
    D0 = 12;
    U0 = 52;
    start = 100;
    C1 = 0.95;
    C2 = 0.95;
    a1 = 17;
    a2 = 15;
    step = 1;


    h20 = Y0;
    V20 = C2 * h20^2;
    h10 = (a2/a1)^2 * h20;
    V10 = C1*h10^2;

    U = U0*ones(1,n);
    D = D0*ones(1,n);
    U(1,start:n) = U0 + step;
    Y = h20*ones(1,n);
    V1 = V10;
    V2 = V20;
    h1 = h10;
    for i = start+1:n
        V1 = V1 + (U(i-1-tau) - U0) + (D(i-1) - D0) - a1/2*h10^-0.5 * (h1 - h10);
        V2 = V2 + a1/2*h10^-0.5 * (h1 - h10) - a2/2*h20^-0.5 * (Y(i-1) - h20);
        h1 = h10 + 1/2*(C1*V10)^-0.5 * (V1 - V10);
        Y(i) = h20 + 1/2*(C2*V20)^-0.5 * (V2 - V20);
    end
    s = (Y(start+1:n)-Y0)/step;

    U = U0*ones(1,n);
    D = D0*ones(1,n);
    D(1,start:n) = D0 + step;
    Y = h20*ones(1,n);
    V1 = V10;
    V2 = V20;
    h1 = h10;
    for i = start+1:n
        V1 = V1 + (U(i-1-tau) - U0) + (D(i-1) - D0) - a1/2*h10^-0.5 * (h1 - h10);
        V2 = V2 + a1/2*h10^-0.5 * (h1 - h10) - a2/2*h20^-0.5 * (Y(i-1) - h20);
        h1 = h10 + 1/2*(C1*V10)^-0.5 * (V1 - V10);
        Y(i) = h20 + 1/2*(C2*V20)^-0.5 * (V2 - V20);
    end
    z = (Y(start+1:n)-h20)/step;
    if draw
        figure
        subplot(2,1,1)
        plot(s)
        ylabel('s')
        xlabel('k')
        subplot(2,1,2)
        plot(z)
        ylabel('z')
        xlabel('k')
    end
end