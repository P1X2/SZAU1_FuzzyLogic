function [a, c, hr20] = ZAD2_model_rozmyty(il, draw, c)


    n = 400;
    C1 = 0.95;
    C2 = 0.95;
    a1 = 17;
    a2 = 15;
    tau = 40;
    start = 100;

    ymax = 150;
    ymin = 5;



    dy = (ymax-ymin)/il;
    a = 3;
    if isempty(c)
        c = ymin+dy:dy:ymax-dy;
    end

    %pkt linearyzacji w kazdym zbiorze rozmytym
    hr20 = ones(1,il);
    hr20(1) = (c(1)+ymin)/2-1;
    hr20(il) = min((ymax+c(il-1))/2+1, ymax);

    if il > 2
        hr20(2:il-1) = (c(2:il-1)+c(1:il-2))./2;
    end
    
    FD0 = 12;
    hr10 = ((a2/a1)^2)*hr20;
    Vr20 = C2*hr20.^2;
    Vr10 = C1*hr10.^2;
    Fr0 = a1*hr10.^0.5-FD0;


    if draw
        figure
        hold on
        for r = 1:il
            if r == 1
                plot(ymin:0.1:ymax,1-1./(1+exp(-a*([ymin:0.1:ymax]-c(1)))))
            elseif r == il
                plot(ymin:0.1:ymax,1./(1+exp(-a*([ymin:0.1:ymax]-c(il-1)))))
            else
                plot(ymin:0.1:ymax,1./(1+exp(-a*([ymin:0.1:ymax]-c(r-1))))-1./(1+exp(-a*([ymin:0.1:ymax]-c(r)))))
            end
        end
        plot(hr20, ones(1,il), 'ko')
        xlabel('h2')
        ylabel('przynależności')
        title('Funkcje przynaleznosci modeli lokalnych w regulacji rozmytej względem wartosci wyjścia')
    end









    F10 = 52;
    FD0 = 12;
    h10 = 14.1730;
    h20 = 18.2044;
    V10 = C1*h10*h10;
    V20 = C2*h20*h20;

    if draw
        figure
        title('Przebiegi wyjścia dla skoku wartości sterowania w chwili 100s.')
        xlabel('czas[t]')
        ylabel('h2[cm]')
        hold on
    end

    for i = 82:-10:22
        F1in = F10 * ones(1,n);
        F1in(start:n) = i;
        F1 = F10 * ones(1,n);
        FD = FD0 * ones(1,n);
        h1 = h10 * ones(2+il+1,n);
        h2 = h20 * ones(2+il+1,n);
        V1 = C1*h10*h10 * ones(2+il+1,n);
        V2 = C2*h20*h20 * ones(2+il+1,n);
        w = ones(1,il);
        for t = tau+1 : n
            %nlin nierozmyty
            F1(t) = F1in(t-tau);
            V1(1,t) = V1(1,t-1) + F1(t-1)+ FD(t-1) - a1*h1(1,t-1)^0.5;
            V2(1,t) = V2(1,t-1) + a1*h1(1,t-1)^0.5 - a2*h2(1,t-1)^0.5;
            h1(1,t) = (V1(1,t)/C1)^0.5;
            h2(1,t) = (V2(1,t)/C2)^0.5;
            %lin nierozmyty
            V1(2,t) = V1(2,t-1) + (F1(t-1) - F10) + (FD(t-1) - FD0) - a1/2*h10^-0.5 * (h1(2,t-1) - h10);
            V2(2,t) = V2(2,t-1) + a1/2*h10^-0.5 * (h1(2,t-1) - h10) - a1/2*h20^-0.5 * (h2(2,t-1) - h20);
            h1(2,t) = h10 + 1/2*(C1*V10)^-0.5 * (V1(2,t) - V10);
            h2(2,t) = h20 + 1/2*(C2*V20)^-0.5 * (V2(2,t) - V20);

            %lin rozmyty
            for r = 1:il
                V1(2+r,t) = V1(2+il+1,t-1) + (F1(t-1) - Fr0(r)) + (FD(t-1) - FD0) - a1/2*hr10(r)^-0.5 * (h1(2+il+1,t-1) - hr10(r));
                V2(2+r,t) = V2(2+il+1,t-1) + a1/2*hr10(r)^-0.5 * (h1(2+il+1,t-1) - hr10(r)) - a1/2*hr20(r)^-0.5 * (h2(2+il+1,t-1) - hr20(r));
                h1(2+r,t) = hr10(r) + 1/2*(C1*Vr10(r))^-0.5 * (V1(2+r,t) - Vr10(r));
                h2(2+r,t) = hr20(r) + 1/2*(C2*Vr20(r))^-0.5 * (V2(2+r,t) - Vr20(r));


                if r == 1
                    w(r) = 1-1/(1+exp(-a*(h2(2+il+1,t-1)-c(1))));
                elseif r == il
                    w(r) = 1/(1+exp(-a*(h2(2+il+1,t-1)-c(il-1))));
                else
                    w(r) = 1/(1+exp(-a*(h2(2+il+1,t-1)-c(r-1)))) - 1/(1+exp(-a*(h2(2+il+1,t-1)-c(r))));
                end


            end
            h2(2+il+1,t) = w*h2(3:2+il, t)/sum(w);
            h1(2+il+1,t) = w*h1(3:2+il, t)/sum(w);
            V1(2+il+1,t) = w*V1(3:2+il, t)/sum(w);
            V2(2+il+1,t) = w*V2(3:2+il, t)/sum(w);

        end
        
        if draw
            plot(start:n,h2(1,start:end),'b')
            plot(start:n,h2(2,start:end),'m')
            plot(start:n,h2(2+il+1,start:end),'g')
            legend('m.nlin', 'm.lin', 'm.rozmyty', Location='northwest')
        end
% 
    end
%     
end