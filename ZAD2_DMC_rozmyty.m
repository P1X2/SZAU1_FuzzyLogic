function[E,e] = ZAD2_DMC_rozmyty(D, N, Nu, DZ, lambda, a, c, draw, pomiar_z)

    il = 15;
    lambda = lambda*ones(1,il);
    ku = zeros(il,D-1);
    kz = zeros(il,DZ);
    ke = zeros(1,il);


    %%
    ymin = 5;
    ymax = 150;
    if isempty(a)
        [a, c, hr20] = ZAD2_model_rozmyty(il, false, []);
    else
        hr20(1) = (c(1)+ymin)/2-1;
        hr20(il) = min((ymax+c(il-1))/2+1, ymax);
        if il > 2
            hr20(2:il-1) = (c(2:il-1)+c(1:il-2))./2;
        end
    end

    
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
        title('Funkcje przynależnoœci regulatorów lokalnych w regulacji rozmytej względem wartości wyjścia')
    end

    %%
    
    for r = 1:il
        [s_ob, s_zak] = ZAD1_odp_skokowa(hr20(r), false);

    % macierz M obiektu
    M=zeros(N,Nu); 
    for i=1:Nu
        M(i:N,i)=s_ob(1:(N-i+1));
    end
    
    
    % macierz Mp obiektu
    MP=zeros(N, D-1);
    po1 = 1; % pomocnicza1
    po2 = 1; % pomocnicza2
    for j=1:D-1
        for i=1:N
    
            if po1+1 >= D %   po pozycji D w macierszy s przyjmujemy ze s=s(D)
                po1 = D-1;
            end
    
            MP(i,j)=(s_ob(po1+1)-s_ob(po2));
            po1 = po1+1;
        end
        po2 = po2+1;
        po1 = j + 1;
    end
    
    % Macierz MZP zakłóceń
    
    MZP = zeros(N, DZ);
    MZP(:,1) = s_zak(1:N)';
    po1 = 2; % pomocnicza1
    po2 = 1; % pomocnicza2
    for i=1:N 
        for j=2:DZ
            MZP(i,j) = s_zak(po1)-s_zak(po2);
            po1 = po1+1;
            po2 = po2+1;
        end
        po1 = i+1;
        po2 = 1;
    end

        I=eye(Nu);
        K=((M'*M+lambda(r)*I)^-1)*M';
        ku(r,:)=K(1,:)*MP;
        kz(r,:)=K(1,:)*MZP;
        ke(r)=sum(K(1,:));
    end

    %%
    deltaup=zeros(D-1,1);
    deltazp=zeros(DZ, 1);
    deltauk = zeros(il,1);
    w = zeros(1,il);

    %% dane
    C1 = 0.95;
    C2 = 0.95;
    a1 = 17;
    a2 = 15;
    n = 12100;
    tau = 40;

    U0 = 52;
    D0 = 12;
    Y0 = 18.2044;
    V2 = C2*Y0^2;
    h1 = 14.1730;
    V1 = C1 * h1^2;
    start = 100;
    Umax = 150;
    Umin = 0;


    U = U0*ones(1,n);
    z = D0*ones(1,n);
    Y = Y0*ones(1,n);


    Yzad(1:start) = 18.2044; Yzad(start:3000) = 50; Yzad(3001:6000) = 10; Yzad(6000:n) = 100;
    z(1500:4500) = 5; z(4501:7000) = 20; z(7001:10000) = 13; z(10001:n) = 25; 
    e = zeros(1,n);


    %%
    hold on
    for k = start:n
        % symulacja
        V1 = V1 + U(k-1-tau) + z(k-1) - a1*h1^0.5;
        V2 = V2 + a1*h1^0.5 - a2*Y(k-1)^0.5;
        h1 = (V1/C1)^0.5;
        Y(k) = (V2/C2)^0.5;
        % uchyb
        e(k) = Yzad(k) - Y(k);

        % Prawo regulacji
        for i = 1:il

            if pomiar_z
                deltauk(i) = ke(i)*e(k)-ku(i,:)*deltaup-kz(i,:)*deltazp;
            else
                deltauk(i) = ke(i)*e(k)-ku(i,:)*deltaup;
            end
            if i == 1
                w(i) = 1-1/(1+exp(-a*(Y(k)-c(1))));
            elseif i == il
                w(i) = 1/(1+exp(-a*(Y(k)-c(il-1))));
            else
                w(i) = 1/(1+exp(-a*(Y(k)-c(i-1)))) - 1/(1+exp(-a*(Y(k)-c(i))));
            end
        end
        % oblicznie deltyuk użytej do sterowania
        DELTAuk = sum(w*deltauk)/sum(w);

        

        
        % "przesuwanie" macierzy deltaup i deltazp
        for i = D-1:-1:2
          deltaup(i) = deltaup(i-1);
        end
        deltaup(1) = DELTAuk;

        for i = DZ:-1:2
           deltazp(i) = deltazp(i-1);
        end
        deltazp(1) = z(k) - z(k-1);
        
        % uwzględnianie ograniczeń
        if U(k-1)+deltaup(1) < Umax || U(k-1)+deltaup(1) > Umin
            U(k) = U(k-1)+deltaup(1);
        else
            U(k) = U(k-1);
            deltaup(1) = 0;
        end

    end

    E = sum(e.^2);

    if draw
        figure
        subplot(3,1,1)
        plot(1:n, Yzad, 'r')
        xlabel('t[s]')
        ylabel('h2[cm]')
        hold on
        plot(1:n, Y, 'b')
        subplot(3,1,2)
        plot(1:n, U)
        xlabel('t[s]')
        ylabel('F1in[cm3/s]')
        subplot(3,1,3)
        plot(1:n, z)
        xlabel('t[s]')
        ylabel('FD[cm3/s]')
    end

end