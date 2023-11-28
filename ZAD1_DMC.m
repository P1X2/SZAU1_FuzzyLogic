function [E] = ZAD1_DMC(D, N, Nu, DZ, lambda, draw, pomiar_z)

    % zmienne i macierze regulatora


    [s_ob, s_zak] = ZAD1_odp_skokowa(18.2044, false);
%%
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
    K=((M'*M+lambda*I)^-1)*M';
    ku=K(1,:)*MP;
    kz=K(1,:)*MZP;
    ke=sum(K(1,:));
    deltaup=zeros(D-1,1);
    deltazp=zeros(DZ, 1);

    % dane
    C1 = 0.95;
    C2 = 0.95;
    a1 = 17;
    a2 = 15;
    n = 12000;
    tau = 40;

    U0 = 52;
    D0 = 12;
    Y0 = 18.2044;
    V2 = C2*Y0^2;
    h1 = 14.1730;
    V1 = C1 * h1^2;
    start = 100;

    U = U0*ones(1,n);
    z = D0*ones(1,n);
    Y = Y0*ones(1,n);

    Yzad(1:start) = 18.2044; Yzad(start:3000) = 50; Yzad(3001:6000) = 10; Yzad(6000:n) = 100;
    z(1500:4500) = 5; z(4501:7000) = 20; z(7001:10000) = 13; z(10001:n) = 25; 
    e = zeros(1,n);

    for k = start:n
        % symulacja
        V1 = V1 + U(k-1-tau) + z(k-1) - a1*h1^0.5;
        V2 = V2 + a1*h1^0.5 - a2*Y(k-1)^0.5;
        h1 = (V1/C1)^0.5;
        Y(k) = (V2/C2)^0.5;
        e(k) = Yzad(k) - Y(k);

        % Prawo regulacji
        if pomiar_z
            deltauk = ke*e(k)-ku*deltaup-kz*deltazp;
        else
            deltauk = ke*e(k)-ku*deltaup;
        end

        % "przesuwanie" macierzy deltaup i deltazp    
        for i = DZ:-1:2
           deltazp(i) = deltazp(i-1);
        end
        deltazp(1) = z(k) - z(k-1);

        for i = D-1:-1:2
          deltaup(i) = deltaup(i-1);
        end

        deltaup(1) = deltauk;
        U(k) = U(k-1)+deltauk;
    end

    if draw
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
    E = sum(e.^2);
end