function [E] = ZAD1_DMC(D, N, Nu, DZ, lambda, draw)

    % zmienne i macierze regulatora
    [s, z] = ZAD1_odp_skokowa(18.2044, false);
    D = min(D, length(s));
    DZ = min(DZ, length(z));
    N = min(min(N,D),DZ);
    Nu = min(Nu, N);

    M=zeros(N,Nu);
    for i=1:N
       for j=1:Nu
          if (i>=j)
             M(i,j)=s(i-j+1);
          end
       end
    end

    MP=zeros(N,D-1);
    for i=1:N
       for j=1:D-1
          if i+j<=D
             MP(i,j)=s(i+j)-s(j);
          else
             MP(i,j)=s(D)-s(j);
          end    
       end
    end

    MZP=zeros(N,DZ);
    for i=1:N
        MZP(i,1) = z(i);
       for j=2:DZ
          if i+j-1<=DZ
             MZP(i,j)=z(i+j-1)-z(j);
          else
             MZP(i,j)=z(DZ)-z(j);
          end      
       end
    end

    I=eye(Nu);
    K=((M'*M+lambda*I)^-1)*M';
    ku=K(1,:)*MP;
    kz=K(1,:)*MZP;
    ke=sum(K(1,:));
    deltaup=zeros(1,D-1);
    deltazp=zeros(1,DZ-1);

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
    dY = [start 50; start+3000 20; start+6000 100; start+9000 100];
    dZ = [start+1500 12; start+4500 5; start+7500 12.5; start+10500 7.5];

    U = U0*ones(1,n);
    Dist = D0*ones(1,n);
    Y = Y0*ones(1,n);
    Yz = Y;
    for i = 1:length(dY)
        Yz(dY(i,1):n) = dY(i,2);
        Dist(dZ(i,1):n) = dZ(i,2);
    end
    e = zeros(1,n);

    for k = start:n
        % symulacja
        V1 = V1 + U(k-1-tau) + Dist(k-1) - a1*h1^0.5;
        V2 = V2 + a1*h1^0.5 - a2*Y(k-1)^0.5;
        h1 = (V1/C1)^0.5;
        Y(k) = (V2/C2)^0.5;
        e(k) = Yz(k) - Y(k);

        %uwzglêdnianie zak³ócenia
        for i = DZ:-1:2
           deltazp(i) = deltazp(i-1);
        end
        deltazp(1) = Dist(k) - Dist(k-1);

        % Prawo regulacji
        deltauk = ke*e(k)-ku*deltaup'-kz*deltazp';

        for i = D-1:-1:2
          deltaup(i) = deltaup(i-1);
        end
        deltaup(1) = deltauk;
        U(k) = U(k-1)+deltaup(1);
    end

    if draw
        subplot(3,1,1)
        plot(1:n, Yz, 'r')
        xlabel('t[s]')
        ylabel('h2[cm]')
        hold on
        plot(1:n, Y, 'b')
        subplot(3,1,2)
        plot(1:n, U)
        xlabel('t[s]')
        ylabel('F1in[cm3/s]')
        subplot(3,1,3)
        plot(1:n, Dist)
        xlabel('t[s]')
        ylabel('FD[cm3/s]')
    end
    E = sum(e.^2);
end