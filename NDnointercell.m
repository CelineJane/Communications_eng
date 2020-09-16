function [MRC ZF MMSE] = NDnointercell(d,link)
    
    disp('Running simulation ....')
    Nt = 4 ;
    N0 = 1; %watt
    pl = 3; %path loss exponent
    shadow_var = 8; %shadow variance
    
    [users num_users] = drop_users();
    Nr = num_users ; %assume all users haveone antenna each

    rho = 0.9;
    Rt = toeplitz(rho.^[0:Nt-1]);     % tx antenna correlation matrix
    Rr = toeplitz(rho.^[0:Nr-1]);
    H_uncorr  = sqrt(1/2)*(randn(Nt,Nr)+1i*randn(Nt,Nr)); %rayleigh fading
    H = sqrtm(Rt)*H_uncorr*sqrtm(Rr); % assume some correlation for non-distributed array

    Pt_ul_dB = linspace(-20,20,19);
    Pt_ul = 10.^(Pt_ul_dB./10); %watts, assume each user has same transmit power for uplink

    Pt_dl_dB = linspace(-20,20,19);
    Pt_dl = 10.^(Pt_dl_dB./10); %watts, assume each user has same transmit power for uplink


    for SNR = 1:length(Pt_dl)
        if contains(link,'ul')
            Pr_ul_dB(:,SNR) = get_Pr_ul(users,num_users,pl,shadow_var,Pt_ul_dB(SNR));
            Pr_ul(:,SNR) = 10.^(Pr_ul_dB(:,SNR)./10);
            
            G = zeros(Nt,num_users) ; % channel matrix for all users
            for k = 1:num_users
                G(:,k) = sqrt(Pr_ul(k,SNR))*H(:,k);
            end
            SINR_ul_mrc(:,SNR) = get_SINR('mrc',H,G,Pt_ul(SNR),Pr_ul(:,SNR),N0);
            R_ul_mrc(:,SNR) = log2(1+SINR_ul_mrc(:,SNR)); %uplink spectral efficiency
            MRC = sum(R_ul_mrc);

            SINR_ul_zf(:,SNR) = get_SINR('zf',H,G,Pt_ul(SNR),Pr_ul(:,SNR),N0);
            R_ul_zf(:,SNR) = log2(1+SINR_ul_zf(:,SNR)); %uplink spectral efficiency
            ZF = sum(R_ul_zf);

            SINR_ul_mmse(:,SNR) = get_SINR('mmse',H,G,Pt_ul(SNR),Pr_ul(:,SNR),N0);
            R_ul_mmse(:,SNR) = log2(1+SINR_ul_mmse(:,SNR)); %uplink spectral efficiency
            MMSE = sum(R_ul_mmse);            
            
        end

        if contains(link,'dl')
            Pr_dl_dB(:,SNR) = get_Pr_ul(users,num_users,pl,shadow_var,Pt_dl_dB(SNR));
            Pr_dl(:,SNR) = 10.^(Pr_dl_dB(:,SNR)./10);
            
            SINR_dl_mrc(:,SNR) = get_SINR_dl('mrc',H,Pt_dl(SNR),Pr_dl(:,SNR),N0);
            R_dl_mrc(:,SNR) = log2(1+SINR_dl_mrc(:,SNR)); %uplink spectral efficiency
            MRC = sum(R_dl_mrc);

            SINR_dl_zf(:,SNR) = get_SINR_dl('zf',H,Pt_dl(SNR),Pr_dl(:,SNR),N0);
            R_dl_zf(:,SNR) = log2(1+SINR_dl_zf(:,SNR)); %uplink spectral efficiency
            ZF = sum(R_dl_zf);
        %     
            SINR_dl_mmse(:,SNR) = get_SINR_dl('mmse',H,Pt_dl(SNR),Pr_dl(:,SNR),N0);
            R_dl_mmse(:,SNR) = log2(1+SINR_dl_mmse(:,SNR)); %uplink spectral efficiency
            MMSE = sum(R_dl_mmse);            
        end

    end
    disp('Done')


    %creates hexagonal cell array (7, one central cell + 6 surrounding cells)
    %and drops users 
    function [users,num_users] = drop_users() 
        angles = linspace(0, 360, 7);
        r = 1 ;%km
        x = r * cosd(angles); % Assign x coordinates
        y = r * sind(angles); % Assign y coordinates
        figure
        hold on
        plot(x, y, 'b-', 'LineWidth', 2);
        plot(x, y+2*max(y), 'b-', 'LineWidth', 2);
        plot(x, y-2*max(y), 'b-', 'LineWidth', 2);
        plot(x-1.5*max(x), y-max(y), 'b-', 'LineWidth', 2);
        plot(x+1.5*max(x), y-max(y), 'b-', 'LineWidth', 2);
        plot(x+1.5*max(x), y+max(y), 'b-', 'LineWidth', 2);
        plot(x-1.5*max(x), y+max(y), 'b-', 'LineWidth', 2);
        axis equal

        a = r*3*sqrt(3)/2; %area of hexagon formula
        density = d*a;

        num_users = poissrnd(density);
        users = r*sqrt((rand(num_users,1))).*exp(j*2*pi*(rand(num_users,1)));

        [in on] = inpolygon(real(users),imag(users),x,y);

        while ~all(in) & ~all(on) 
            users(~in) = r*sqrt((rand(length(users(~in)),1))).*exp(j*2*pi*(rand(length(users(~in)),1)));
            [in on] = inpolygon(real(users),imag(users),x,y);

        end
        plot(users(in),'.')
        hold off

    end

    % received power at BS from kth user
    function Pr_ul = get_Pr_ul(users,num_users,pl,var,Pt) 
    %     PL = 10*pl*log10(abs(users)) + sqrt(var)*randn(num_users,1)
         PL = -10*pl*log10(abs(users)) + sqrt(var)*randn(num_users,1)   ;
         Pr_ul = Pt - abs(PL);
    end

    function SINR = get_SINR(precoder,H,G,Pt,Pr,N0)
        if contains(precoder,'mrc')
            a = H'; 
            intfrnce = zeros(length(H)-1, length(H));

            for n = 1:length(H) % user that is receiveing inteference from others
                for k = 1:length(H) %users giving interference
                    if n ~= k
                        intfrnce(k,n) = Pr(k)*(abs(H(:,n)'*H(:,k))^2); 
                    end
                end
                SINR(n) = (Pt*Pr(n)*(norm(H(:,n))^4))/((N0*(norm(H(:,n))^2)+(Pt*sum(intfrnce(:,n)))));
            end

        elseif contains(precoder,'zf')
    %         a = (inv(G'*G)*G')';
            a = (G'*inv(G*G'))'; %TAKE RIGHT INVERSE NOT LEFT !!!
            intfrnce = zeros(length(H)-1, length(H));

            for n = 1:length(H) % user that is receiveing inteference from others
                for k = 1:length(H) %users giving interference
                    if n ~= k
                        intfrnce(k,n) = (abs(a(:,n)'*G(:,k)))^2;
                    end
                end
                (Pt*(norm(a(:,n).*G(:,n)))^2);
                SINR(n) = (Pt*(abs(a(:,n)'*G(:,n)))^2)/((N0*(norm(a(:,n))^2)+(Pt*sum(intfrnce(:,n)))));
            end
            
        elseif contains(precoder,'mmse')

            for n = 1:length(H)
                intfrnce = zeros(size(H,1),size(H,1));
                for k = 1:length(H)
                    if n ~= k
                        intfrnce = intfrnce + H(:,k)*H(:,k)';
                    end   

                end
                SINR(n) = Pt.*Pr(n).*H(:,n)'*(inv((Pt.*Pr(n).*intfrnce)+eye(size(H,1))))*H(:,n);
            end
        end

    end


    function SINR = get_SINR_dl(precoder,H,Pt,Pr,N0)
        if contains(precoder,'mrc')
            w = H;
        elseif contains(precoder,'zf')
    %         w = (inv(H'*H)*H')';
            w = (H'*inv(H*H'))';
        elseif contains(precoder,'mmse')
            for x = 1:length(H)
                w(:,x) = H(:,x)*inv(H(:,x)'*H(:,x) + (size(H(:,x),2)/(Pr(x)*Pt/N0)));
    %         w = H*inv(H'*H + (size(H,2)/(Pt/N0)));
            end
        end
        nu = 1/(mean(trace(w*w'))) ;

        intfrnce = zeros(length(H)-1, length(H));
        for n = 1:length(H) % user that is receiveing inteference from others
            for k = 1:length(H) %users giving interference
                if n ~= k
                    intfrnce(k,n) = (abs(H(:,n)'*w(:,k))^2); 
                end
            end

            SINR(n) = (Pt*nu*Pr(n)*(norm(H(:,n)'*w(:,n))^2))/(((N0^2)+(Pt*nu*Pr(n)*sum(intfrnce(:,n)))));
        end

    end


    
end
