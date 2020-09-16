function [MRC ZF MMSE] = NDintercell(d,link)
  
    disp('Running simulation ....')
    Nt = 4 ;
    N0 = 1; %watt
    pl = 3; %path loss exponent
    shadow_var = 8; %shadow variance

    [users num_users, ant] = drop_users_ant(Nt);
    Nr = num_users ; %assume all users haveone antenna each

    rho = 1; % can assume no correlation between antennas as they are spacced far apart
    Rt = toeplitz(rho.^[0:Nt-1]);     % tx antenna correlation matrix
    Rr = toeplitz(rho.^[0:Nr-1]);
    H_uncorr  = sqrt(1/2)*(randn(Nt,Nr,length(ant))+1i*randn(Nt,Nr,length(ant))); %rayleigh fading
    H = H_uncorr; %since rho = 1

    Pt_ul_dB = linspace(-20,20,19);
    Pt_ul = 10.^(Pt_ul_dB./10); %watts, assume each user has same transmit power for uplink

    Pt_dl_dB = linspace(0,60,19);
    Pt_dl = 10.^(Pt_dl_dB./10); %watts, assume each user has same transmit power for uplink


    for SNR = 1:length(Pt_dl)
        if contains(link,'ul')
            Pr_ul_dB(:,:,:,SNR) = get_Pr_ul(users,num_users,pl,shadow_var,Pt_ul_dB(SNR),ant);
            Pr_ul(:,:,:,SNR) = 10.^(Pr_ul_dB(:,:,:,SNR)./10);            

            G = zeros(Nt,num_users,length(ant)) ; % channel matrix for all users
            for L = 1:length(ant)
                for n = 1:size(ant,1)
                    G(n,:,L) = sqrt(Pr_ul(L,:,n,SNR)).*H(n,:,L);
                end
            end        

            SINR_ul_mrc(:,SNR) = get_SINR_ul('mrc',H,G,Pt_ul(SNR),Pr_ul(:,:,:,SNR),N0);
            R_ul_mrc(:,SNR) = log2(1+SINR_ul_mrc(:,SNR)); %uplink spectral efficiency
            MRC = sum(R_ul_mrc);

            SINR_ul_zf(:,SNR) = get_SINR_ul('zf',H,G,Pt_ul(SNR),Pr_ul(:,:,:,SNR),N0);
            R_ul_zf(:,SNR) = log2(1+SINR_ul_zf(:,SNR)); %uplink spectral efficiency
            ZF = sum(R_ul_zf);

            SINR_ul_mmse(:,SNR) = get_SINR_ul('mmse',H,G,Pt_ul(SNR),Pr_ul(:,:,:,SNR),N0);
            R_ul_mmse(:,SNR) = log2(1+SINR_ul_mmse(:,SNR)); %uplink spectral efficiency
            MMSE = sum(R_ul_mmse);            
        end
        
        if contains(link,'dl')
            Pr_dl_dB(:,:,:,SNR) = get_Pr_dl(users,num_users,pl,shadow_var,Pt_dl_dB(SNR),ant);
            Pr_dl(:,:,:,SNR) = 10.^(Pr_dl_dB(:,:,:,SNR)./10);            
            
            SINR_dl_mrc(:,SNR) = get_SINR_dl('mrc',H,Pt_dl(SNR),Pr_dl(:,:,:,SNR),N0);
            R_dl_mrc(:,SNR) = log2(1+SINR_dl_mrc(:,SNR)); %uplink spectral efficiency
            MRC = sum(R_dl_mrc);

            SINR_dl_zf(:,SNR) = get_SINR_dl('zf',H,Pt_dl(SNR),Pr_dl(:,:,:,SNR),N0);
            R_dl_zf(:,SNR) = log2(1+SINR_dl_zf(:,SNR)); %uplink spectral efficiency
            ZF = sum(R_dl_zf);
        %     
            SINR_dl_mmse(:,SNR) = get_SINR_dl('mmse',H,Pt_dl(SNR),Pr_dl(:,:,:,SNR),N0);
            R_dl_mmse(:,SNR) = log2(1+SINR_dl_mmse(:,SNR)); %uplink spectral efficiency
            MMSE = sum(R_dl_mmse);
        end
 
    end
    disp('Done')

    %creates hexagonal cell array (7, one central cell + 6 surrounding cells)
    %and drops users 
    function [users,num_users, ant] = drop_users_ant(Nt) 
        angles = linspace(0, 360, 7);
        r = 1 ;%km
        r_BS = 0.5;
        x = r * cosd(angles); % Assign x coordinates
        y = r * sind(angles); % Assign y coordinates

        hex_x = [x; x; x; x-1.5*max(x); x+1.5*max(x); x+1.5*max(x); x-1.5*max(x)];
        hex_y = [y ;y+2*max(y) ;y-2*max(y) ;y-max(y); y-max(y); y+max(y); y+max(y)];

        figure
        hold on
        for k = 1:length(hex_x)
            plot(hex_x(k,:), hex_y(k,:), 'b-', 'LineWidth', 2);
        end
        axis equal

        a = r*3*sqrt(3)/2; %area of hexagon formula
        density = d*a;

        num_users = poissrnd(density);
        users = r*sqrt((rand(num_users,7))).*exp(j*2*pi*(rand(num_users,7)));
        midx = (max(hex_x')+min(hex_x'))/2;
        midy = (max(hex_y')+min(hex_y'))/2;

        ant = (r_BS*sqrt((rand(Nt,7))).*exp(j*2*pi*(rand(Nt,7))));

        for k = 1:size(users,2)
            [in on] = inpolygon(real(users(:,k)),imag(users(:,k)),hex_x(1,:),hex_y(1,:));
            while ~all(in) & ~all(on) 
                users(~in,k) = r*sqrt((rand(length(users(~in,k)),1))).*exp(j*2*pi*(rand(length(users(~in,k)),1)));
                [in on] =  inpolygon(real(users(:,k)),imag(users(:,k)),hex_x(1,:),hex_y(1,:));
            end
            users(:,k) = users(:,k) + midx(k) + j*midy(k);
            ant(:,k) = ant(:,k) + midx(k) + j*midy(k);
            plot(users(in,k),'.')
            plot(ant(:,k),'p','MarkerSize',10,'MarkerFaceColor','r')
        end 

        hold off



    end
    % 
    % received power at BS from kth user
    function Pr_ul = get_Pr_ul(users,num_users,pl,var,Pt,ant) 

     for cell = 1:size(ant,1)
        for k = 1:num_users
            PL = -10*pl*log10(abs(ant(cell,:) - users(k,:))) + sqrt(var)*randn(1,1) ;
    %          PL =  -10*pl*log10(abs(users(:,k))) + sqrt(var)*randn(1,1);
             Pr_ul(:,k,cell) = Pt - abs(PL) ;

        end


     end
    end

    function Pr_dl = get_Pr_dl(users,num_users,pl,var,Pt,ant) 
    %     PL = 10*pl*log10(abs(users)) + sqrt(var)*randn(num_users,1)
    for cell = 1:size(ant,1)
        for k = 1:num_users

             PL = -10*pl*log10(abs(ant(cell,:) - users(k,1))) + sqrt(var)*randn(1,1) ;
    %          pause(2)
             Pr_dl(:,k,cell) = Pt - abs(PL)   ; 
        end
    end

    end


    function SINR = get_SINR_ul(precoder,H,G,Pt,Pr,N0)
        if contains(precoder,'mrc')
            for L = 1:size(H,3)
                a(:,:,L) = H(:,:,L)';
            end



        elseif contains(precoder,'zf')
            for L = 1:size(H,3)
                a(:,:,L) = (G(:,:,L)'*inv(G(:,:,L)*G(:,:,L)'));
    %             a(:,:,L) = (G(:,:,L)'*inv(G(:,:,L)*G(:,:,L)'))';
            end

        elseif contains(precoder,'mmse')
            for L = 1:size(H,3)
                for k = 1:length(G)
                    a(k,:,L) =  sqrt(Pt)*inv(Pt*G(:,:,L)*G(:,:,L)' + eye(size(H,1)))*G(:,k,L);
                end
            end

        end

        intfrnce = zeros(length(H), length(H), size(H,3));
        for i = 1:size(H,3)
            for n = 1:length(H) % user that is receiveing inteference from others
                for k = 1:length(H) %users giving interference           
                    if i == 1 & n == k %if the same user 
                        intfrnce(k,n,i) = 0;                                  
                    else
                        intfrnce(k,n,i) = (abs(a(n,:,1)*G(:,k,i))^2);                  
                    end
                end      
                sum_int_BS = sum(intfrnce);
    %             pause(5)
            end
            %sum_int_BS(:,:,i) = sum_int_BS(:,:,i);        
        end
       sum_total = sum(sum_int_BS,3);


       for n = 1:length(H)

           SINR(n) = ((Pt*(norm(a(n,:,1)*G(:,n,1))^2))/((((N0^2)*(norm(a(n,:,1)))^2)+(Pt*sum_total(n)))));
       end



    end
    % 
    function SINR = get_SINR_dl(precoder,H,Pt,Pr,N0)
        if contains(precoder,'mrc')
            w = H;
        elseif contains(precoder,'zf')
    %         w = (inv(H'*H)*H')';
            for k = 1:size(H,3)
%                 w(:,:,k) = (H(:,:,k)'*inv(H(:,:,k)*H(:,:,k)'))';
                 w(:,:,k) = ((inv(H(:,:,k)'*H(:,:,k)))*H(:,:,k)')';
            end
        elseif contains(precoder,'mmse')
            for k = 1:size(H,3)
                for x = 1:length(H)
                    w(:,x,k) =  H(:,x,k)*inv(H(:,x,k)'*H(:,x,k) + (size(H(:,x,k),2)/(Pr(k,x)*Pt/N0)));
    %         w = H*inv(H'*H + (size(H,2)/(Pt/N0)));
                end
            end
        end

        for x = 1:size(H,3)
            nu(x) = 1/(mean(trace(w(:,:,x)*w(:,:,x)')));
        end

        intfrnce = zeros(length(H), length(H), size(H,3));
        for i = 1:size(H,3)
            for n = 1:length(H) % user that is receiveing inteference from others
                for k = 1:length(H) %users giving interference           
                    if i == 1 & n == k %if the same user 
                        intfrnce(k,n,i) = 0;                                  
                    else
                        P(:,k,i) =  Pr(i,k,:);
                        intfrnce(k,n,i) = (abs(H(:,n,1)'*(P(:,k,i).*w(:,k,i)))^2);
                    end
                end      
                sum_int_BS = sum(intfrnce);           
            end

            sum_int_BS(:,:,i) = nu(i).*sum_int_BS(:,:,i);        
        end
       sum_total = sum(sum_int_BS,3);

       for n = 1:length(H)
           SINR(n) = ((Pt*nu(1)*(norm(H(:,n,1)'*(P(:,n,1).*w(:,n,1)))^2))/(((N0^2)+(Pt*sum_total(n)))))';
       end

    end

end

