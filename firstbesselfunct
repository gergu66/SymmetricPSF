function [ FI1 ] = firstbesselfunct(theta1)
    %FIRSTBESSELFUNCT: this part of the code is in-charge of the integral with the second bessel function
    %all this things depend on theta
    %all this things depend on theta
    global M n h k0 k rp thetap;
    theta(:,1) = theta1(:);
    for j = 2:M
        theta(:,j) = asin(n(1)*sin(theta(:,1))/n(j));
    end
    for j = 1:M
        p(:,j,1) = n(j).*cos(theta(:,j));
        p(:,j,2) = cos(theta(:,j))./n(j);
    end
    t(:,1:M-1,:) = 2*p(:,1:M-1,:)./(p(:,1:M-1,:)+p(:,2:M,:));
    r(:,1:M-1,:) = (p(:,1:M-1,:) - p(:,2:M,:))./(p(:,1:M-1,:) + p(:,2:M,:));
    for j = 2:M-1
        beta(:,j) = k(j).*(h(j-1)-h(j)).*cos(theta(:,j)); %calculate all constants beta
    end
    beta(:,1) = 0; %this will not be used
    for l = 1:2
        D(:,l) = 1 + r(:,1,l).*(r(:,2,l).*exp(2*1i*beta(:,2))+r(:,3,l).*exp(2*1i*(beta(:,2)+beta(:,3)))+r(:,4,l).*exp(2*1i*(beta(:,2)+beta(:,3)+beta(:,4))))+r(:,2,l).*(r(:,3,l).*exp(2*1i*beta(:,3))+r(:,4,l).*exp(2*1i*(beta(:,3)+beta(:,4))))+r(:,3,l).*r(:,4,l).*exp(2*1i*beta(:,4))+r(:,1,l).*r(:,2,l).*r(:,3,l).*r(:,4,l).*exp(2*1i*(beta(:,2)+beta(:,4)));
    end
    TN1 = ones(numel(theta1),2);
    for j = 1:1:M-2
        for l = 1:2
            TN1(:,l) = TN1(:,l).* t(:,j,l).*exp(1i*beta(:,j+1));
        end
    end
    %for j = 1:1:M-2
        for l = 1:2
            TN1(:,l) = TN1(:,l).* t(:,M-1,l)./D(:,l); %I don't remember why t(1)
        end
    %end
    psy = h(M-1)*n(M)*cos(theta(:,M)) - h(1)*n(1)*cos(theta1(:));
    %this are the functions to integrate
    FI1 = sqrt(cos(theta1(:))).*sin(theta1(:)).*exp(1i*k0*psy(:)).*TN1(:,2).*sin(theta(:,M)).*besselj(1,k(1).*sin(theta1(:))*rp*sin(thetap)).*exp(1i*k(M)*rp*cos(thetap).*cos(theta(:,M)));
end
