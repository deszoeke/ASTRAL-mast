function [fincombined, fobscombined, fcr, fin1, fobs1] = Collins_solutions(U,thetar,g,fobs)

if thetar <=90  %case steming into waves
    u=U.*cos(deg2rad(thetar));
    k1=[g + 4*pi.*fobs.*u - sqrt(g.^2 + 8*pi.*g.*fobs.*u)]./[2*u.^2];
    k1(abs(k1)>1e5)=0;   %fix some spikes due to problems with dividing by 0 (some seems ok, some not)
%     k1=despike2(k1);   %fix some spikes due to problems with dividing by 0 (some seems ok, some not)
    fin1=fobs - k1.*u./(2*pi); %intrinsic frequency
    
    fincombined=fin1;
    fobscombined=fobs;
%     fcr=NaN;
    fcr=[max(fin1) fobs(fin1==(max(fin1)))];
    fobs1=fobs;
elseif thetar>90  %case of steaming with waves
    
    u=U.*cos(deg2rad(thetar));
    %solution for long, fast waves (Cp faster than twice the ship speed)
    k1=[g + 4*pi.*fobs.*u - sqrt(g.^2 + 8*pi.*g.*fobs.*u)]./[2*u.^2];
    fin1=fobs - k1.*u./(2*pi); %intrinsic frequency
    %solution for intermediated speeds
    k2=[g + 4*pi.*fobs.*u + sqrt(g.^2 + 8*pi.*g.*fobs.*u)]./[2*u.^2];
    fin2=fobs - k2.*u./(2*pi); %intrinsic frequency
    %solution for vessel overtaking the waves
    k3=[g - 4*pi.*fobs.*u + sqrt(g.^2 - 8*pi.*g.*fobs.*u)]./[2*u.^2];
    fin3=-fobs - k3.*u./(2*pi); %intrinsic frequency
    
    if length(fin1)>1
        toto1=find(imag(fin1)==0);  %find indices of real values and not imaginary numbers
        toto2=find(imag(fin2)==0);
        if length(toto1)==length(fin1)  %case mostly when SOG~0 and all freuencies are in
            fincombined=real([fin1(toto1(1):toto1(end)) flip(fin2(toto2(1):toto2(end))) fin3]);
            fobscombined=[fobs(toto1(1):toto1(end)) flip(fobs(toto2(1):toto2(end))) fobs];
            fcr=[real(fin1(toto1(end))) fobs(toto1(end))]; %critical intrinsic frequency where there is no ambiguity for fobs > fcr , and its associated fobs
        else
            fincombined=real([fin1(toto1(1):toto1(end)+1) flip(fin2(toto2(1):toto2(end)+1)) fin3]);
            fobscombined=[fobs(toto1(1):toto1(end)+1) flip(fobs(toto2(1):toto2(end)+1)) fobs];
    %         fcr=real(fin1(toto1(end)+1)); %critical frequency where there is no ambiguity for fobs > fcr 
            fcr=[real(fin1(toto1(end)+1)) fobs(toto1(end)+1)]; %critical intrinsic frequency where there is no ambiguity for fobs > fcr , and its associated fobs
        end    
        fin1=real(fin1(toto1));%save the intrinsic frequencies od loweset frequency branch only as the ones used to correct DC
        fobs1=fobs(toto1);
    else
        fincombined=real([fin1 fin2 fin3]);
        fobscombined=[fobs fobs fobs];
        fcr=[NaN NaN];
        fin1=real(fin1);
        fobs1=fobs;
    end

end    