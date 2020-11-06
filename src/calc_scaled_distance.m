function scaled_distance = calc_scaled_distance(params_em)

numpos2 = size(params_em,2);

S_11 = params_em(5,:,:,:);
S_12 = params_em(6,:,:,:);
S_13 = params_em(7,:,:,:);
S_22 = params_em(8,:,:,:);
S_23 = params_em(9,:,:,:);
S_33 = params_em(10,:,:,:);

detS =    S_11.*S_22.*S_33 ...
    + 2*S_12.*S_13.*S_23 ...
    -   S_11.*S_23.*S_23 ...
    -   S_12.*S_12.*S_33 ...
    -   S_13.*S_22.*S_13;

detS_inv = 1./detS;

Sinv_11 = detS_inv.*(S_22.*S_33 - S_23.*S_23);
Sinv_22 = detS_inv.*(S_11.*S_33 - S_13.*S_13);
Sinv_33 = detS_inv.*(S_11.*S_22 - S_12.*S_12);
Sinv_12 = detS_inv.*(S_13.*S_23 - S_12.*S_33);
Sinv_13 = detS_inv.*(S_12.*S_23 - S_13.*S_22);
Sinv_23 = detS_inv.*(S_12.*S_13 - S_11.*S_23);


scaled_distance = zeros(numpos2,numpos2);
for p=1:numpos2
    
    xd2 = bsxfun(@minus,params_em(2,:),params_em(2,p));
    yd2 = bsxfun(@minus,params_em(3,:),params_em(3,p));
    zd2 = bsxfun(@minus,params_em(4,:),params_em(4,p));
    
    xp12 = Sinv_11(p)*xd2 + Sinv_12(p)*yd2 + Sinv_13(p)*zd2;
    yp12 = Sinv_12(p)*xd2 + Sinv_22(p)*yd2 + Sinv_23(p)*zd2;
    zp12 = Sinv_13(p)*xd2 + Sinv_23(p)*yd2 + Sinv_33(p)*zd2;
    xp21 = Sinv_11  .*xd2 + Sinv_12  .*yd2 + Sinv_13  .*zd2;
    yp21 = Sinv_12  .*xd2 + Sinv_22  .*yd2 + Sinv_23  .*zd2;
    zp21 = Sinv_13  .*xd2 + Sinv_23  .*yd2 + Sinv_33  .*zd2;
    
    scaled_distance(:,p) =   (xp12.*xd2 + yp12.*yd2 + zp12.*zd2)...
        .*(xp21.*xd2 + yp21.*yd2 + zp21.*zd2);
    
end

end
