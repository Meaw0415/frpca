function [Um, Sm,Vm] = frpca_subspace_merge2(U1, S1, U2, S2)
        [Qp, Rp] = qr([(U1*S1),U2*S2], 0);
        [Ur, Sr, Vm] = svd(Rp);
        qq = Qp;
        uu = Ur;
        Um = qq*uu;
        Sm = Sr;


        
