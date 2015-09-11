clear variables
Npol = 8;
Ntor = 8;

Lt = 0.001;
Lp = 0.001;

Pp = 20;

N = Npol*Ntor;
A = Lt*Lp;

Tij = zeros(N,N);

tracker_param

for i2=1:Npol
    i2
    for j2=1:Ntor
         for k2 = 1:Pp
             k2;
            thomp
            
            p = (i2-1+rand)*Lp;
            t = (j2-1+rand)*Lt;
            r = [p,1e-6,t];
            Ka;

            tracker
           
            r;
            if ((r(1) > 0) && (r(3)>0)) 
            pind = fix(r(1)/Lp) + 1;
            tind = fix(r(3)/Lt) + 1;
            
            Tij((i2-1)*Npol +j2, (pind-1)*Npol +tind) = Tij((i2-1)*Npol +j2, (pind-1)*Npol +tind)+1;
            
            end

           
         end
    end
end
colormap('winter')
imagesc(Tij)
colorbar

map2 = transpose(reshape(Tij(37,:), [Npol,Ntor]));

imagesc(map2)
colorbar