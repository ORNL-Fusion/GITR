           function [Xrefined, Yrefined, Zrefined] = refineXYZ(X,Y,Z,n)
            
            for j=1:n
                Xrefined = zeros(2*length(X),3);
                Yrefined = zeros(2*length(X),3);
                Zrefined = zeros(2*length(X),3);
                for i=1:length(X)
                    A = [X(i,1) Y(i,1) Z(i,1)];
                    B = [X(i,2) Y(i,2) Z(i,2)];
                    C = [X(i,3) Y(i,3) Z(i,3)];
                    
                    AB = B-A;
                    AC = C-A;
                    BC = C-B;
                    
                    norms =[norm(AB) norm(BC) norm(AC)];
                    [maxVal maxInd] = max(norms);
                    if maxInd ==1
                        midPtAB = A + 0.5*AB;
                        Xrefined(2*i-1,1) = A(1);
                        Xrefined(2*i,1) = midPtAB(1);
                        Xrefined(2*i-1,2) = midPtAB(1);
                        Xrefined(2*i,2) = B(1);
                        Xrefined(2*i-1,3) = C(1);
                        Xrefined(2*i,3) = C(1);
                        Yrefined(2*i-1,1) = A(2);
                        Yrefined(2*i,1) = midPtAB(2);
                        Yrefined(2*i-1,2) = midPtAB(2);
                        Yrefined(2*i,2) = B(2);
                        Yrefined(2*i-1,3) = C(2);
                        Yrefined(2*i,3) = C(2);
                        Zrefined(2*i-1,1) = A(3);
                        Zrefined(2*i,1) = midPtAB(3);
                        Zrefined(2*i-1,2) = midPtAB(3);
                        Zrefined(2*i,2) = B(3);
                        Zrefined(2*i-1,3) = C(3);
                        Zrefined(2*i,3) = C(3);
                    elseif maxInd ==2
                        midptBC = B + 0.5*BC;
                        Xrefined(2*i-1,1) = A(1);
                        Xrefined(2*i,1) = A(1);
                        Xrefined(2*i-1,2) = B(1);
                        Xrefined(2*i,2) = midptBC(1);
                        Xrefined(2*i-1,3) = midptBC(1);
                        Xrefined(2*i,3) = C(1);
                        Yrefined(2*i-1,1) = A(2);
                        Yrefined(2*i,1) = A(2);
                        Yrefined(2*i-1,2) = B(2);
                        Yrefined(2*i,2) = midptBC(2);
                        Yrefined(2*i-1,3) = midptBC(2);
                        Yrefined(2*i,3) = C(2);
                        Zrefined(2*i-1,1) = A(3);
                        Zrefined(2*i,1) = A(3);
                        Zrefined(2*i-1,2) = B(3);
                        Zrefined(2*i,2) = midptBC(3);
                        Zrefined(2*i-1,3) = midptBC(3);
                        Zrefined(2*i,3) = C(3);
                    elseif maxInd ==3
                        midptAC = A + 0.5*AC;
                        Xrefined(2*i-1,1) = A(1);
                        Xrefined(2*i,1) = midptAC(1);
                        Xrefined(2*i-1,2) = B(1);
                        Xrefined(2*i,2) = B(1);
                        Xrefined(2*i-1,3) = midptAC(1);
                        Xrefined(2*i,3) = C(1);
                        Yrefined(2*i-1,1) = A(2);
                        Yrefined(2*i,1) = midptAC(2);
                        Yrefined(2*i-1,2) = B(2);
                        Yrefined(2*i,2) = B(2);
                        Yrefined(2*i-1,3) = midptAC(2);
                        Yrefined(2*i,3) = C(2);
                        Zrefined(2*i-1,1) = A(3);
                        Zrefined(2*i,1) = midptAC(3);
                        Zrefined(2*i-1,2) = B(3);
                        Zrefined(2*i,2) = B(3);
                        Zrefined(2*i-1,3) = midptAC(3);
                        Zrefined(2*i,3) = C(3);
                    end
                end
                X = Xrefined;
                Y = Yrefined;
                Z = Zrefined;
            end
            end
