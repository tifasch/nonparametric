function [z1,z2,Faktoren,Splinegrad,Knopt,t] = FactorScoresEst(X,Y,AnzManVar,DimUnVar,AnzahlMessungenX)


% In dieser Funktion werden die eingegebenen Daten nicht normiert. Sollte
% dies für Ihre Zwecke notwendig sein, normieren Sie sie bitte vorher.

%Der Datensatz wird in einer Datei eingelesen und in X und Y vom Programm
%aufgeteilt

%AnzahlMessungen ist im Paper n und stimmt in X und Y überein
dimY = AnzManVar(end);


%% %Conditions fulfilled of the Paper?
for it = 1:DimUnVar
    if AnzManVar(it) < 3
        disp('Dimension of one latent variable is too small');
        return
    end
end
if dimY < 3
    disp('Dimension of Y is too small');
    return
end

%% 
AnzManVar_1 = [0 AnzManVar];
a = zeros(sum(AnzManVar(1:(end-1))),1);
for it = 1:DimUnVar
    a(sum(AnzManVar_1(1:it))+1,1) = 1;
    a(sum(AnzManVar_1(1:it))+2,1) = ((1/AnzManVar_1(it+1))*sum(X{it}(2,:)*X{it}(3,:)'))/((1/AnzManVar_1(it+1))*sum(X{it}(1,:)*X{it}(3,:)'));
    for indexa=3:AnzManVar(it)
        a(sum(AnzManVar_1(1:it))+indexa,1) = ((1/AnzManVar_1(it+1))*sum(X{it}(2,:)*X{it}(indexa,:)'))/((1/AnzManVar_1(it+1))*sum(X{it}(1,:)*X{it}(2,:)'));
    end
end
b = zeros(dimY,1);
b(1,1) = 1;
b(2,1) = ((1/dimY)*sum(Y(2,:)*Y(3,:)'))/((1/dimY)*sum(Y(1,:)*Y(3,:)'));
for indexb=3:dimY
    b(indexb,1) = ((1/dimY)*sum(Y(2,:)*Y(indexb,:)'))/((1/dimY)*sum(Y(1,:)*Y(2,:)'));
end
%%
 
Nn = round(AnzahlMessungenX^(1/3));
rho(1:Nn,1) = 1/Nn;
alpha(1:Nn,1:(DimUnVar+1)) = randn(Nn,(DimUnVar+1));
beta = randn(Nn,sum(AnzManVar(1:(end-1))));
gamma = randn(Nn,dimY);

x0 = zeros((DimUnVar+1)*AnzahlMessungenX,1);
for it = 1:DimUnVar
    x0((it-1)*AnzahlMessungenX+1:(it)*AnzahlMessungenX,1) = (X{it}'*(1./a(sum(AnzManVar_1(1:it))+1:sum(AnzManVar_1(1:(it+1))),1)))/AnzManVar(it);
end
x0(DimUnVar*AnzahlMessungenX+1:(DimUnVar+1)*AnzahlMessungenX) = (Y'*(1./b))/AnzManVar(end);

options = optimset('Algorithm','interior-point','MaxFunEvals',4500,'Display','off');

[z,fval] = fmincon(@myfun,x0,[],[],[],[],[],[],@mycon,options);

%%
function [c,ceq]=mycon(x)
    %c ist ein Vektor von Skalaren an der Stelle x ausgewertet
    for it1 = 1:DimUnVar
        c(it) = (1/AnzahlMessungenX)*sum(((x((it-1)*AnzahlMessungenX+1:it*AnzahlMessungenX)).^2)) - 1 -...
            (1/AnzahlMessungenX)*sum(((X{it}(1,:)).^2));
    end
    c(it+1) = (1/AnzahlMessungenX)*sum(((x((it*AnzahlMessungenX+1):...
        (it+1)*AnzahlMessungenX)).^2)) - 1 - ...
        (1/AnzahlMessungenX)*sum(((Y(1,:)).^2));
    c(it+2) = (1/AnzahlMessungenX)*sum(((x((it*AnzahlMessungenX+1):...
        (it+1)*AnzahlMessungenX)).^4)) - 1 - (64/AnzahlMessungenX)*...
        sum(((Y(1,:)).^4)) - 72*(((1/AnzahlMessungenX)*sum(Y(1,:)))^4);
    ceq = [];        
end
%%

function Zielfunktionswert = myfun(x)
    Teil1 = 0;
    MatrixAb = zeros(AnzahlMessungenX,sum(AnzManVar_1));
    MatrixAc = zeros(AnzahlMessungenX,sum(dimY));
for r = 1:Nn
    TeilA = 0;
    TeilAa = ones(AnzahlMessungenX,1);
    for i = 1:AnzahlMessungenX
        for it2 = 1:DimUnVar+1
            TeilAa(i) = TeilAa(i)*1/(1+exp(AnzahlMessungenX*(x(AnzahlMessungenX*(it2-1)+i,1)-alpha(r,it2))));
        end
        for j1 = 1:DimUnVar
            for j2 = 1:AnzManVar(j1)
                MatrixAb(i,sum(AnzManVar_1(1:j1))+j2) = 1/(1+exp(AnzahlMessungenX*((X{j1}(j2,i)-a(sum(AnzManVar_1(1:j1))+j2,1)*x(AnzahlMessungenX*(j1-1)+i,1))-beta(r,sum(AnzManVar_1(1:j1))+j2)))); 
            end
        end
        for k = 1:dimY % #ok<ALIGN>
            MatrixAc(i,k) = 1/(1+exp(AnzahlMessungenX*((Y(k,i)-b(k,1)*x(end-AnzahlMessungenX+i,1))-gamma(r,k))));
        end
        TeilAb = prod(MatrixAb(i,:));
        TeilAc = prod(MatrixAc(i,:));
        TeilA = TeilA + TeilAa(i) * TeilAb * TeilAc;
    end
 
    TeilA = TeilA/AnzahlMessungenX;
    TeilBb = 1;
    TeilBc = 1;
    
    for j = 1:sum(AnzManVar(1:(end-1))) % #ok<ALIGN>
       TeilBb = TeilBb * (1/AnzahlMessungenX) * sum(MatrixAb(:,j)); 
    end
    for k = 1:dimY % #ok<ALIGN>
       TeilBc = TeilBc * (1/AnzahlMessungenX) * sum(MatrixAc(:,k));
    end
    
    TeilBa = sum(TeilAa)/AnzahlMessungenX;
    TeilB = TeilBa * TeilBb * TeilBc;
    Teil1a = (abs(TeilA-TeilB))^2;
    Teil1 = Teil1 + Teil1a*rho(r);

    clear MatrixAb MatrixAc
end
    
Teil2 = 0;
Teil3 = 0;

for j1 = 1:DimUnVar
    for j2 = 1:AnzManVar(j1)
    Teil2a = 0;
    for i = 1:AnzahlMessungenX
        Teil2a = Teil2a + (X{j1}(j2,i)-a(sum(AnzManVar_1(1:j1))+j2,1)*x(AnzahlMessungenX*(j1-1)+i,1));
    end
    Teil2 = Teil2 + (Teil2a*(1/AnzahlMessungenX))^2;
    end
end

for k = 1:dimY
    Teil3a = 0;
    for i = 1:AnzahlMessungenX
        Teil3a = Teil3a + (Y(k,i) - b(k,1)*x(end-AnzahlMessungenX + i,1));
    end
    Teil3 = Teil3 + (Teil3a*(1/AnzahlMessungenX))^2;
end

Teil = Teil1 + Teil2 + Teil3;


Zielfunktionswert = Teil;
end


clear alpha beta gamma
Zielfunktionswerte(1,1) = fval;

%%


alpha(1:Nn,1:(DimUnVar+1)) = randn(Nn,(DimUnVar+1));
beta = randn(Nn,sum(AnzManVar(1:(end-1))));
gamma = randn(Nn,dimY);
Zielfunktionswerte(2,1) = myfun(z);

Veraenderung = Zielfunktionswerte(1,1) - Zielfunktionswerte(2,1);
Kriterium = Zielfunktionswerte(1,1)/500;

AnzahlIterationen = 2;

while AnzahlIterationen <= 10 && Veraenderung > Kriterium        
    
    [z,fval] = fmincon(@myfun,z,[],[],[],[],[],[],@mycon,options);
    
    Zielfunktionswerte(AnzahlIterationen+1,1) = fval;
        
    alpha(1:Nn,1:(DimUnVar+1)) = randn(Nn,(DimUnVar+1));
    beta = randn(Nn,sum(AnzManVar(1:(end-1))));
    gamma = randn(Nn,dimY);
    
    Zielfunktionswerte(AnzahlIterationen,2) = myfun(z);
    
    Veraenderung = abs(Zielfunktionswerte(AnzahlIterationen,1) - ...
        Zielfunktionswerte(AnzahlIterationen+1,1));
    Kriterium = Zielfunktionswerte(AnzahlIterationen,1)/500;
    
    AnzahlIterationen = AnzahlIterationen + 1;
end
    


z1 = cell(1,DimUnVar);
for it3 = 1:DimUnVar
    z1{it3} = z((AnzahlMessungenX*(it3-1)+1):(AnzahlMessungenX*it3),1);
end
z2 = z((end-AnzahlMessungenX+1):end,1);


%%
for ind = 1:DimUnVar
    for indexj = 1:AnzManVar(ind)
        for indexi = 1:AnzahlMessungenX
            epsilon(sum(AnzManVar_1(1:ind))+indexj, indexi) = X{ind}(indexj,indexi) - a(sum(AnzManVar_1(1:ind))+indexj,1)*z1{ind}(indexi,1);
        end
    end
end
for indexk = 1:dimY
    for indexi = 1:AnzahlMessungenX
        delta(indexk, indexi) = Y(indexk,indexi)-b(indexk,1)*z2(indexi,1);
    end
end
%%

[Faktoren,Splinegrad,Knopt,t] = BSplineEst(z1,z2,AnzManVar,DimUnVar,AnzahlMessungenX);


end



