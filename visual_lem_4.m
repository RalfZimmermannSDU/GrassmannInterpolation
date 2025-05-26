clear
close all
rng(123123,'twister')
n = 1000;
p = 20;

M = matrix_tools();

U = M.RandG(n,p);
Delta = M.vectorG(U);
Delta = Delta/(norm(Delta,'fro')*0.5);

B = M.LocalCoordG(U,n,p);

man = [];
eucl = [];
loc_eucl = [];
for i = 1:100
    t = (i)/100 * 50 * 1/2;
    man(i) = t;
    V = M.ExpG(U,Delta,t);
    Btilde = M.LocalCoordG(V,n,p);
    eucl(i) = norm(U-V,'fro');
    loc_eucl(i) = norm(B - Btilde,'fro');
end

f = figure;
f.Position = [40,800,1200*5/6,650*5/6];

subplot(1,2,1)
semilogy(man,loc_eucl,'LineWidth',2)
grid on
hold on
xlabel("Manifold distance")
ylabel("Distance between local coords")

title("Manifold dist vs. local coord dist.")

[seucl,I] = sort(eucl,'ascend');

subplot(1,2,2)
semilogy(seucl,loc_eucl(I),'LineWidth',2);
grid on
xlabel("Euclidean disrance")
ylabel("Distance between local coords")
title("Euclidean dist vs. local coord dist.")

fontsize(f,18,"pixels")

exportgraphics(f,"lemma_4.png","Resolution",300);