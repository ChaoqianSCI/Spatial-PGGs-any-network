
%% star graph, n leafs
n=9;
N=n+1;
nei=zeros(N,N);
numnei=zeros(N,1);

nei(1,:)=[1:N];
numnei(1,1)=n;

for i=2:N
    nei(i,1)=i;
    nei(i,2)=1;
    numnei(i,1)=1;
end

w_ij=zeros(N,N);
w_i=zeros(N,1);
for i=1:N
    for j=2:numnei(i,1)+1
        w_ij(i,nei(i,j))=1;
    end
    w_i(i,1)=sum(w_ij(i,:));
end

[r_PC,r_DB,r_BD,r_PC_accu,r_DB_accu,r_BD_accu,bc_PC,bc_DB,bc_BD,bc_PC_accu,bc_DB_accu,bc_BD_accu]=everything_r_bc_accu(N,w_ij)