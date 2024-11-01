function [r_PC,r_DB,r_BD,r_PC_accu,r_DB_accu,r_BD_accu,bc_PC,bc_DB,bc_BD,bc_PC_accu,bc_DB_accu,bc_BD_accu]=everything_r_bc_accu(N,w_ij)

w_i=zeros(N,1);
for i=1:N
    w_i(i,1)=sum(w_ij(i,:));
end
p_ij=zeros(N,N);
for i=1:N
    p_ij(i,:)=w_ij(i,:)/w_i(i,1);
end


locanetoij=zeros(N*(N-1)/2,2);%i,j
locaijtone=zeros(N,N);
ne=0;
for i=1:N
    for j=i+1:N
        ne=ne+1;
        locanetoij(ne,1)=i;
        locanetoij(ne,2)=j;
        locaijtone(i,j)=ne;
    end
end

A=zeros(N*(N-1)/2,N*(N-1)/2);x=zeros(N*(N-1)/2,1);b_=zeros(N*(N-1)/2,1);
for i=1:N
    for j=i+1:N
        for k=1:N
            if k~=j
                J=min(k,j);K=max(k,j);
                A(locaijtone(i,j),locaijtone(J,K))=A(locaijtone(i,j),locaijtone(J,K))+p_ij(i,k)/2;
            end
            if k~=i
                I=min(k,i);K=max(k,i);
                A(locaijtone(i,j),locaijtone(I,K))=A(locaijtone(i,j),locaijtone(I,K))+p_ij(j,k)/2;
            end
        end
        A(locaijtone(i,j),locaijtone(i,j))=A(locaijtone(i,j),locaijtone(i,j))-1;
        b_(locaijtone(i,j),1)=1;
    end
end

x=A\(-b_);

tau_ij=zeros(N,N);
for i=1:N
    for j=i+1:N
        tau_ij(i,j)=x(locaijtone(i,j),1);
    end
end
for i=1:N
    for j=1:i-1
        tau_ij(i,j)=tau_ij(j,i);
    end
end






p2_ij=zeros(N,N);
for i=1:N
    for j=1:N
        for l=1:N
            p2_ij(i,j)=p2_ij(i,j)+p_ij(i,l)*p_ij(l,j);
        end
    end
end

p3_ij=zeros(N,N);
for i=1:N
    for j=1:N
        for l=1:N
            for ll=1:N
                p3_ij(i,j)=p3_ij(i,j)+p_ij(i,l)*p_ij(l,ll)*p_ij(ll,j);
            end
        end
    end
end

nume=0;
for i=1:N
    for j=1:N
        nume=nume+w_i(i,1)*p_ij(i,j)*tau_ij(i,j);
    end
end

deno=0;
for i=1:N
    for j=1:N
        part1=(1/(w_i(i,1)+1)^2)*tau_ij(i,j);
        part2=0;
        for l=1:N
            part2=part2+(w_i(i,1)/(w_i(i,1)+1)^2*p_ij(i,l))*(-tau_ij(i,l)+tau_ij(j,l));
        end
        part3=0;
        for l=1:N
            part3_1=0;
            for ll=1:N
                part3_1=part3_1+w_i(l,1)/(w_i(l,1)+1)*p_ij(l,ll)*(-tau_ij(i,ll)+tau_ij(j,ll));
            end
            part3=part3+(w_i(i,1)/(w_i(i,1)+1)*p_ij(i,l))*(1/(w_i(l,1)+1)*(-tau_ij(i,l)+tau_ij(j,l))+part3_1);
        end
        deno=deno+w_i(i,1)*p_ij(i,j)*(part1+part2+part3);
    end
end

r_PC=nume/deno;

nume=0;
for i=1:N
    for j=1:N
        nume=nume+w_i(i,1)*(w_i(i,1)+1)*p_ij(i,j)*tau_ij(i,j);
    end
end

deno=0;
for i=1:N
    for j=1:N
        part1=(1/(w_i(i,1)+1)^2)*tau_ij(i,j);
        part2=0;
        for l=1:N
            part2=part2+(w_i(i,1)/(w_i(i,1)+1)^2*p_ij(i,l))*(-tau_ij(i,l)+tau_ij(j,l));
        end
        part3=0;
        for l=1:N
            part3_1=0;
            for ll=1:N
                part3_1=part3_1+w_i(l,1)/(w_i(l,1)+1)*p_ij(l,ll)*(-tau_ij(i,ll)+tau_ij(j,ll));
            end
            part3=part3+(w_i(i,1)/(w_i(i,1)+1)*p_ij(i,l))*(1/(w_i(l,1)+1)*(-tau_ij(i,l)+tau_ij(j,l))+part3_1);
        end
        deno=deno+w_i(i,1)*(w_i(i,1)+1)*p_ij(i,j)*(part1+part2+part3);
    end
end

r_PC_accu=nume/deno;

nume=0;
for i=1:N
    for j=1:N
        nume=nume+w_i(i,1)*p_ij(i,j)*tau_ij(i,j);
    end
end

deno=0;
for i=1:N
    for j=1:N
        deno=deno+w_i(i,1)*p2_ij(i,j)*tau_ij(i,j)-w_i(i,1)*p_ij(i,j)*tau_ij(i,j);
    end
end

bc_PC=nume/deno;


nume=0;
for i=1:N
    for j=1:N
        nume=nume+w_i(i,1)^2*p_ij(i,j)*tau_ij(i,j);
    end
end

deno=0;
for i=1:N
    for j=1:N
        for l=1:N
            deno=deno+w_i(i,1)^2*p_ij(i,j)*p_ij(i,l)*(tau_ij(j,l)-tau_ij(i,l));
        end
    end
end

bc_PC_accu=nume/deno;



nume=0;
for i=1:N
    for j=1:N
        nume=nume+w_i(i,1)*p2_ij(i,j)*tau_ij(i,j);
    end
end

deno=0;
for i=1:N
    for j=1:N
        part1=(1/(w_i(i,1)+1)^2)*tau_ij(i,j);
        part2=0;
        for l=1:N
            part2=part2+(w_i(i,1)/(w_i(i,1)+1)^2*p_ij(i,l))*(-tau_ij(i,l)+tau_ij(j,l));
        end
        part3=0;
        for l=1:N
            part3_1=0;
            for ll=1:N
                part3_1=part3_1+w_i(l,1)/(w_i(l,1)+1)*p_ij(l,ll)*(-tau_ij(i,ll)+tau_ij(j,ll));
            end
            part3=part3+(w_i(i,1)/(w_i(i,1)+1)*p_ij(i,l))*(1/(w_i(l,1)+1)*(-tau_ij(i,l)+tau_ij(j,l))+part3_1);
        end
        deno=deno+w_i(i,1)*p2_ij(i,j)*(part1+part2+part3);
    end
end

r_DB=nume/deno;


nume=0;
for i=1:N
    for j=1:N
        nume=nume+w_i(i,1)*(w_i(i,1)+1)*p2_ij(i,j)*tau_ij(i,j);
    end
end

deno=0;
for i=1:N
    for j=1:N
        part1=(1/(w_i(i,1)+1)^2)*tau_ij(i,j);
        part2=0;
        for l=1:N
            part2=part2+(w_i(i,1)/(w_i(i,1)+1)^2*p_ij(i,l))*(-tau_ij(i,l)+tau_ij(j,l));
        end
        part3=0;
        for l=1:N
            part3_1=0;
            for ll=1:N
                part3_1=part3_1+w_i(l,1)/(w_i(l,1)+1)*p_ij(l,ll)*(-tau_ij(i,ll)+tau_ij(j,ll));
            end
            part3=part3+(w_i(i,1)/(w_i(i,1)+1)*p_ij(i,l))*(1/(w_i(l,1)+1)*(-tau_ij(i,l)+tau_ij(j,l))+part3_1);
        end
        deno=deno+w_i(i,1)*(w_i(i,1)+1)*p2_ij(i,j)*(part1+part2+part3);
    end
end

r_DB_accu=nume/deno;

nume=0;
for i=1:N
    for j=1:N
        nume=nume+w_i(i,1)*p2_ij(i,j)*tau_ij(i,j);
    end
end

deno=0;
for i=1:N
    for j=1:N
        deno=deno+w_i(i,1)*p3_ij(i,j)*tau_ij(i,j)-w_i(i,1)*p_ij(i,j)*tau_ij(i,j);
    end
end


bc_DB=nume/deno;

nume=0;
for i=1:N
    for j=1:N
        nume=nume+w_i(i,1)^2*p2_ij(i,j)*tau_ij(i,j);
    end
end

deno=0;
for i=1:N
    for j=1:N
        for l=1:N
        deno=deno+w_i(i,1)^2*p2_ij(i,j)*p_ij(i,l)*(tau_ij(j,l)-tau_ij(i,l));
        end
    end
end


bc_DB_accu=nume/deno;

% BD

A1=zeros(N*(N-1)/2,N*(N-1)/2);x1=zeros(N*(N-1)/2,1);b_1=zeros(N*(N-1)/2,1);
for i=1:N
    for j=i+1:N
        denosum=0;
        for k=1:N
            if k~=j
                J=min(k,j);K=max(k,j);
                A1(locaijtone(i,j),locaijtone(J,K))=A1(locaijtone(i,j),locaijtone(J,K))+p_ij(k,i);
            end
            if k~=i
                I=min(k,i);K=max(k,i);
                A1(locaijtone(i,j),locaijtone(I,K))=A1(locaijtone(i,j),locaijtone(I,K))+p_ij(k,j);
            end
            denosum=denosum+p_ij(k,i)+p_ij(k,j);
        end
        denosum=1/denosum;
        A1(locaijtone(i,j),:)=A1(locaijtone(i,j),:)*denosum;
        A1(locaijtone(i,j),locaijtone(i,j))=A1(locaijtone(i,j),locaijtone(i,j))-1;
        b_1(locaijtone(i,j),1)=denosum;
    end
end

x1=A1\(-b_1);

tau1_ij=zeros(N,N);
for i=1:N
    for j=i+1:N
        tau1_ij(i,j)=x1(locaijtone(i,j),1);
    end
end
for i=1:N
    for j=1:i-1
        tau1_ij(i,j)=tau1_ij(j,i);
    end
end


nume=0;
for i=1:N
    for j=1:N
        nume=nume+w_ij(i,j)/(w_i(i,1)*w_i(j,1))*tau1_ij(i,j);
    end
end

deno=0;
for i=1:N
    for j=1:N
        part1=(1/(w_i(i,1)+1)^2)*tau1_ij(i,j);
        part2=0;
        for l=1:N
            part2=part2+(w_i(i,1)/(w_i(i,1)+1)^2*p_ij(i,l))*(-tau1_ij(i,l)+tau1_ij(j,l));
        end
        part3=0;
        for l=1:N
            part3_1=0;
            for ll=1:N
                part3_1=part3_1+w_i(l,1)/(w_i(l,1)+1)*p_ij(l,ll)*(-tau1_ij(i,ll)+tau1_ij(j,ll));
            end
            part3=part3+(w_i(i,1)/(w_i(i,1)+1)*p_ij(i,l))*(1/(w_i(l,1)+1)*(-tau1_ij(i,l)+tau1_ij(j,l))+part3_1);
        end
        deno=deno+w_ij(i,j)/(w_i(i,1)*w_i(j,1))*(part1+part2+part3);
    end
end

r_BD=nume/deno;

nume=0;
for i=1:N
    for j=1:N
        nume=nume+w_ij(i,j)/(w_i(i,1)*w_i(j,1))*(w_i(i,1)+1)*tau1_ij(i,j);
    end
end

deno=0;
for i=1:N
    for j=1:N
        part1=(1/(w_i(i,1)+1)^2)*tau1_ij(i,j);
        part2=0;
        for l=1:N
            part2=part2+(w_i(i,1)/(w_i(i,1)+1)^2*p_ij(i,l))*(-tau1_ij(i,l)+tau1_ij(j,l));
        end
        part3=0;
        for l=1:N
            part3_1=0;
            for ll=1:N
                part3_1=part3_1+w_i(l,1)/(w_i(l,1)+1)*p_ij(l,ll)*(-tau1_ij(i,ll)+tau1_ij(j,ll));
            end
            part3=part3+(w_i(i,1)/(w_i(i,1)+1)*p_ij(i,l))*(1/(w_i(l,1)+1)*(-tau1_ij(i,l)+tau1_ij(j,l))+part3_1);
        end
        deno=deno+w_ij(i,j)/(w_i(i,1)*w_i(j,1))*(w_i(i,1)+1)*(part1+part2+part3);
    end
end

r_BD_accu=nume/deno;


nume=0;
for i=1:N
    for j=1:N
        nume=nume+w_ij(i,j)/(w_i(i,1)*w_i(j,1))*tau1_ij(i,j);
    end
end

deno=0;
for i=1:N
    for j=1:N
        for l=1:N
            deno=deno+w_ij(i,j)*w_ij(i,l)/(w_i(i,1)^2*w_i(j,1))*(tau1_ij(j,l)-tau1_ij(i,l));
        end
    end
end

bc_BD=nume/deno;


nume=0;
for i=1:N
    for j=1:N
        nume=nume+p_ij(i,j)*tau1_ij(i,j);
    end
end

deno=0;
for i=1:N
    for j=1:N
        for l=1:N
            deno=deno+p_ij(j,i)*p_ij(i,l)*(tau1_ij(j,l)-tau1_ij(i,l));
        end
    end
end

bc_BD_accu=nume/deno;

end