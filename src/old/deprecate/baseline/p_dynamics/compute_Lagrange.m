function L_e = compute_Lagrange(K,N,t,t_nodes)

L_e = ones(N+1,K+1); %up to order N

for k=0:N
    for m = 0:N
        if m~=k
            L_e(k+1,:) = L_e(k+1,:).*((t-t_nodes(m+1))/(t_nodes(k+1)-t_nodes(m+1)));
        end
    end
    
end


end
