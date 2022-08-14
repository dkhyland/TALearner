function P_in = reduceMC(P_in)
    % pre-process
    i = 2;
    while i <= size(P_in, 2)
        if P_in(:,i) == 0
            P_in(:,i) = [];
            P_in(i,:) = [];
        else
            i=i+1;
        end
    end
    cutoff=1e-3;
    i = 1;
    while i < size(P_in, 2)
        j = i+1;
        temp = zeros(size(P_in,2), 1);
        while j <= size(P_in, 2)
            if ( (P_in(:,i) > cutoff) ==  (P_in(:,j) > cutoff) )
                temp(j) = 1;
            end
            j = j+1;
        end
        for j=i+1:size(P_in, 2)
            if temp(j)
                P_in(:,i) = P_in(:,i) + P_in(:,j);
            end
        end
        for j=i+1:size(P_in, 2)
            if temp(j)
                P_in(i,:) = P_in(i,:) + P_in(j,:);
            end
        end
        for j=size(P_in, 2):-1:i+1
            if temp(j)
                P_in(j,:) = [];
                P_in(:,j) = [];
            end
        end
        P_in(i,:) = P_in(i,:)./sum(P_in(i,:));
        i = i +1;
    end
end