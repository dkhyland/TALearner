function [] = PMC2NFA(ProdMC, L, accept_states, file_name)
    
    % Convert Markov chain matrix into unweighted digraph adjacency matrix
    cutoff = 1e-3;
    dim = length(ProdMC);
    
    % eliminate out-going edges from states with no in-coming edges
    for i=1:dim
        %disp(sum(ProdMC(:,i)))
        if all(ProdMC(:,i) < cutoff)
            disp(i)
            ProdMC(i,:) = zeros(1,dim);
        end
    end
    
    ProdMC = ProdMC > cutoff;
    
    Alphabet = 'abcdefghijklmnopqrstuvwxyz'; % only up to 26 labels for now
    
    fid = fopen( file_name, 'wt' );
    fprintf( fid, '%s\n', char(string(dim)) ); 
    fprintf( fid, '%s\n', Alphabet(unique(L)) ); 
    fprintf( fid, '%s ', char(string(length(accept_states))) ); 
    for i=1:length(accept_states)
        fprintf( fid, '%s ', char(string(accept_states(i)-1)) ); 
    end
    fprintf(fid, '\n');
    fprintf(fid, '0\n');
    
    maxL = max(L);
    for i =1:accept_states(1)-1
        temp = zeros(1,maxL);
        for j=1:dim
            if ProdMC(i,j) == 1
                fprintf( fid, '%s %s %s\n', char(string(i-1)), Alphabet(L(j)), char(string(j-1)) ); 
                temp(L(j)) = 1;
            end
        end
        % If you want to add self-loops for all impossible next-labels:
%         for a=1:maxL
%             if temp(a) == 0
%                 fprintf( fid, '%s %s %s\n', char(string(i-1)), Alphabet(a), char(string(i-1)) ); 
%             end
%         end
    end
    for i = accept_states
        for a =1:maxL
            fprintf( fid, '%s %s %s\n', char(string(i-1)), Alphabet(a), char(string(i-1)) ); 
        end
    end
    
    fclose(fid);
    
end