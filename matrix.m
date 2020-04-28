classdef matrix
    
    properties
        
        N = 0;            % number of rows
        M = 0;            % number of columns
        N_n = zeros(0,1)  % number of elements in row n
        M_n = zeros(0,1)  % number of elements in column n
        
        content = zeros(0,0);   % content of matrix
                                % 0 means non-existence of an element
                                
        metric;
    end
    
    methods
        
        function obj = matrix(content, metric)
           % class constructor
           if(nargin > 0)
             obj.content = content;
             obj.metric  = metric;
             obj.N = size(content,1);
             obj.M = size(content,2);
             [obj.N_n, obj.M_n] = count_elements(obj);
           end
        end
           
        function [N_n M_n] = count_elements(obj)
            % count elements per row and column
            A = obj.content;
            N_n = zeros(size(A,1),1);
            M_n = zeros(size(A,2),1);
            for i = 1:length(N_n)
                N_n(i) = sum(A(i,:) ~= 0);
            end
            for i = 1:length(M_n)
                M_n(i) = sum(A(:,i) ~= 0);
            end
        end
        
        function print_content(obj)
            for i = 1:obj.N
                for j = 1:obj.M
                    fprintf(1, obj.metric.print_element(obj.content(i,j)));
                end
                fprintf(1, '\n');
            end
        end
        
        function [new] = content2matrix(obj, S)
            tmp = [];
            cont = [];
            row = 1;
            i = 1;
            while i <= length(S)
                cache = S(i);
                if cache == '\'
                    cont(row, 1:length(tmp)) = tmp;
                    row = row + 1;
                    i = i + 1;
                    if S(i) == 'n'
                        i = i + 1;
                    end
                    tmp = [];
                else
                    tmp = [tmp, obj.metric.get_element(S(i))];
                    i = i + 1;
                end
            end
            cont(row, 1:length(tmp)) = tmp;
            new = matrix(cont, obj.metric);
            
        end
        
        function [bool] = check_definition(obj)
            % checks if a logical mask of the matrix containing the position of elements
            % is conform with the definition

            mask = obj.content ~= 0;
            bool = 1;
            
            top_neighbours = find( mask(2:end, :) - mask(1:end-1, :) == 1, 1); 
            if ~isempty(top_neighbours)
                bool = 0;
            end
            
            right_neighbours = find( mask(:, 2:end) - mask(:, 1:end-1) == 1, 1); 
            if ~isempty(right_neighbours)
                bool = 0;
            end
            
        end
        
        
        function [matrix_out cost error] = row_insert(obj, a, i, j)
            % row-add element a at position (i,j) within the content matrix
            matrix_out = obj;
            error = 0;
            
            % check if element exists
            if matrix_out.content(i, j) == 0
                error = 1;
            end
            
            if error == 0
              matrix_out.content(i,j+1:matrix_out.N_n(i)+1) = matrix_out.content(i,j:matrix_out.N_n(i));
              matrix_out.content(i,j) = a;
              [matrix_out.N, matrix_out.M] = size(matrix_out.content);
              [matrix_out.N_n, matrix_out.M_n] = count_elements(matrix_out);
              
              cost = matrix_out.metric.cost_add;
            end
            
            % check if matrix is conform with the definition
            if ~matrix_out.check_definition
                error = 2;
                fprintf(1, 'Operation is not conform with the definition of a matrix!\n');
            end
            
            if error == 1
                fprintf(1, 'Row-insertion at non-existent index not possible!\n');
            end
            

            if error ~= 0
                cost = [];
            end
        end
        
        
        function [matrix_out cost error] = column_insert(obj, a, i, j)
            % column-add element a at position (i,j) within the content matrix
            matrix_out = obj;
            error = 0;
            
            % check if element exists
            if matrix_out.content(i, j) == 0
                error = 1;
            end
            
            if error == 0
              matrix_out.content(i+1:matrix_out.M_n(j)+1, j) = matrix_out.content(i:matrix_out.M_n(j), j);
              matrix_out.content(i,j) = a;
              [matrix_out.N, matrix_out.M] = size(matrix_out.content);
              [matrix_out.N_n, matrix_out.M_n] = count_elements(matrix_out);
              
              cost = matrix_out.metric.cost_add;
            end
            
            % check if matrix is conform with the definition
            if ~matrix_out.check_definition
                error = 2;
                fprintf(1, 'Operation is not conform with the definition of a matrix!\n');
            end
            
            if error == 1
                fprintf(1, 'Column-insertion at non-existent index not possible!\n');
            end
            
            if error ~= 0
                cost = [];
            end
        end
        
        
        function [matrix_out cost error] = row_delete(obj, i, j)
            % row-delete element at position (i,j) within the content matrix
            matrix_out = obj;
            error = 0;
            
            % check if element exists
            if matrix_out.content(i, j) == 0
                error = 1;
            end
            
            
            % check if row-deletion generates holes in column
%             neighbours = [i-1, i+1];
%             bool = 1;
%             if i == 1
%                 neighbours = i+1;
%             elseif i == matrix_out.N
%                 neighbours = i-1;
%                 bool = 0;
%             end
%             for n = neighbours
%                 bool = (matrix_out.N_n(i) - 1 < matrix_out.N_n(n)) && bool;
%             end
%             if bool
%                 error = 2;
%             end
%             
            if error == 0
              matrix_out.content(i, j:matrix_out.M-1) = matrix_out.content(i, j+1:matrix_out.M);
              matrix_out.content(i, matrix_out.M) = 0;
              [matrix_out.N, matrix_out.M] = size(matrix_out.content);
              [matrix_out.N_n, matrix_out.M_n] = count_elements(matrix_out);
              
              cost = matrix_out.metric.cost_del;
            end
            
            % check if matrix is conform with the definition
            if ~matrix_out.check_definition
                error = 2;
                fprintf(1, 'Operation is not conform with the definition of a matrix!\n');
            end
            
            if error == 1
                fprintf(1, 'Row-deletion at non-existent index not possible!\n');
            end
            
%             if error == 2
%                 fprintf(1, 'Row-deletion not possible: no holes in columns allowed!\n');
%             end
            
            if error ~= 0
                cost = [];
            end
        end
        
        
        
        function [matrix_out cost error] = column_delete(obj, i, j)
            % row-delete element at position (i,j) within the content matrix
            matrix_out = obj;
            error = 0;
            
            % check if element exists
            if matrix_out.content(i, j) == 0
                error = 1;
            end
            
            % check if row-deletion generates holes in column
%             neighbours = [j-1, j+1];
%             bool = 1;
%             if j == 1
%                 neighbours = j+1;
%             elseif j == matrix_out.M
%                 neighbours = j-1;
%                 bool = 0;
%             end
%             for n = neighbours
%                 bool = (matrix_out.M_n(j) - 1 < matrix_out.M_n(n)) && bool;
%             end
%             if bool
%                 error = 2;
%             end

            
            if error == 0
              matrix_out.content(i:matrix_out.N-1, j) = matrix_out.content(i+1:matrix_out.N, j);
              matrix_out.content(matrix_out.N, j) = 0;
              [matrix_out.N, matrix_out.M] = size(matrix_out.content);
              [matrix_out.N_n, matrix_out.M_n] = count_elements(matrix_out);
              
              cost = matrix_out.metric.cost_del;
            end
            
            % check if matrix is conform with the definition
            if ~matrix_out.check_definition
                error = 2;
                fprintf(1, 'Operation is not conform with the definition of a matrix!\n');
            end
            
            if error == 1
                fprintf(1, 'Column-deletion at non-existent index not possible!\n');
            end
            
%             if error == 2
%                 fprintf(1, 'Column-deletion not possible: no holes in columns allowed!\n');
%             end
            
            if error ~= 0
                cost = [];
            end
                
        end
        
        
        function [matrix_out cost error] = substitute(obj, a, i, j)
            % substitute element position (i,j) within element a
            matrix_out = obj;
            error = 0;
            
            % check if element exists
            if matrix_out.content(i, j) == 0
                error = 1;
            end
            
            if error == 0
                cost = matrix_out.metric.distance(matrix_out.content(i,j), a);
                matrix_out.content(i,j) = a;  
            end
            
            if error == 1
                fprintf(1, 'Substitution at non-existent index not possible!\n');
            end
            
            if error ~= 0
                cost = [];
            end
            
        end
        
        
         function [matrix_out cost error] = clone(obj)
            % no operation: output identical matrix
            error = 0;
            cost = 0;
            matrix_out = obj;
         end
         
         
         function draw_matrix(obj)
             
%              obj.content = content;
%              obj.metric  = metric;
%              obj.N = size(content,1);
%              obj.M = size(content,2);
             
             
            %figure
            %axes
            set(gca, 'XTick', 0.5:obj.M+0.5);
            set(gca, 'XTickLabel', 1:obj.M);
            set(gca, 'YTick', 0.5:obj.N+0.5);
            set(gca, 'YTickLabel', 1:obj.N);
            set(gca, 'XLim', [0,obj.M+1]);
            set(gca, 'YLim', [0,obj.N+1]);
            set(gca, 'YDir', 'reverse');
            grid on
            
            for i = 1:obj.N
                for j = 1:obj.M
                    string = obj.metric.print_element(obj.content(i,j));
                    text(j,i,string);
                end
            end
                        
            
        end
        
               
    end
    

end