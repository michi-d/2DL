
classdef boundary
    
    properties
        
        parent;  
        elements;
        children;
        step_max;
        
    end
    
    methods
        
        function obj = boundary(start_matrix, step_max)
            
            obj.parent = start_matrix;
            
            obj.step_max = step_max;
            obj.children = {};
            
            obj.elements = zeros(0,0); % [i_1, j_1;
                                       %  i_2, j_2;
                                       %  ...
                                       %  i_n, j_n]
            
        end
        
        
        function [new_boundary] = extend_boundary(obj, position, bij)
           
            % extends the boundary to a new position [i,j]
            
            new_boundary = boundary(obj.parent, obj.step_max);
            new_boundary.elements = [obj.elements; [position(1), position(2)]];
            
            % first check if new position is valid
            if ~(isempty(obj.elements))
                    columns = obj.elements(:,2);
                    rows    = obj.elements(:,1);
                if sum((position(1) == rows) & (position(2) == columns))
                    error('position is not valid');
                end
                if ~(  ((position(1)-1 == rows) & (position(2) == columns)) | ...
                    ((position(1) == rows) & (position(2)-1 == columns)) )
                
                    error('position is not valid');
               end
            else 
                if ~(position(1) == 1 && position(2) == 1)
                    error('position is not valid');
                else
                    columns = [];
                    rows    = [];
                end
            end
            
            % For every children state look in the environment of staircase
            % walks of length step_max in A' if there is an element equal to b_ij
            % that could be put in position via a staircase walk of
            % deletions.
            %
            % - Create a new children for every possibility.
            % 
            % - Furthermore, create two new children by:
            %     1.) substituting the current element.
            %  or 2.) insert b_ij.
            %
            % (maybe one should consider also substituting elements nearby which
            %  are close to b_ij in terms of substitution).
           
           
            if isempty(obj.elements)
                obj.children{1} = state(obj.parent);
            end
            
            for c = 1:length(obj.children);
                
            % deletion staircase walks   
            
            [square_x, square_y] = meshgrid(1:obj.step_max + 1, 1:obj.step_max + 1);
            square_x = square_x - 1;
            square_y = square_y - 1;
            walk_distance = square_x(:)  + square_y(:);
            possible_x = square_x(walk_distance <= obj.step_max);
            possible_y = square_y(walk_distance <= obj.step_max);
            possible_x = possible_x(2:end) + position(1);
            possible_y = possible_y(2:end) + position(2);
            possible_x = possible_x(possible_x <= size(obj.parent.content, 1));
            possible_y = possible_y(possible_x <= size(obj.parent.content, 1));
            possible_x = possible_x(possible_y <= size(obj.parent.content, 2));
            possible_y = possible_y(possible_y <= size(obj.parent.content, 2));
                
            for k = 1:length(possible_x)
                
                if obj.children{c}.now.content(possible_x(k), possible_y(k)) == bij
                    % possible element found
                    %fprintf('%i   %i   \n', possible_x(k), possible_y(k));
                    
                    % find all possible staircase walks...
                    diff = [possible_x(k), possible_y(k)] - position;
                    
                    % distribute diff(1) vertical walking steps to
                    % sum(diff) steps to be done
                    step_pos = nchoosek(1:sum(diff), diff(1)); % choose when to make vertical step
                    step_pos = step_pos + 1;
                    if diff(1) == 0
                        step_pos = zeros(1, sum(diff));
                    end
                  
                    
                    possible_walks = zeros(2, sum(diff) + 1, size(step_pos,1));
                    possible_walks(1,1,:) = possible_x(k);
                    possible_walks(2,1,:) = possible_y(k);
                    for i = 1:size(possible_walks,3)
                        for n = 2:sum(diff) + 1
                            if sum(step_pos(i,:) == n) % if vertical step occurs in staircase walk at step n
                                possible_walks(1, n, i) = possible_walks(1, n-1, i) - 1; 
                                possible_walks(2, n, i) = possible_walks(2, n-1, i); 
                            else  % if vertical step occurs in staircase walk at step n
                                possible_walks(1, n, i) = possible_walks(1, n-1, i); 
                                possible_walks(2, n, i) = possible_walks(2, n-1, i) - 1; 
                            end
                        end
                    end
                    
                    % check if staircase walks contains only A-elements 
                    
                    delete_walk = [];
                    for i = 1:size(possible_walks, 3)
                        for j = 2:size(possible_walks, 2) % 2, because last element (which gets shifted) does not have to be an A-element necessarily
                            
                           % if any element on the path is not an A-element
                           % erase the path
                           if obj.children{c}.index_pos_sub(possible_walks(1,j,i), possible_walks(2,j,i)) == 0
                               delete_walk = [delete_walk, i];
                               %possible_walks(:,:,i) = [];
                               fprintf(1, 'possible path deleted because there are B-elements on it');
                           end
                                 
                        end
                    end
                    possible_walks(:,:,delete_walk) = [];

                    
                    % walk every remaining possible path and create a new
                    % children state for it
                    
                    for i = 1:size(possible_walks, 3)
                        seq = cell(0,0);
                        for j = 2:size(possible_walks, 2)
                            diff = possible_walks(:,j-1,i) - possible_walks(:,j,i);
                            
                            if diff == [0; 1]  % horizontal step
                                orig_position_lin = obj.children{c}.index_pos_sub(possible_walks(1,j,i), possible_walks(2,j,i));
                                [ii, jj] = ind2sub(obj.children{c}.size0, orig_position_lin);
                                orig_position = [ii, jj];
                                seq{j-1} = {3, orig_position};
                            elseif diff == [1; 0]  % vertical step   
                                orig_position_lin = obj.children{c}.index_pos_sub(possible_walks(1,j,i), possible_walks(2,j,i));
                                [ii, jj] = ind2sub(obj.children{c}.size0, orig_position_lin);
                                orig_position = [ii, jj];
                                seq{j-1} = {4, orig_position};
                            end
                            
                        end
                        % walk staircase walk and create new children state
                        % for new boundary
                        temp = obj.children{c}.apply_sequence(seq);
                        if ~isempty(temp)
                            new_boundary.children{end + 1} = temp;
                        end
                    end
                     
                end
             
            end
            
            
            % substituting current element
            
            if obj.children{c}.index_pos_sub(position(1), position(2)) == 0
                % errorcheck if element is an A-element
                fprintf(1, 'substitution not possible because element is a B-element');
            else
                orig_position_lin = obj.children{c}.index_pos_sub(position(1), position(2));
                [ii, jj] = ind2sub(obj.children{c}.size0, orig_position_lin);
                orig_position = [ii, jj];
                
                seq = { {5, orig_position, bij} };
                new_boundary.children{end + 1} = obj.children{c}.apply_sequence(seq);
            end
            
            % inserting bij at current element
            
            if obj.children{c}.index_pos_sub(position(1), position(2)) == 0
                % errorcheck if element is an A-element
                fprintf(1, 'insertion not possible because element is a B-element');
            else
                orig_position_lin = obj.children{c}.index_pos_sub(position(1), position(2));
                [ii, jj] = ind2sub(obj.children{c}.size0, orig_position_lin);
                orig_position = [ii, jj];
                
                % row-insertion
                seq = { {1, orig_position, bij} };
                new_boundary.children{end + 1} = obj.children{c}.apply_sequence(seq);
                
                % column-insertion
                seq = { {2, orig_position, bij} };
                new_boundary.children{end + 1} = obj.children{c}.apply_sequence(seq);
            end
            
                      
            end
            
            
            % delete empty children
            
            delete_children = [];
            for c = 1:length(new_boundary.children)
                if isempty(new_boundary.children{c})
                    delete_children = [delete_children, c];
                end
            end
            new_boundary.children(delete_children) = [];
            
            
        end
        
        
        function [new_boundary] = crop_boundary(obj)
            
            % deletes all elements which are outside the boundary elements
            
            new_boundary = obj;
            
            ne = max(obj.elements(:,1));
            me = max(obj.elements(:,2));
            elements_boundary = zeros(ne, me);
            for i = 1:size(obj.elements, 1)
                elements_boundary(obj.elements(i,1), obj.elements(i,2)) = 1;
            end
            
            for c = 1:length(obj.children)
                
                proz = c/length(obj.children) * 100;
                proz_before = (c-1)/length(obj.children) * 100;
                
                if ceil(proz_before/10) < ceil(proz/10)
                    fprintf(1, '%f %% \n', proz);
                end
                
                mask = zeros(size(obj.children{c}.now.content));
                mask(1:ne, 1:me) = elements_boundary;
                
                existing_elements = obj.children{c}.now.content ~= 0;
                outside_elements = existing_elements .* ~mask;
                
                [row, col] = find(outside_elements);
                dist = row+col;
                [a, order] = sort(dist, 'descend'); 
                
                seq = cell(1, length(order));
                for n = 1:length(order)
                    orig_position_lin = new_boundary.children{c}.index_pos_sub(row(order(n)), col(order(n)));
                    [ii, jj] = ind2sub(new_boundary.children{c}.size0, orig_position_lin);
                    orig_position = [ii, jj];
                    
                    seq{n} = {3, orig_position};
                end
                
                new_boundary.children{c} = new_boundary.children{c}.apply_sequence(seq);
            end
            
            
        end
        
        
        
        function draw_boundary(obj)
            
            obj.parent.draw_matrix;
            
            columns = obj.elements(:,2);
            rows    = obj.elements(:,1);
            
            max_columns = max(columns);
            
            points = [0,0];
            for n = 1:max_columns
                rows_nth_column = max(rows(columns == n));
                
                points = [points; [n-1, rows_nth_column] ];
                points = [points; [n, rows_nth_column] ];
            end
            points = [points; [n, 0] ];
            points = [points; [0, 0] ]; 
            
            points = points + 0.5;
            
            b_box = fill(points(:,1), points(:,2),  'r', 'LineWidth', 2, 'edgecolor', 'b', 'FaceAlpha', 0);
            obj.parent.draw_matrix;
        end
        

        
    end
    
    
end