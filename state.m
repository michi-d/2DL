
classdef state
    
    properties
        
        origin;
        now;
        sequence;
        index_pos_lin;
        index_pos_sub;
        size0;
        cost;
        ID;
        
    end
        
    methods
        
        function obj = state(start_matrix)
            obj.origin = start_matrix;
            obj.now    = start_matrix;
            obj.cost   = 0;
            obj.sequence   = cell(0);

            obj.size0  = size(start_matrix.content);
            obj.index_pos_sub = zeros(obj.size0);
            obj.index_pos_lin  = cell(numel(start_matrix.content), 1);
            
            [is, js] = ind2sub(size(start_matrix.content), 1:numel(start_matrix.content));
            for s = 1:numel(start_matrix.content)
                if start_matrix.content(s) ~= 0
                    obj.index_pos_lin{s} = [is(s), js(s)];
                    obj.index_pos_sub(ind2sub(obj.size0, s)) = s;
                end
            end
            
        end
        
        function new = apply_sequence(obj, seq)
            
            % CODE:
            % 1: row_insert
            % 2: column_insert
            % 3: row_delete
            % 4: column_delete
            % 5: substitution
            % 6: no operations
            %
            % seq : cell(L, 1) : { {CODE, position, element} ; ... } 
            
            new = obj;
            
            for s = 1:length(seq)
                
                switch seq{s}{1}
                    
                    case 1 % row insertion
                        i_lin = sub2ind(new.size0, seq{s}{2}(1), seq{s}{2}(2)); 
                        ij = new.index_pos_lin{i_lin};
                        [M_, co, error] = new.now.row_insert(seq{s}{3}, ij(1), ij(2));
                        if error == 0
                            [new_index_pos_lin, new_index_pos_sub] = index_row_insert(new, ij);
                            new.index_pos_lin = new_index_pos_lin;
                            new.index_pos_sub = new_index_pos_sub;
                            new.sequence = [new.sequence, {seq{s}}];
                            new.now = M_;
                            new.cost = new.cost + co;
                        end
                        
                    case 2 % column insertion
                        i_lin = sub2ind(new.size0, seq{s}{2}(1), seq{s}{2}(2)); 
                        ij = new.index_pos_lin{i_lin};
                        [M_, co, error] = new.now.column_insert(seq{s}{3}, ij(1), ij(2));
                        if error == 0
                            [new_index_pos_lin, new_index_pos_sub] = index_column_insert(new, ij);
                            new.index_pos_lin = new_index_pos_lin;
                            new.index_pos_sub = new_index_pos_sub;
                            new.sequence = [new.sequence, {seq{s}}];
                            new.now = M_;
                            new.cost = new.cost + co;
                        end
                        
                    case 3 % row deletion
                        i_lin = sub2ind(new.size0, seq{s}{2}(1), seq{s}{2}(2)); 
                        ij = new.index_pos_lin{i_lin};
                        [M_, co, error] = new.now.row_delete(ij(1), ij(2));
                        if error == 0
                            [new_index_pos_lin, new_index_pos_sub] = index_row_delete(new, ij);
                            new.index_pos_lin = new_index_pos_lin;
                            new.index_pos_sub = new_index_pos_sub;
                            new.sequence = [new.sequence, {seq{s}}];
                            new.now = M_;
                            new.cost = new.cost + co;
                        end
                        
                    case 4 % column deletion
                        i_lin = sub2ind(new.size0, seq{s}{2}(1), seq{s}{2}(2)); 
                        ij = new.index_pos_lin{i_lin};
                        [M_, co, error] = new.now.column_delete(ij(1), ij(2));
                        if error == 0
                            [new_index_pos_lin, new_index_pos_sub] = index_column_delete(new, ij);
                            new.index_pos_lin = new_index_pos_lin;
                            new.index_pos_sub = new_index_pos_sub;
                            new.sequence = [new.sequence, {seq{s}}];
                            new.now = M_;
                            new.cost = new.cost + co;
                        end
                        
                    case 5 % substitution
                        i_lin = sub2ind(new.size0, seq{s}{2}(1), seq{s}{2}(2)); 
                        ij = new.index_pos_lin{i_lin};
                        [M_, co, error] = new.now.substitute(seq{s}{3}, ij(1), ij(2));
                        if error == 0
                            new.sequence = [new.sequence, {seq{s}}];
                            new.now = M_;
                            new.cost = new.cost + co;
                        end
                        
                    case 6 % no operation
                        [M_, co, error] = new.now.clone();
                        if error == 0
                            new.sequence = [new.sequence, {seq{s}}];
                            new.now = M_;
                            new.cost = new.cost + co;
                        end
 
                end
                
                if error ~= 0
                   fprintf(1, 'Forbidden operation in sequence: STOPPED.\n');
                   new = [];
                   break
                end
               
            end
        end
        
        
        function [new_index_pos_lin new_index_pos_sub] = index_row_insert(obj, ij)
          new_index_pos_lin = obj.index_pos_lin;
          new_index_pos_sub = obj.index_pos_sub;
          
          if new_index_pos_sub(ij(1),end) ~= 0
              new_index_pos_sub(ij(1),end+1) = 0;
          end
          new_index_pos_sub(ij(1), ij(2)+1:end) =  new_index_pos_sub(ij(1), ij(2):end-1);
          new_index_pos_sub(ij(1), ij(2)) = 0;
          
          lins = new_index_pos_sub(ij(1), ij(2):end);
          lins = lins(lins ~=0);
          
          for s = 1:length(lins)
              new_index_pos_lin{lins(s)} = new_index_pos_lin{lins(s)} + [0 1];
          end
        end
        
        
        function [new_index_pos_lin new_index_pos_sub] = index_column_insert(obj, ij)
          new_index_pos_lin = obj.index_pos_lin;
          new_index_pos_sub = obj.index_pos_sub;
          
          if new_index_pos_sub(end, ij(2)) ~= 0
              new_index_pos_sub(end+1, ij(2)) = 0;
          end
          new_index_pos_sub(ij(1)+1:end, ij(2)) = new_index_pos_sub(ij(1):end-1, ij(2));
          new_index_pos_sub(ij(1), ij(2)) = 0;
          
          lins = new_index_pos_sub(ij(1):end, ij(2));
          lins = lins(lins ~=0);
          
          for s = 1:length(lins)
              new_index_pos_lin{lins(s)} = new_index_pos_lin{lins(s)} + [1 0];
          end
        end
        
        function [new_index_pos_lin new_index_pos_sub] = index_row_delete(obj, ij)
          new_index_pos_lin = obj.index_pos_lin;
          new_index_pos_sub = obj.index_pos_sub;
          i_lin = new_index_pos_sub(ij(1), ij(2));
          
          new_index_pos_sub(ij(1), ij(2):end-1) = new_index_pos_sub(ij(1), ij(2)+1:end); 
          new_index_pos_sub(ij(1), end) = 0;
          
          lins = new_index_pos_sub(ij(1), ij(2):end);
          lins = lins(lins ~=0);
          
          for s = 1:length(lins)
              new_index_pos_lin{lins(s)} = new_index_pos_lin{lins(s)} - [0 1];
          end
          new_index_pos_lin{i_lin} = [];
        end
        
        function [new_index_pos_lin new_index_pos_sub] = index_column_delete(obj, ij)
          new_index_pos_lin = obj.index_pos_lin;
          new_index_pos_sub = obj.index_pos_sub;
          i_lin = new_index_pos_sub(ij(1), ij(2));
          
          new_index_pos_sub(ij(1):end-1, ij(2)) = new_index_pos_sub(ij(1)+1:end, ij(2)); 
          new_index_pos_sub(end, ij(2)) = 0;
          
          lins = new_index_pos_sub(ij(1):end, ij(2));
          lins = lins(lins ~=0);
          
          for s = 1:length(lins)
              new_index_pos_lin{lins(s)} = new_index_pos_lin{lins(s)} - [1 0];
          end
          new_index_pos_lin{i_lin} = [];
        end
        
        
        function show_positions(obj)
             [N, M] = size(obj.index_pos_sub);
             for i = 1:N
                for j = 1:M
                    if obj.index_pos_sub(i,j) ~= 0
                      [i_, j_] = ind2sub(obj.size0, obj.index_pos_sub(i,j));
                      str = ['(', num2str(i_), ', ', num2str(j_), ')'];
                      str = [str, repmat(' ', 1, 9 - length(str))];
                    else
                      if obj.now.content(i,j) == 0
                        str = '    -    ';
                      else 
                        str = '    X    ';
                      end
                    end
                    fprintf(1, str);
                end
                fprintf(1, '\n');
             end
        end
        
        function new = step_n(obj, n)
            new = state(obj.origin);
            new = new.apply_sequence(obj.sequence(1:n));
        end
        
        function animate_life(obj, dt)
            for s = 0:length(obj.sequence)
              clc
              t_ = obj.step_n(s);
              t_.now.print_content
              pause(dt)
            end
        end
            

        
      
    end
        
        
        
    
end