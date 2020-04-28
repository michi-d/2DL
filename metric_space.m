
classdef metric_space
    
    properties        
        type;
        cost;
        cost_add;
        cost_del;
        name;
    end
    
    methods
        
        function obj = metric_space(type, cost, cost_add, cost_del)
            obj.type = type;
            obj.cost = cost;
            obj.cost_add = cost_add;
            obj.cost_del = cost_del;
            if type == 0
                obj.name = 'letters';
            end
            if type == 1
                obj.name = 'numbers';
            end
        end
        
        function d = distance(obj, a,b)
            if obj.type == 0
                d = obj.cost * double(a ~= b);
            elseif obj.type == 1
                d = obj.cost * double(abs(a-b));
            end
        end
        
        function [out] = print_element(obj, a)
            if obj.type == 0
                out = char(a);
                if a == 0
                    out = ' ';
                end
            elseif obj.type == 1
                out = num2str(a);
            end
        end
        
        function [out] = get_element(obj, a)
            if obj.type == 0;
                out = double(a);
            end
            if obj.type == 1;
                out = str2double(a);
            end
        end
        
        
    end
    
end