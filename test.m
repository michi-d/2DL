
% define metric space
m = metric_space(0, 1, 1, 1);

% define original matrix
S = 'ABCC\nBBB\nACB';
a = matrix([], m);
a = a.content2matrix(S);
  
% print original matrix
a.print_content;

% generate a state which start with a as original matrix
s = state(a);

% define transformation sequence
% CODE:
% 1: row_insert
% 2: column_insert
% 3: row_delete
% 4: column_delete
% 5: substitution
% 6: no operations
% structure: {CODE, position, element}
% element is given in ASCII code so 80 means 'P' for example
seq = { {1, [1,1], 80 }, {3, [1,1]} , {3, [1,2]}, {4, [1,4]}, {4, [2,2]} };

% apply transformation sequence
t = s.apply_sequence(seq);

% print current state after all transformations
t.now.print_content;

% show current positions of elements which were in the original matrix
t.show_positions;

% animate life of the matrix with a time constant of 1 sec
%t.animate_life(1)

b = boundary(a, 2);
b = b.extend_boundary([1,1], 66);
b = b.extend_boundary([1,2], 66);
b = b.extend_boundary([2,1], 66);
b = b.extend_boundary([1,3], 66);
b = b.extend_boundary([2,2], 66);
b = b.extend_boundary([2,3], 66);
%b = b.extend_boundary([1,4], 66);
%b = b.extend_boundary([3,1], 66);
b.draw_boundary;

b = b.crop_boundary;


costs = [];
for i = 1:length(b.children)
    costs = [costs, b.children{i}.cost];
end
[min_cost i] = min(costs);
minima = find(costs == min_cost);

