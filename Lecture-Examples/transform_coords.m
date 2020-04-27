function [B1, B2, B3] = transform_coords(A1, A2, A3, Ma2b)

B1 = Ma2b(1,1)*A1 + Ma2b(1,2)*A2 + Ma2b(1,3)*A3 + Ma2b(1,4);
B2 = Ma2b(2,1)*A1 + Ma2b(2,2)*A2 + Ma2b(2,3)*A3 + Ma2b(2,4);
B3 = Ma2b(3,1)*A1 + Ma2b(3,2)*A2 + Ma2b(3,3)*A3 + Ma2b(3,4);

end