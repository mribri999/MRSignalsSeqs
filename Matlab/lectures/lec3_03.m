% Lecture 3 - Example 3
%
% Compact animation by splitting [A,B] propagation
%
[A,B] = relax(5,2,1.5);      % Relaxation over 5sec with T1=2sec and T2=1.5sec.
M=abanim([1;0;0],40,A,B);    % Note we can't "split" a 360-degree (or more) rotation, due to ambiguity.


