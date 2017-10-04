%	Simple spin echo example, with M = 0.81 at 1st spin echo.
%
%	Display state at each step (after trim)
%
T = 1;
T2 = 9.4877;			% exp(-2T/T2)=0.81
T1 = 19.4932;			% exp(-T/T1)=0.95
Q = [0;0;1]
Q = epg_rf(Q,pi/2,pi/2)
Q = epg_grelax(Q,T1,T2,T);
Q = epg_trim(Q,.001)
Q = epg_rf(Q,pi,0);
Q = epg_trim(Q,.001)
Q = epg_grelax(Q,T1,T2,T);
Q = epg_trim(Q,.001)
Q = epg_grelax(Q,T1,T2,T);
Q = epg_trim(Q,.001)
Q = epg_rf(Q,pi,0);
Q = epg_trim(Q,.001)
Q = epg_grelax(Q,T1,T2,T);
Q = epg_trim(Q,.001)

