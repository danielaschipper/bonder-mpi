
bonder : Analize.cpp fill.cpp iface.cpp LLi.cpp main.cpp readwfn.cpp stdafx.cpp output.cpp minion.cpp line.cpp intersect.cpp
	mpic++ -Wall -o mpibond Analize.cpp fill.cpp iface.cpp LLi.cpp main.cpp readwfn.cpp stdafx.cpp output.cpp minion.cpp line.cpp intersect.cpp -I .   -std=c++11 -pthread -static-libgcc -Ofast

