clc
clear
//~~~Given Data~~~`
g=9.81; //m/s^(-2)
global('g')
eps_max=1e-6;
betav=2; //exponent of discharge.
pipen=8; //Number of pipes.
nmaxpl=4;//Maximum number of pipes connected to a particular loop.
junn=5; //Number of junctions.
ploop=1; //Number of pseudo loop.
iloop=2; //Number of interior loop.
loopn=ploop+iloop; //Total number of loops.
Qi=[0.319 0.134 0.062 0.019 0.043 0.185 0.035 0.072]; //initial discharge
head_diff=[20]; //m , in pseudo loop
loop_con=[4 -1 -2 -3 -4
          4  2  8 -7 -6
          3  3  5 -8  0];
//In 'loop_con' marix, 1st column represents the number of pipes connected to that loop. Other shows the pipes connected to that. + sign means flow in clockwise direction.          
pipe_con=zeros(pipen,loopn);
//'pipe_con' should be a 8*3 matrix.

for l=1:loopn
//'l' stands for loop numbering.(i.e, l=1:3).
        for i=1:loop_con(l,1)
//"loop_con(l,1)" represents the number of pipes connected to l-th loop.
            pipe_con(abs(loop_con(l,i+1)),l)=l;
//"for i=1:loop_con(l,1), abs(loop_con(l,i+1))" covers all the 'pipe nos' connected to lth loop. We can see, each row of "loop_con"  matrix is dedicated to a particular pipe no. It shows, in which loop a particular pipe is connected. Output: pipe_con=
//   1.   0.   0.
//   1.   2.   0.
//   1.   0.   3.
//   1.   0.   0.
//   0.   0.   3.
//   0.   2.   0.
//   0.   2.   0.
//   0.   2.   3.
        end
end
Kv=[100 500 200 100 400 300 400 300];
//'Kv' means "Ki(cap)", which considers the total loss thing.

//~~~~~~Problem dependent parameters~~~~~~~~
Qv=zeros(loopn,nmaxpl);
//'Qv' is a (3*4) matrix.

//~~~~~General Identification Matrix~~~~
for l=1:loopn
    for i=1:loop_con(l,1)
        Qv(l,i)=sign(loop_con(l,i+1))*Qi(abs(loop_con(l,i+1)));
//'loop_con(l,i+1)' because channel numbering starts from 2nd column. For each 'pipe loop(l)', i runs for 'loop_con(l,1)' times. Because, 'loop_con(l,1)' is the number of pipes connected to l-th loop. Here, Qv considers the sign of discharge also. Because, in loop_con matrix we considered the pipe numberings with appropriate sign (Clockwise + ,Anticlockwise -)
    end
end

count=0;
rmse=1;

//~~~~~~~~Space Loop~~~~~~~~~~~~

while rmse>eps_max
    rmse=0;
    delQ=zeros(loopn,1);
//We are considering a particular delQ value for all pipes in each loop.

for l=1:loopn
    nr=0; 
    dr=0;
//numerator and denominator initialized as zero values.
    if (l<=ploop)then
//Because, we have assigned the pseudo loop as loop no 1. As ploop=1, therefore, l<=1 implies the pseudo loop , where a particular head difference is present between two reservoirs.
        nr=nr+head_diff(l);
    end
    for i=1:loop_con(l,1)
//"loop_con(l,1)" represents number of pipes connected to l-th loop.
        nr=nr+Kv(abs(loop_con(l,i+1)))*Qv(l,i)*abs(Qv(l,i))^(betav-1);
// In "Kv(abs(loop_con(l,i+1)))" i+1 is taken because,numbering of connected pipes are started from 2nd column of 'loop_con' matrix. 
        dr=dr+betav*Kv(abs(loop_con(l,i+1)))*abs(Qv(l,i))^(betav-1);
//See equation (412) for "nr" and "dr" calculation.
    end
    delQ(l)=-(nr/dr);
end
//From output, delQ= [0.0000010; 0.0000007; 0.0000009]
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for l=1:loopn 
//Considers loops
    for j=1:loop_con(l,1) 
//Considers all pipes connected to that loop.
        for k=1:loopn 
//Related to 'pipe_con' matrix. Considers a pipe connected to how many loops.
            if (pipe_con(abs(loop_con(l,j+1)),k)<>0)then
                if (pipe_con(abs(loop_con(l,j+1)),k)==l)then //this condition is important.
                    Qv(l,j)=Qv(l,j)+delQ(pipe_con(abs(loop_con(l,j+1)),k));
                else
                    Qv(l,j)=Qv(l,j)-delQ(pipe_con(abs(loop_con(l,j+1)),k));
                end
            end
        end
    end
end
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Explanation of the above loop: 
//when l=1,j=1,k=1,   pipe_con(abs(loop_con(1,2)),1)=pipe_con(1,1)=1;    Qv(1,1)=Qv(1,1)+delQ(pipe_con(abs(loop_con(1,2)),1))=Qv(1,1)+delQ(1);
//when l=1,j=1,k=2,   pipe_con(abs(loop_con(1,2)),2)=pipe_con(1,2)=0; -->stop
//when l=1,j=1,k=3,   pipe_con(abs(loop_con(1,2)),3)=pipe_con(1,3)=0; -->stop
//when l=1,j=2,k=1,   pipe_con(abs(loop_con(1,3)),1)=pipe_con(2,1)=1;    Qv(1,2)=Qv(1,2)+delQ(pipe_con(abs(loop_con(1,3)),1))=Qv(1,2)+delQ(1);
//when l=1,j=2,k=2,   pipe_con(abs(loop_con(1,3)),2)=pipe_con(2,2)=2;    Qv(1,2)=Qv(1,2)+delQ(pipe_con(abs(loop_con(1,3)),2))=Qv(1,2)-delQ(2);
//when l=1,j=2,k=3,   pipe_con(abs(loop_con(1,3)),3)=pipe_con(2,2)=0;--->stop
//when l=1,j=3,k=1,   pipe_con(abs(loop_con(1,4)),1)=pipe_con(3,1)=1;    Qv(1,3)=Qv(1,3)+delQ(pipe_con(abs(loop_con(1,4)),1))=Qv(1,3)+delQ(1);
//when l=1,j=3,k=2,   pipe_con(abs(loop_con(1,4)),2)=pipe_con(3,2)=0; --> stop
//when l=1,j=3,k=3,   pipe_con(abs(loop_con(1,4)),3)=pipe_con(3,3)=3;    Qv(1,3)=Qv(1,3)+delQ(pipe_con(abs(loop_con(1,4)),3))=Qv(1,3)-delQ(3);

////when l=2,j=1,k=1,   pipe_con(abs(loop_con(2,2)),1)=pipe_con(2,1)=1;    Qv(2,1)=Qv(2,1)+delQ(pipe_con(abs(loop_con(2,2)),1))=Qv(2,1)-delQ(1);
//Similarly other iterations will be there. I.T ---> (As 'l' stands for loop and 'j' stands for pipe, so from "loop_con" matrix, we can say that Qv(1,2)*,Qv(2,1)** both gives discharge of pipe 2. Similarly Qv(1,3),Qv(3,1) both gives discharge of pipe 3. Qv(2,2),Qv(2,3) both gives discharge of pipe 8. Those pairs should be added together to get the discharges of pipe no 2,3 & 8.)Wrong Idea,Discharge value of a pipe for any loop has same value. Because, for any loop it considers all effects.
// * In loop_con matrix, row 1,2 ,3 stands foe loop 1, loop2, loop3 respectively. So, Qv(1,2) means 2nd pipe of row 1 (i.e, pipe represented by (1,3)th element. As pipes are started from 2nd column, 1st column represents the total number of pipes connected to that particular loop).
// ** Similarly, Qv(2,1) means 1st pipe of row 2 (i.e, pipe represented by (2,2)th element.
for l=1:loopn
    rmse=rmse+delQ(l)^2;
end
    
rmse=sqrt(rmse/loopn);
count=count+1;
disp([count rmse])
end

//Print output
for l=1:loopn
    disp(['Loop Number:' string(l)])
    disp('Pipe Number  Discharge(m^3/s)')
    for i=1:loop_con(l,1)
        disp([abs(loop_con(l,i+1)) Qv(l,i)])
        end
    end















































































