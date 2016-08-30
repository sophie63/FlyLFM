function A=byteestobit(B)

if B>=128
    A(1)=1;
    B=B-128;
else 
    A(1)=0;
end

if B>=64
    A(2)=1;
    B=B-64;
else 
    A(2)=0;
end

if B>=32
    A(3)=1;
    B=B-32;
else 
    A(3)=0;
end

if B>=16
    A(4)=1;
    B=B-16;
else 
    A(4)=0;
end

if B>=8
    A(5)=1;
    B=B-8;
else 
    A(5)=0;
end

if B>=4
    A(6)=1;
    B=B-4;
else 
    A(6)=0;
end

if B>=2
    A(7)=1;
    B=B-2;
else 
    A(7)=0;
end

if B>=1
    A(8)=1;
    B=B-1;
else 
    A(8)=0;
end

