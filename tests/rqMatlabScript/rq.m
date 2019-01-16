A = [1 2 3 5; 5 12 1 0; 1 -1 1 -3]

[m n]=size(A);
if m>n
    error('RQ: Number of rows must be smaller than column');
end

[Q R]=qr(flipud(A).');
R=flipud(R.');
R(:,1:m)=R(:,m:-1:1);
Q=Q.';
Q(1:m,:)=Q(m:-1:1,:);

% Bruno



[R, Q] = rq(A);

