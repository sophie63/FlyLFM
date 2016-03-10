% SEP goes once through the scrambled mixed speech signals, x 
% (which is of length P), in batch blocks of size B, adjusting weights,
% w, at the end of each block.
%
% I suggest a learning rate L, of 0.01 at least for 2->2 separation.
% But this will be unstable for higher dimensional data. Test it.
% Use smaller values. After convergence at a value for L, lower
% L and it will fine tune the solution.
%
% NOTE: this rule is the rule in our NC paper, but multiplied by w^T*w,
% as proposed by Amari, Cichocki & Yang at NIPS '95. This `natural
% gradient' method speeds convergence and avoids the matrix inverse in the 
% learning rule.

sweep=sweep+1; t=1;
noblocks=fix(P/B);
BI=B*Id;
for t=t:B:t-1+noblocks*B,
  u=w*x(:,t:t+B-1); 
  w=w+L*(BI+(1-2*(1./(1+exp(-u))))*u')*w;
end;
sepout