function [val] = derivace_nn(State,Input,Noise,Time)
nh = 20;
ni = 9;
nt = length(State);
W1 = reshape(State(1:nh*ni),ni,nh)';
W2 = State((nh*ni+1):nt)';
y1 = [tanh(W1*Input);1];

val(nt-nh:nt) = y1';
for j = 1:nh,
  tmp = W2(j)*(1-y1(j).*y1(j));
  val((j-1)*ni+1:j*ni) = tmp*Input';
end

