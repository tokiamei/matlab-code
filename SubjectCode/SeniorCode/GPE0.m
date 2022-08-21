function [X,beta,phi] = fringes(N)
% calculate set of detection points x_1, x_2, ... x_N for 
% N atoms in 2 BEC clouds (N/2 in each) following Javanainen's PRL 1996 paper. 
% INPUT: 
% N - number of atoms in 2 BEC clouds (even integer)
% OUTPUT
% X - length N vector of detection points.
% beta, phi (lengh N vectors) - parameters of conditional probability
% distirbution for different numbers of detected atoms. 
%INITIAL SETUP
% N must be even
if mod(N,2)> 0
 error('N should be even')
end
% allocate arrays of detection coordinates X, joint probabilities jprob,
% betas and phis
X = zeros(N,1);
jprob = zeros(size(X));
beta = zeros(size(X));
phi = zeros(size(X));
% initial state vector of the system, an (N/2+1) x (N/2+1) matrix with 
% v_in{i,j} = A_{n_{+} = i-1,n_{-} = j-1} being amplitudes of the 
% number states |n_{+},n_{-}> (0 <= i,j, <= N/2). 
v_in = zeros(N/2+1,N/2+1); v_in(N/2+1,N/2+1) = 1;
% START THE RUN
m_already_detected = 0;
% generate the first detection coordinate
beta(1) = 0; phi(1) = 0;
X(1) = detect_atom(beta(1),phi(1));
% calculate state vector after the first detection
v_in = apply_field_operator(v_in,X(1),m_already_detected);
% find probability of detecting 1 atom at coordinate X(1)
jprob(1) = v_in(:)'*v_in(:); % we know it must be one
results(30) = 0; 
binVals(31) = 0;
binVals(1) = 0;
for i=2:31;
 binVals(i) = binVals(i-1)+ 1;
end
for i=1:30
 if (binVals(i) <= 30*X(1) & 30*X(1) < binVals(i+1));
 results(i) = results(i) + 1;
 found = i;
 break
 end
end
% continue
for m = 2:N
 [beta(m),phi(m)] = find_cond_prob(v_in,jprob(m-1),m - 1);
 % generate new detection coordinate
 X(m) = detect_atom(beta(m),phi(m));
 % calculate new state vector after the detection of the m-th atom 
 v_in = apply_field_operator(v_in,X(m),m-1);
 % find new joint probability
 jprob(m) = v_in(:)'*v_in(:);
 
 for i=2:31;
 binVals(i) = binVals(i-1)+ 1;
end
for i=1:30
 if (binVals(i) <= 30*X(m) & 30*X(m) < binVals(i+1));
 results(i) = results(i) + 1;
 found = i;
 break
 end
end
 
end
plot(results,'ok')
%--------------------------------------------------------------------------------------
function [beta,phi] = find_cond_prob(v_in,previous_joint_prob,m_already_detected)
%calculate parameters beta, phi of conditional probability
% p(x) = 1 + beta*cos(2*pi*x + phi)
 xs(1) = 0;
 v_tmp = apply_field_operator(v_in,xs(1),m_already_detected);
 ps(1) = v_tmp(:)'*v_tmp(:);
 xs(2) = 1/4;
 v_tmp = apply_field_operator(v_in,xs(2),m_already_detected);
 ps(2) = v_tmp(:)'*v_tmp(:);
 cond_probs = ps/previous_joint_prob;
 [beta,phi] = find_beta_phi(xs,cond_probs);
 function [beta,phi] = find_beta_phi(xs,ps)
% finds parameters of the distribution p(x) = 1 + beta cos(2pi x + phi)
% given two values of x and corresponding values of p(x)
% Find phi from the equation A*cos(phi) = B*sin(phi)
A = (ps(1) -1)*cos(2*pi*xs(2)) - (ps(2) - 1)*cos(2*pi*xs(1));
B = (ps(1) -1)*sin(2*pi*xs(2)) - (ps(2) - 1)*sin(2*pi*xs(1));
phi = atan2(A,B);
beta = (ps(1) - 1)/cos(2*pi*xs(1) + phi);
%--------------------------------------------------------------------------------------
function v_out = apply_field_operator(v_in,x,m_already_detected)
% v_out = apply_field_operator(v_in,x,m_already_detected)
% Applies field operator at the coordinate x to the state vector v_in
% of a system with m atoms already detected
% INPUT: 
% v_in - input state vector of the system (n x n matrix)
% m_already_detected - number of atoms already detected
% OUTPUT: 
% v_out - output state vector of the system (n x n matrix)
% Vector v_in is a square matrix with the entries (0 <= i,j <= N/2)
% v_in{i,j} = A_{n_{+} = i-1,n_{-} = j-1} being amplitudes of the 
% number states |n_{+},n_{-}>. 
n = size(v_in,1); % n = N/2 + 1
% form (n x n) matrix for the annihilation operator b_{+}
sqrt_line = sqrt(1:n)'; % column vector [sqrt(1),sqrt(2),... sqrt(n)]^{T};
b = repmat(sqrt_line,1,n);
% circularly shift v_in and zero the row with n_{+} = n
tmp = circshift(v_in,[-1,0]);
tmp(end,:) = 0;
% introduce phase multipilcand
phase_factor = exp(1i*pi*x);
% find result of action on v_in by operator exp(i*pi*x)*b_{+}
v_out = phase_factor*tmp.*b;
% find result of action on v_in by operator exp(-i*pi*x)*b_{-}
tmp = circshift(v_in,[0,-1]);
tmp(:,end) = 0;
v_out = v_out + conj(phase_factor)*tmp.*b';
% dividing by (N-m)^{1/2}
v_out = (2*n - 2 - m_already_detected)^(-1/2)*v_out;
%--------------------------------------------------------------------------------------
function x = detect_atom(beta,phi)
% generates random number x in the range 0 <= x <= 1 with the probability
% p(x) = 1 + beta*cos(2*pi*x + phi) using rejection method (see "Numerical Recipes"). 
% Comparison function f(x) = 1 + |beta|.
A = 1 + abs(beta); 
while 1
 x = rand(1);
 y = A*rand(1);
 if y < 1 + beta*cos(2*pi*x + phi)
 break
 end
end