function samples = generateHaltonSamples(dim,n,varargin)
  sf = [];

  if nargin == 3
    sf = varargin{1};
  end
  
  if (nargin == 2 || numel(sf) ~= dim)
    sf = ones(dim);
  end

  samples = zeros(dim, n);
  primeNums = primes(1000);
  for idx = 1:dim
    base_val = primeNums(idx);
    if sf(idx) > 0
      samples(idx,:) = sf(idx)*generateHaltonSequence(n,base_val);
    else
      samples(idx,:) = -2*sf(idx)*(generateHaltonSequence(n,base_val)-0.5);
    end
  end
end

function halton_sequence = generateHaltonSequence(n,base) 
  halton_sequence = zeros(n, 1);
  for idx = 1:n 
    halton_sequence(idx) = localHaltonSingleNumber(idx,base);
  end
end

function hn = localHaltonSingleNumber(n,b)
  n0 = n;
  hn = 0;
  f = 1/b;

  while(n0 > 0)
    n1 = floor(n0/b);
    r = n0 - b*n1;
    hn = hn + f*r;
    f = f/b;
    n0 = n1;
  end
end
