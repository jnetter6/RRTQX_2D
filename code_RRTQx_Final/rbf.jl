# this file contains radial basis functions
# not using this file for Journal '20

function nu(t, Tf, m, n)
  Nu = zeros(Int((n+m)*(n+m+1)/2),Int((n+m)*(n+m+1)/2))
  for i = 1:Int((n+m)*(n+m+1)/2)
    for j = 1:Int((n+m)*(n+m+1)/2)
      if ((n+m)*(n+m+1)/2-i-j+1) >= 0
        Nu[i,j] = (t-Tf)^((n+m)*(n+m+1)/2-i-j+1)
      else
        Nu[i,j] = (t-Tf)^((n+m)*(n+m+1)-i-j+1)
      end
    end
  end
  Nu
end


function mu(t, Tf, n)
  Mu = zeros(n,n)
  for i = 1:n
    for j = 1:n
      if (n-i-j+1) >= 0
        Mu[i,j] = (t-Tf)^(n-i-j+1)
      else
        Mu[i,j] = (t-Tf)^(2n-i-j+1)
      end
    end
  end
  Mu
end
