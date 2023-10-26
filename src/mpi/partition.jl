fences(N,pieces) = [round(Int,i*N/pieces) for i in 0:pieces]
starts(N,pieces) = 1 .+ fences(N,pieces)[1:end-1]
ends(N,pieces) = fences(N,pieces)[2:end]
partitions(N,pieces) = [ s:e for (s,e) in zip(starts(N,pieces), ends(N,pieces))]

